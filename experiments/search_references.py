"""Bounded QUBO sampling and exhaustive lattice-path reference searches.

The path search optimizes only physically valid Hamiltonian paths. It is a
separate reference search, not a change to the exported QUBO formulation.
"""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time

from qubo_bench.io import read_qubo
from qubo_bench.model import energy, validate


def bits_for_path(path, residues, sites):
    bits = ["0"] * (residues * sites)
    for i, site in enumerate(path):
        bits[i * sites + site] = "1"
    return "".join(bits)


def record_candidate(bits, coefficients, offset, metadata):
    valid, contact, path = validate(bits, metadata)
    return {
        "bits": bits,
        "qubo_energy": energy(coefficients, offset, bits),
        "valid": valid,
        "contact_energy": contact,
        "path": path,
    }


def enumerate_paths(args, coefficients, n, offset, metadata):
    sequence = metadata["sequence"]
    sites = len(metadata["coordinates"])
    if len(sequence) != sites:
        raise ValueError("complete Hamiltonian-path enumeration requires one residue per site")
    neighbors = [[] for _ in range(sites)]
    for u, v in metadata["lattice_edges"]:
        neighbors[u].append(v)
        neighbors[v].append(u)
    indices = {aa: i for i, aa in enumerate(metadata["contact_alphabet"])}
    contact = metadata["contact_values"]
    sign = -1 if metadata["contact_model"] == "toy" else 1
    values = [[sign * contact[indices[a]][indices[b]] for b in sequence]
              for a in sequence]
    occupied = [-1] * sites
    path = []
    count = nodes = 0
    best = float("inf")
    best_path = None
    started = time.perf_counter()
    deadline = started + args.seconds

    def visit(depth, current_energy):
        nonlocal count, nodes, best, best_path
        nodes += 1
        if nodes % 100000 == 0 and time.perf_counter() >= deadline:
            raise TimeoutError
        if depth == sites:
            count += 1
            if current_energy < best - 1e-9:
                best, best_path = current_energy, tuple(path)
            return
        options = range(sites) if depth == 0 else neighbors[path[-1]]
        for site in options:
            if occupied[site] != -1:
                continue
            increment = sum(values[occupied[other]][depth]
                            for other in neighbors[site]
                            if occupied[other] != -1 and occupied[other] != depth - 1)
            occupied[site] = depth
            path.append(site)
            visit(depth + 1, current_energy + increment)
            path.pop()
            occupied[site] = -1

    complete = True
    try:
        visit(0, 0.0)
    except TimeoutError:
        complete = False
    result = {
        "method": "exhaustive_hamiltonian_paths",
        "seed": None,
        "instance": args.instance.as_posix(),
        "instance_sha256": hashlib.sha256(args.instance.read_bytes()).hexdigest(),
        "complete": complete,
        "oriented_paths_evaluated": count,
        "search_nodes": nodes,
        "runtime_seconds": time.perf_counter() - started,
        "time_limit_seconds": args.seconds,
        "best_path": best_path,
        "best_contact_energy": best if best_path is not None else None,
        "certainty": "proven_feasible_fold_minimum" if complete else "empirical_feasible_upper_bound",
    }
    if best_path is not None:
        result["best_candidate"] = record_candidate(
            bits_for_path(best_path, len(sequence), sites), coefficients, offset, metadata)
        if (not result["best_candidate"]["valid"]
                or abs(result["best_candidate"]["contact_energy"] - best) > 1e-6
                or abs(result["best_candidate"]["qubo_energy"] - best) > 1e-6):
            raise AssertionError("independent validator disagrees with path search")
    return result


def sample_qubo(args, coefficients, n, offset, metadata):
    import dimod
    from dwave.samplers import SimulatedAnnealingSampler, TabuSampler

    linear = {i: coefficients.get((i, i), 0.0) for i in range(n)}
    quadratic = {(i, j): value for (i, j), value in coefficients.items() if i != j}
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, "BINARY")
    settings = {"num_reads": args.reads, "seed": args.seed}
    if args.sampler == "sa":
        sampler = SimulatedAnnealingSampler()
        settings.update(num_sweeps=args.sweeps, beta_range=(1 / args.temp_start, 1 / args.temp_end),
                        beta_schedule_type="geometric")
    else:
        sampler = TabuSampler()
        settings.update(timeout=args.timeout_ms, tenure=args.tenure)
    started = time.perf_counter()
    samples = sampler.sample(bqm, **settings)
    elapsed = time.perf_counter() - started
    digest = hashlib.sha256(args.instance.read_bytes()).hexdigest()
    records = []
    for trial, datum in enumerate(samples.data(sorted_by=None)):
        bits = "".join(str(datum.sample[i]) for i in range(n))
        candidate = record_candidate(bits, coefficients, offset, metadata)
        if abs(candidate["qubo_energy"] - datum.energy) > 1e-6:
            raise AssertionError("sampler energy disagrees with exported QUBO")
        records.append({
            "method": "dwave_" + args.sampler,
            "instance": args.instance.as_posix(),
            "instance_sha256": digest,
            "trial": trial,
            "sampler_package": "dwave-samplers",
            "sampler_version": importlib.metadata.version("dwave-samplers"),
            "dimod_version": importlib.metadata.version("dimod"),
            "seed": args.seed,
            "settings": settings,
            "attempted_updates": args.sweeps * n if args.sampler == "sa" else None,
            "certainty": "empirical_qubo_sample",
            "batch_runtime_seconds": elapsed,
            **candidate,
        })
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    paths = subparsers.add_parser("paths", help="enumerate all Hamiltonian paths, subject to a time limit")
    paths.add_argument("instance", type=Path)
    paths.add_argument("--seconds", type=float, default=30)
    paths.add_argument("--output", type=Path, required=True)
    sample = subparsers.add_parser("sample", help="sample the unmodified QUBO")
    sample.add_argument("instance", type=Path)
    sample.add_argument("--sampler", choices=("sa", "tabu"), required=True)
    sample.add_argument("--reads", type=int, default=8)
    sample.add_argument("--seed", type=int, default=20261102)
    sample.add_argument("--sweeps", type=int, default=5000)
    sample.add_argument("--temp-start", type=float, default=10)
    sample.add_argument("--temp-end", type=float, default=0.1)
    sample.add_argument("--timeout-ms", type=int, default=100)
    sample.add_argument("--tenure", type=int, default=20)
    sample.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "paths" and args.seconds <= 0:
        parser.error("--seconds must be positive")
    if args.mode == "sample" and (args.reads < 1 or args.sweeps < 1 or args.timeout_ms < 1
                                  or args.tenure < 1 or args.temp_start <= 0 or args.temp_end <= 0):
        parser.error("sampling budgets, tenure, and temperatures must be positive")
    coefficients, n, offset, metadata = read_qubo(args.instance)
    result = enumerate_paths(args, coefficients, n, offset, metadata) if args.mode == "paths" else sample_qubo(
        args, coefficients, n, offset, metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as stream:
        if args.mode == "paths":
            stream.write(json.dumps(result) + "\n")
            print(f"paths: complete={result['complete']} count={result['oriented_paths_evaluated']} "
                  f"best={result['best_contact_energy']} seconds={result['runtime_seconds']:.2f}")
        else:
            for record in result:
                stream.write(json.dumps(record) + "\n")
            print(f"{args.sampler}: reads={len(result)} best_qubo={min(r['qubo_energy'] for r in result):.2f} "
                  f"valid={sum(r['valid'] for r in result)}")


if __name__ == "__main__":
    main()
