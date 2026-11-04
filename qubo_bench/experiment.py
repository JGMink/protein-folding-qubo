"""Run comparable serial trials and independently validate saved candidates."""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import time

from .io import read_qubo
from .model import energy, validate


def annotate(record, coefficients, n, offset, metadata):
    bits = record["best_bits"] if record["solver"] == "cpp" else record["final_bits"]
    if len(bits) != n or abs(energy(coefficients, offset, bits) - record["best_energy" if record["solver"] == "cpp" else "final_energy"]) > 1e-6:
        raise ValueError("reported energy disagrees with exported QUBO")
    valid, contact, path = validate(bits, metadata)
    record.update(valid=valid, contact_energy=contact, path=path)
    if record["solver"] == "cpp":
        final_bits = record["final_bits"]
        if abs(energy(coefficients, offset, final_bits) - record["final_energy"]) > 1e-6:
            raise ValueError("reported final energy disagrees with exported QUBO")
        final_valid, final_contact, final_path = validate(final_bits, metadata)
        record.update(final_valid=final_valid, final_contact_energy=final_contact, final_path=final_path)
    reference = metadata.get("reference")
    if reference:
        category, value = reference["category"], reference["energy"]
        record["reference_category"] = category
        record["reference_energy"] = value
        if "status" in reference:
            record["reference_status"] = reference["status"]
            if reference["status"] == "reported_exhaustively_established":
                observed = None if not valid else round(100 * contact)
                published = round(100 * value)
                record["matches_published_reference_energy"] = observed == published
                record["below_published_reference_investigate"] = (
                    observed is not None and observed < published)
                return record
        if category == "complete_qubo_global_minimum":
            record["reached_proven_qubo_minimum"] = record["best_energy" if record["solver"] == "cpp" else "final_energy"] <= value + 1e-6
        elif category == "feasible_fold_minimum":
            record["reached_proven_feasible_minimum"] = valid and contact <= value + 1e-6
        else:
            record["at_or_below_known_witness"] = valid and contact <= value + 1e-6
    return record


def run_cpp(args, instance, coefficients, n, offset, metadata):
    command = [str(args.solver.resolve()), str(instance), "--sweeps", str(args.sweeps),
               "--trials", str(args.trials), "--seed", str(args.seed), "--schedule", args.schedule,
               "--temp-start", str(args.temp_start), "--temp-end", str(args.temp_end),
               "--checkpoint", str(args.checkpoint)]
    if args.target is not None:
        command += ["--target", str(args.target)]
    if args.initial:
        command += ["--initial", args.initial]
    output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    records = []
    for line in output.splitlines():
        record = json.loads(line)
        record["solver"] = "cpp"
        record["instance"] = instance.name
        records.append(annotate(record, coefficients, n, offset, metadata))
    if len(records) != args.trials:
        raise ValueError("solver returned the wrong number of trials")
    return records


def run_neal(args, instance, coefficients, n, offset, metadata):
    if args.schedule != "geometric" or args.initial:
        raise ValueError("neal comparison supports geometric schedule and random initial states")
    import dimod
    try:
        import neal
        sampler = neal.SimulatedAnnealingSampler()
        package_name = "dwave-neal"
    except ImportError:
        from dwave.samplers import SimulatedAnnealingSampler
        sampler = SimulatedAnnealingSampler()
        package_name = "dwave-samplers"
    linear = {i: coefficients.get((i, i), 0.0) for i in range(n)}
    quadratic = {(i, j): value for (i, j), value in coefficients.items() if i != j}
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, vartype="BINARY")
    started = time.perf_counter()
    samples = sampler.sample(bqm, num_reads=args.trials, num_sweeps=args.sweeps,
                             beta_range=(1 / args.temp_start, 1 / args.temp_end), seed=args.seed,
                             beta_schedule_type="geometric")
    runtime = time.perf_counter() - started
    package_version = importlib.metadata.version(package_name)
    records = []
    for trial, datum in enumerate(samples.data(sorted_by=None)):
        bits = "".join(str(datum.sample[i]) for i in range(n))
        record = {"solver": "neal", "instance": instance.name, "trial": trial,
                  "sampler_package": package_name, "sampler_version": package_version,
                  "seed": args.seed, "schedule": "geometric",
                  "temperature_start": args.temp_start, "temperature_end": args.temp_end,
                  "sweeps": args.sweeps, "attempted_updates": args.sweeps * n,
                  "final_bits": bits, "final_energy": float(datum.energy),
                  "batch_runtime_seconds": runtime}
        records.append(annotate(record, coefficients, n, offset, metadata))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instances", nargs="+", type=Path)
    parser.add_argument("--solver", type=Path, help="standalone C++ solver executable")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sweeps", type=int, default=1000)
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--schedule", choices=("linear", "geometric", "fixed"), default="geometric")
    parser.add_argument("--temp-start", type=float, default=10)
    parser.add_argument("--temp-end", type=float, default=0.1)
    parser.add_argument("--checkpoint", type=int, default=100)
    parser.add_argument("--target", type=float)
    parser.add_argument("--initial")
    parser.add_argument("--neal", action="store_true", help="also run neal at the same sweeps/read budget")
    parser.add_argument("--neal-only", action="store_true", help="record a Python sampler baseline alone")
    args = parser.parse_args()
    if args.trials < 1 or args.sweeps < 0 or args.checkpoint < 0:
        parser.error("trials must be positive; sweeps and checkpoint must be nonnegative")
    if not args.neal_only and args.solver is None:
        parser.error("--solver is required unless --neal-only is selected")
    if args.neal and args.neal_only:
        parser.error("choose --neal or --neal-only")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as output:
        for instance in args.instances:
            coefficients, n, offset, metadata = read_qubo(instance)
            digest = hashlib.sha256(instance.read_bytes()).hexdigest()
            records = [] if args.neal_only else run_cpp(args, instance, coefficients, n, offset, metadata)
            if args.neal or args.neal_only:
                records += run_neal(args, instance, coefficients, n, offset, metadata)
            for record in records:
                record["instance_sha256"] = digest
                output.write(json.dumps(record, separators=(",", ":")) + "\n")
            for solver in sorted({record["solver"] for record in records}):
                group = [record for record in records if record["solver"] == solver]
                key = "best_energy" if solver == "cpp" else "final_energy"
                print(f"{instance.name} {solver}: {len(group)} trials, min {min(r[key] for r in group):.4f}, "
                      f"valid {sum(r['valid'] for r in group)}/{len(group)}")


if __name__ == "__main__":
    main()
