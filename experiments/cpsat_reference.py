"""Time-limited CP-SAT search for a full-lattice protein-fold reference.

This model searches physical Hamiltonian paths directly. It does not modify or
optimize the exported QUBO formulation. Requires ortools (tested with 9.15.6755).
Run from the repository root with ``python -m experiments.cpsat_reference``.
"""

import argparse
from decimal import Decimal
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import time

from qubo_bench.io import read_qubo
from qubo_bench.model import energy, validate


def encoded_candidate(path, coefficients, offset, metadata):
    sites = len(metadata["coordinates"])
    bits = ["0"] * (len(path) * sites)
    for residue, site in enumerate(path):
        bits[residue * sites + site] = "1"
    bits = "".join(bits)
    valid, contact, checked_path = validate(bits, metadata)
    qubo = energy(coefficients, offset, bits)
    if not valid or checked_path != list(path) or abs(contact - qubo) > 1e-6:
        raise AssertionError("CP-SAT candidate failed independent feasibility or energy check")
    return {"path": list(path), "bits": bits, "validated_contact_energy": contact,
            "validated_qubo_energy": qubo, "valid": valid}


def integer_contact_table(metadata):
    """Exact integer cents for each residue pair; chain bonds contribute zero."""
    sequence = metadata["sequence"]
    indices = {aa: i for i, aa in enumerate(metadata["contact_alphabet"])}
    values = metadata["contact_values"]
    sign = -1 if metadata["contact_model"] == "toy" else 1
    table = []
    for i, aa in enumerate(sequence):
        for j, bb in enumerate(sequence):
            scaled = Decimal(str(values[indices[aa]][indices[bb]])) * 100 * sign
            if scaled != scaled.to_integral_value():
                raise ValueError("contact coefficient is not an exact hundredth")
            table.append(0 if abs(i - j) <= 1 else int(scaled))
    return table


def solve(args):
    from ortools.sat.python import cp_model

    coefficients, n, offset, metadata = read_qubo(args.instance)
    sequence = metadata["sequence"]
    sites = len(metadata["coordinates"])
    if sites != len(sequence) or n != sites * sites or len(metadata["lattice_dimensions"]) != 2:
        raise ValueError("this reference search requires a fully occupied 2D lattice")
    edges = [tuple(edge) for edge in metadata["lattice_edges"]]
    directed_edges = edges + [(v, u) for u, v in edges]
    contact_table = integer_contact_table(metadata)
    lo, hi = min(contact_table), max(contact_table)
    reference = metadata.get("reference", {})
    reference_path = reference.get("path")
    reference_candidate = (encoded_candidate(reference_path, coefficients, offset, metadata)
                           if reference_path is not None else None)

    started = time.perf_counter()
    model = cp_model.CpModel()
    # position[i] is the site occupied by residue i; occupant[u] is the residue
    # at site u. AddInverse makes them permutations, enforcing non-overlap.
    position = [model.new_int_var(0, sites - 1, f"position_{i}") for i in range(sites)]
    occupant = [model.new_int_var(0, sites - 1, f"occupant_{u}") for u in range(sites)]
    model.add_inverse(position, occupant)
    for i in range(sites - 1):
        model.add_allowed_assignments([position[i], position[i + 1]], directed_edges)

    # Every occupied lattice edge is either a chain bond (cost zero) or a
    # nonbonded contact. The pair table uses the existing signed coefficients.
    edge_costs = []
    for edge_index, (u, v) in enumerate(edges):
        cost = model.new_int_var(lo, hi, f"edge_cost_{edge_index}")
        model.add_element(occupant[u] * sites + occupant[v], contact_table, cost)
        edge_costs.append(cost)
    model.minimize(sum(edge_costs))

    if reference_path is not None:
        inverse = [0] * sites
        for i, u in enumerate(reference_path):
            inverse[u] = i
            model.add_hint(position[i], u)
        for u, i in enumerate(inverse):
            model.add_hint(occupant[u], i)
        for edge_index, (u, v) in enumerate(edges):
            model.add_hint(edge_costs[edge_index], contact_table[inverse[u] * sites + inverse[v]])
    build_seconds = time.perf_counter() - started

    class ValidatingIncumbents(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.count = 0
            self.history = []
            self.best = None

        def on_solution_callback(self):
            path = [self.value(var) for var in position]
            candidate = encoded_candidate(path, coefficients, offset, metadata)
            if abs(candidate["validated_contact_energy"] * 100 - self.objective_value) > 1e-5:
                raise AssertionError("CP-SAT objective disagrees with independent validator")
            self.count += 1
            self.best = candidate
            self.history.append({"energy": candidate["validated_contact_energy"],
                                 "solver_seconds": self.wall_time})

    incumbents = ValidatingIncumbents()
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = args.workers
    solver.parameters.random_seed = args.seed
    solve_started = time.perf_counter()
    status = solver.solve(model, incumbents)
    solve_seconds = time.perf_counter() - solve_started
    status_name = solver.status_name(status)
    candidate = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        candidate = encoded_candidate([solver.value(var) for var in position],
                                      coefficients, offset, metadata)
        if abs(candidate["validated_contact_energy"] * 100 - solver.objective_value) > 1e-5:
            raise AssertionError("final CP-SAT objective disagrees with validator")
    bound = solver.best_objective_bound
    bound = bound / 100 if math.isfinite(bound) else None
    result = {
        "method": "lattice_aware_cp_sat",
        "instance": args.instance.as_posix(),
        "instance_sha256": hashlib.sha256(args.instance.read_bytes()).hexdigest(),
        "ortools_version": importlib.metadata.version("ortools"),
        "seed": args.seed,
        "workers": args.workers,
        "time_limit_seconds": args.seconds,
        "model_build_seconds": build_seconds,
        "solver_runtime_seconds": solve_seconds,
        "solver_reported_wall_time_seconds": solver.wall_time,
        "solver_status": status_name,
        "solutions_validated": incumbents.count,
        "incumbent_history": incumbents.history,
        "best_objective_bound": bound,
        "candidate": candidate,
        "reference_candidate": reference_candidate,
        "certainty": ("proven_feasible_fold_minimum" if status == cp_model.OPTIMAL
                      else "time_limited_feasible_incumbent" if candidate else "no_solver_incumbent"),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instance", type=Path)
    parser.add_argument("--seconds", type=float, default=20)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20261102)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.seconds <= 0 or args.workers < 1 or args.seed < 0:
        parser.error("seconds and workers must be positive; seed must be nonnegative")
    result = solve(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    candidate = result["candidate"]
    print(f"status={result['solver_status']} candidate={candidate['validated_contact_energy'] if candidate else None} "
          f"bound={result['best_objective_bound']} seconds={result['solver_runtime_seconds']:.2f} "
          f"validated={result['solutions_validated']}")


if __name__ == "__main__":
    main()
