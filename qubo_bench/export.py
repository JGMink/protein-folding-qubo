"""Export initial benchmarks or any sequence/lattice supported by the model."""

import argparse
from pathlib import Path
from .io import write_qubo
from .model import build, validate
from .mj85_benchmarks import BENCHMARKS, definition


PRESETS = {
    "toy16": ("HPCH", (2, 2), "row-major", "toy", (1, 1, 1),
              "complete_qubo_global_minimum", -1.0, (0, 1, 3, 2)),
    "fold64": ("CFIKDLWV", (2, 4), "row-major", "signed20", (3, 4, 4),
               "feasible_fold_minimum", -2.84, (0, 1, 2, 3, 7, 6, 5, 4)),
    "fold576": ("ENDVTMDMDPSPCLFRIHNLPRAH", (4, 6), "row-major", "signed20", (3, 4, 4),
                "feasible_fold_minimum", -10.05,
                (1, 0, 6, 12, 18, 19, 20, 21, 22, 23, 17, 16, 15, 14, 13, 7, 8, 2, 3, 9, 10, 11, 5, 4)),
    "witness1024": ("ENDVTMDMDPSPCLFRIHNLPRAHSFDRFGWH", (8, 4), "column-major",
                    "signed20", (1.5, 2, 2), "known_valid_configuration", -14.25,
                    tuple((site % 4) * 8 + site // 4 for site in
                          (9, 10, 14, 18, 22, 21, 25, 26, 30, 29, 28, 24, 20, 16, 17, 13,
                           12, 8, 4, 0, 1, 5, 6, 2, 3, 7, 11, 15, 19, 23, 27, 31))),
}


def export(path, sequence, dims, order, model, penalties, reference=None, benchmark=None):
    coefficients, n, offset, metadata = build(sequence, dims, penalties, model, order)
    if benchmark is not None:
        metadata["benchmark_id"] = benchmark["id"]
        metadata["topology"] = benchmark["topology"]
        metadata["full_occupancy"] = benchmark["full_occupancy"]
        metadata["reference"] = benchmark["reference"]
    if reference is not None:
        category, value, path_witness = reference
        metadata["reference"] = {"category": category, "energy": value}
        if path_witness is not None:
            metadata["reference"]["path"] = path_witness
            bits = ["0"] * n
            for i, site in enumerate(path_witness):
                bits[i * len(metadata["coordinates"]) + site] = "1"
            valid, direct, _ = validate("".join(bits), metadata)
            if not valid or abs(direct - value) > 1e-8:
                raise ValueError("reference witness disagrees with exported model")
    write_qubo(path, coefficients, n, offset, metadata)
    return n, len(coefficients)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--all", action="store_true", help="export the initial benchmark collection")
    choice.add_argument("--preset", choices=tuple(PRESETS) + tuple(BENCHMARKS))
    choice.add_argument("--sequence", help="export an additional protein sequence")
    parser.add_argument("--dims", nargs="+", type=int, help="2D or 3D lattice dimensions")
    parser.add_argument("--order", choices=("row-major", "column-major", "x-fastest"))
    parser.add_argument("--model", choices=("signed20", "toy", "MJ85_TableVI"), default="signed20")
    parser.add_argument("--penalties", nargs=3, type=float, default=(3, 4, 4))
    parser.add_argument("--output", type=Path, required=True, help="file, or directory with --all")
    args = parser.parse_args()
    if args.all:
        args.output.mkdir(parents=True, exist_ok=True)
        for name, (sequence, dims, order, model, penalties, category, value, path) in PRESETS.items():
            target = args.output / f"{name}.qubo"
            n, terms = export(target, sequence, dims, order, model, penalties,
                              (category, value, path))
            print(f"{target}: {n} variables, {terms} terms")
    elif args.preset:
        if args.preset in BENCHMARKS:
            benchmark = definition(args.preset)
            n, terms = export(args.output, benchmark["sequence"], benchmark["lattice_dimensions"],
                              benchmark["site_order"], benchmark["contact_model"], args.penalties,
                              benchmark=benchmark)
        else:
            sequence, dims, order, model, penalties, category, value, path = PRESETS[args.preset]
            n, terms = export(args.output, sequence, dims, order, model, penalties,
                              (category, value, path))
        print(f"{args.output}: {n} variables, {terms} terms")
    else:
        if not args.dims or len(args.dims) not in (2, 3):
            parser.error("--sequence requires --dims with two or three integers")
        order = args.order or ("x-fastest" if len(args.dims) == 3 else "row-major")
        n, terms = export(args.output, args.sequence, args.dims, order, args.model, args.penalties)
        print(f"{args.output}: {n} variables, {terms} terms")


if __name__ == "__main__":
    main()
