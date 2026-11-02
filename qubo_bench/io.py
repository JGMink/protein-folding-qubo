"""Versioned, line-oriented QUBO1 serialization shared with the C++ solver."""

import json
import math
from pathlib import Path
from .model import energy


def write_qubo(path, coefficients, n, offset, metadata=None):
    path = Path(path)
    if n < 1 or not isinstance(n, int):
        raise ValueError("n must be a positive integer")
    if not math.isfinite(offset):
        raise ValueError("offset must be finite")
    lines = ["QUBO1", f"n {n}", f"offset {float(offset):.17g}",
             "meta " + json.dumps(metadata or {}, sort_keys=True, separators=(",", ":"))]
    for (i, j), value in sorted(coefficients.items()):
        if not (0 <= i <= j < n) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("invalid QUBO coefficient")
        if value:
            lines.append(f"q {i} {j} {float(value):.17g}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("\n".join(lines) + "\n")


def read_qubo(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if len(lines) < 4 or lines[0] != "QUBO1":
        raise ValueError("expected QUBO1 header")
    if not lines[1].startswith("n ") or not lines[2].startswith("offset ") or not lines[3].startswith("meta "):
        raise ValueError("invalid QUBO1 header fields")
    n = int(lines[1][2:])
    offset = float(lines[2][7:])
    metadata = json.loads(lines[3][5:])
    if n < 1 or not math.isfinite(offset):
        raise ValueError("n and offset must be finite and valid")
    coefficients = {}
    for line in lines[4:]:
        fields = line.split()
        if len(fields) != 4 or fields[0] != "q":
            raise ValueError("invalid coefficient record")
        i, j, value = int(fields[1]), int(fields[2]), float(fields[3])
        if not 0 <= i <= j < n or (i, j) in coefficients or not math.isfinite(value):
            raise ValueError("invalid or duplicate coefficient index")
        coefficients[i, j] = value
    return coefficients, n, offset, metadata
