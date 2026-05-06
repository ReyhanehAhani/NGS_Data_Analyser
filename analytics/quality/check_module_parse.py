#!/usr/bin/env python3
"""Ensure shipped Python modules parse (lightweight CI for genomic tooling repos)."""
from __future__ import annotations

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    failures = []
    for pattern in ("filter_*.py",):
        for path in ROOT.glob(pattern):
            try:
                ast.parse(path.read_text(encoding="utf-8"))
            except SyntaxError as exc:
                failures.append(f"{path}: {exc}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("OK: filter_*.py files parse cleanly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
