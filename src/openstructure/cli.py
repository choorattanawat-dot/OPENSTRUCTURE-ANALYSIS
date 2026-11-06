"""Command-line entry point.

``openstructure run model.yaml --out results`` reads a model, runs the linear
analysis, writes results.json, report.md, and force-diagram figures.
"""

from __future__ import annotations

import argparse
import sys

from openstructure import __version__
from openstructure.errors import OpenStructureError


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="openstructure",
        description="Linear-static 2D frame analysis and RC member design (v0.1).",
    )
    parser.add_argument("--version", action="version", version=f"openstructure {__version__}")
    sub = parser.add_subparsers(dest="command")
    run = sub.add_parser("run", help="Analyse a YAML or JSON model and write the report.")
    run.add_argument("model", help="Path to the model file (.yaml, .yml, or .json).")
    run.add_argument(
        "--out",
        default="results",
        help="Output directory for results.json, report.md, and diagrams (default: results).",
    )
    run.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip the matplotlib force diagrams.",
    )
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if args.command == "run":
        from openstructure.pipeline import run_file

        try:
            output = run_file(args.model, args.out, plots=not args.no_plots)
        except OpenStructureError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        print(f"Wrote {args.out}/results.json and {args.out}/report.md")
        for packet in output.designs:
            print(f"  {packet.kind} [{packet.code_id}] {packet.check_id}: {packet.status} — {packet.summary}")
        return 0
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
