"""Command-line entry point."""

import argparse
import json
from pathlib import Path
import sys

from .build import build


def main(argv=None):
    parser = argparse.ArgumentParser(prog="sbasecouter", description="Build verified human spatial datasets")
    command = parser.add_subparsers(dest="command", required=True).add_parser("build")
    command.add_argument("--report", type=Path, required=True, action="append", help="SRAgent-spatial research batch.json")
    command.add_argument("--out", type=Path, required=True, help="New output directory")
    command.add_argument("--sample-id", action="append", default=[], help="Exact report sample ID. Repeat to select several")
    command.add_argument("--normalize-total", type=float, help="Also store a log1p library-size layer, for example 10000")
    args = parser.parse_args(argv)
    try:
        catalog = build(args.report,args.out,sample_ids=args.sample_id,normalize_total=args.normalize_total)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    totals = catalog["totals"]
    print(json.dumps(totals))
    return 0 if totals["built"] and not totals["failed"] else 1
