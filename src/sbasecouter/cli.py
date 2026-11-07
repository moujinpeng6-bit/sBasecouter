"""Command-line entry point."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

from .build import build
from .context import build_context
from .recount import recount


def main(argv=None):
    parser = argparse.ArgumentParser(prog="sbasecouter", description="Build verified human spatial datasets")
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser("build")
    command.add_argument("--report", type=Path, required=True, action="append", help="SRAgent-spatial research batch.json")
    command.add_argument("--out", type=Path, required=True, help="New output directory")
    command.add_argument("--sample-id", action="append", default=[], help="Exact report sample ID. Repeat to select several")
    command.add_argument("--normalize-total", type=float, help="Also store a log1p library-size layer, for example 10000")
    context = subparsers.add_parser("context", help="Add a cell neighborhood and derived expression layers")
    context.add_argument("--input", type=Path, required=True)
    context.add_argument("--out", type=Path, required=True)
    context.add_argument("--cell-types", type=Path, help="Optional same-slice slice_id,barcode,cell_type CSV")
    context.add_argument("--cell-types-source", help="Reference for the supplied annotation table, saved in the output")
    context.add_argument("--neighbors", type=int, default=8)
    context.add_argument("--radius-um", type=float, default=50)
    context.add_argument("--alpha", type=float, default=0.25, help="Neighbor share for cells with measured RNA")
    context.add_argument("--same-type-weight", type=float, default=2.0)
    raw = subparsers.add_parser("recount", help="Run a vendor spatial recount from a sample manifest")
    raw.add_argument("--manifest", type=Path, required=True)
    raw.add_argument("--out", type=Path, required=True)
    raw.add_argument("--engine", help="Path to the pinned Space Ranger or Xenium Ranger executable")
    raw.add_argument("--cores", type=int)
    raw.add_argument("--memory-gb", type=int)
    raw.add_argument("--plan", action="store_true", help="Check the declared inputs and print the command")
    args = parser.parse_args(argv)
    try:
        if args.command == "recount":
            result = recount(args.manifest, args.out, executable=args.engine, cores=args.cores,
                             memory_gb=args.memory_gb, plan_only=args.plan)
            print(json.dumps(result))
            return 0
        if args.command == "context":
            result = build_context(args.input, args.out, neighbors=args.neighbors, radius=args.radius_um,
                                   alpha=args.alpha, same_type_weight=args.same_type_weight, label_file=args.cell_types,
                                   label_source=args.cell_types_source)
            print(json.dumps(result))
            return 0
        catalog = build(args.report,args.out,sample_ids=args.sample_id,normalize_total=args.normalize_total)
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    totals = catalog["totals"]
    print(json.dumps(totals))
    return 0 if totals["built"] and not totals["failed"] else 1
