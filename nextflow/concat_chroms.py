#!/usr/bin/env python3
"""
concat_chroms.py

Concatenate per-chromosome outputs from hic2structure into a single file.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate per-chromosome hic2structure outputs."
    )
    parser.add_argument(
        "--indir",
        required=True,
        type=Path,
        help="Directory containing per-chromosome output subdirectories."
    )
    parser.add_argument(
        "--chroms",
        required=True,
        type=Path,
        help="File with one chromosome ID per line (order preserved)."
    )
    parser.add_argument(
        "--resolution",
        required=True,
        type=int,
        help="Resolution used for hic2structure runs."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory where the final concatenated output should be written."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.indir.exists():
        sys.exit(f"[concat_chroms] ERROR: indir does not exist: {args.indir}")
    if not args.chroms.exists():
        sys.exit(f"[concat_chroms] ERROR: chroms file does not exist: {args.chroms}")

    chroms = [line.strip() for line in args.chroms.read_text().splitlines() if line.strip()]
    if not chroms:
        sys.exit("[concat_chroms] ERROR: chromosome list is empty")

    args.outdir.mkdir(parents=True, exist_ok=True)
    final_out = args.outdir / "structure.csv"

    print(f"[concat_chroms] Concatenating {len(chroms)} chromosomes")
    print(f"[concat_chroms] indir      : {args.indir}")
    print(f"[concat_chroms] resolution : {args.resolution}")
    print(f"[concat_chroms] outdir     : {args.outdir}")

    # TODO: implement actual concatenation logic here
    for s in chroms:
        struct_df = pd.read_csv(args.indir / s / "structure.csv")
        struct_df.loc[:,'chromosome'] = str(chrid)
        struct_df.loc[:,'id'] = ((struct_df.loc[:,'id']-1) * $args.resolution) + (args.resolution/2)
        structs.append(struct_df)
    comb_structs = pd.concat(structs)
    comb_structs.to_csv(final_out , index =  False)


if __name__ == "__main__":
    main()
