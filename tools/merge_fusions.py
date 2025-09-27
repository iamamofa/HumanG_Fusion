#!/usr/bin/env python3
"""
Merge fusion caller outputs (Arriba, STAR-Fusion, FusionCatcher).
Creates a unified TSV report with caller support flags.
"""

import argparse
import pandas as pd

def load_table(path, caller):
    if not path or not path.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()
    df["caller"] = caller
    return df

def main():
    parser = argparse.ArgumentParser(description="Merge fusion caller outputs")
    parser.add_argument("--arriba", help="Arriba results TSV")
    parser.add_argument("--starfusion", help="STAR-Fusion results TSV")
    parser.add_argument("--fusioncatcher", help="FusionCatcher results TXT/TSV")
    parser.add_argument("--gtf", help="Optional: GTF annotation", default=None)
    parser.add_argument("--out", required=True, help="Output TSV")

    args = parser.parse_args()
    frames = []

    for caller, path in [("arriba", args.arriba),
                         ("starfusion", args.starfusion),
                         ("fusioncatcher", args.fusioncatcher)]:
        if path:
            frames.append(load_table(pd.Path(path), caller))

    if not frames:
        open(args.out, "w").write("caller\tfusion\n")
        return

    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    main()
