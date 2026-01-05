#!/usr/bin/env python3
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Fusion recurrence analysis")
parser.add_argument("--input", required=True, help="Merged fusion TSV")
parser.add_argument("--output", required=True, help="Recurrence output TSV")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")

# Create unique fusion key
df["fusion"] = df["geneA"] + "--" + df["geneB"]

# Count recurrence
recurrence = (
    df.groupby("fusion")
      .agg(
          geneA=("geneA", "first"),
          geneB=("geneB", "first"),
          samples_detected=("sample_id", "nunique")
      )
      .reset_index()
)

total_samples = df["sample_id"].nunique()
recurrence["recurrence_frequency"] = recurrence["samples_detected"] / total_samples

# Rank by frequency
recurrence = recurrence.sort_values(
    by="samples_detected", ascending=False
).reset_index(drop=True)

recurrence["rank"] = recurrence.index + 1

recurrence.to_csv(args.output, sep="\t", index=False)

print(f"[OK] Recurrence analysis written to {args.output}")
