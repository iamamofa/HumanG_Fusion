#!/usr/bin/env python3
import pandas as pd
import argparse

ONCOGENES = {
    "ALK", "BRAF", "EGFR", "RET", "NTRK1", "NTRK2", "NTRK3",
    "ROS1", "ABL1", "MET", "FGFR1", "FGFR2", "FGFR3"
}

parser = argparse.ArgumentParser(description="Clinical fusion annotation")
parser.add_argument("--input", required=True, help="Recurrence TSV")
parser.add_argument("--output", required=True, help="Clinically annotated TSV")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")

def clinical_flag(row):
    if row["geneA"] in ONCOGENES or row["geneB"] in ONCOGENES:
        return "oncogenic_candidate"
    return "unknown_significance"

df["clinical_significance"] = df.apply(clinical_flag, axis=1)

df.to_csv(args.output, sep="\t", index=False)

print("[OK] Clinical annotation complete")
