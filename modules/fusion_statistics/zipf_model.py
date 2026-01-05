#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from scipy.stats import linregress

parser = argparse.ArgumentParser(description="Zipf's Law modeling")
parser.add_argument("--input", required=True, help="Recurrence TSV")
parser.add_argument("--output", required=True, help="Zipf model TSV")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")

# Use rank and frequency
df = df[df["samples_detected"] > 0]

log_rank = np.log10(df["rank"])
log_freq = np.log10(df["samples_detected"])

slope, intercept, r_value, p_value, std_err = linregress(log_rank, log_freq)

df["zipf_predicted_log_freq"] = intercept + slope * log_rank
df["zipf_predicted_freq"] = 10 ** df["zipf_predicted_log_freq"]

df["zipf_alpha"] = -slope
df["zipf_r2"] = r_value ** 2
df["zipf_pvalue"] = p_value

df.to_csv(args.output, sep="\t", index=False)

print("Zipf Model Results:")
print(f"  Alpha (scaling exponent): {-slope:.4f}")
print(f"  R²: {r_value**2:.4f}")
print(f"  p-value: {p_value:.4e}")
