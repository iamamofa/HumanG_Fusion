#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from scipy.stats import chisquare

parser = argparse.ArgumentParser(description="Benford's Law test")
parser.add_argument("--input", required=True, help="Recurrence TSV")
parser.add_argument("--output", required=True, help="Benford test TSV")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")

# Extract first digit
df = df[df["samples_detected"] > 0]
df["first_digit"] = df["samples_detected"].astype(str).str[0].astype(int)

observed = df["first_digit"].value_counts().sort_index()
digits = np.arange(1, 10)

# Expected Benford distribution
expected = np.log10(1 + 1 / digits)
expected = expected * observed.sum()

chi2, p = chisquare(
    f_obs=observed.reindex(digits, fill_value=0),
    f_exp=expected
)

benford_df = pd.DataFrame({
    "digit": digits,
    "observed": observed.reindex(digits, fill_value=0).values,
    "expected": expected
})

benford_df["chi2_stat"] = chi2
benford_df["p_value"] = p

benford_df.to_csv(args.output, sep="\t", index=False)

print(f"[OK] Benford test completed (χ²={chi2:.3f}, p={p:.4e})")
