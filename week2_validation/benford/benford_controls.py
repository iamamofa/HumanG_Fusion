"""
Benford Control Data Generators
Week 2: Data Integrity & Statistical Validation (Preparation Only)

This script generates:
1. Benford-compliant synthetic data (positive control)
2. Non-Benford (uniform) synthetic data (negative control)

No real data is used.
No external connections are made.
"""

import math
import random
from collections import Counter


def first_significant_digit(value):
    """
    Extract the first significant digit from a positive number.
    """
    value = abs(value)
    if value == 0:
        return None
    return int(str(value).lstrip("0.")[0])


def generate_benford_compliant_data(n):
    """
    Generate Benford-compliant data using inverse transform sampling.
    """
    data = []
    for _ in range(n):
        u = random.random()
        value = 10 ** u
        data.append(value)
    return data


def generate_uniform_data(n):
    """
    Generate uniform random data that does NOT follow Benford's Law.
    """
    return [random.uniform(1, 10) for _ in range(n)]


def compute_fsd_distribution(data):
    """
    Compute first significant digit distribution.
    """
    fsd_counts = Counter()
    for x in data:
        d = first_significant_digit(x)
        if d is not None:
            fsd_counts[d] += 1
    return fsd_counts


def expected_benford_distribution(n):
    """
    Expected Benford distribution for sample size n.
    """
    return {d: n * math.log10(1 + 1 / d) for d in range(1, 10)}


if __name__ == "__main__":
    N = 1000  # control size (safe, arbitrary)

    benford_data = generate_benford_compliant_data(N)
    uniform_data = generate_uniform_data(N)

    benford_fsd = compute_fsd_distribution(benford_data)
    uniform_fsd = compute_fsd_distribution(uniform_data)

    expected = expected_benford_distribution(N)

    print("Positive Control (Benford-compliant):")
    print(benford_fsd)
    print("\nNegative Control (Uniform):")
    print(uniform_fsd)
    print("\nExpected Benford Distribution:")
    print(expected)
