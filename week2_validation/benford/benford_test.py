"""
Chi-squared Goodness-of-Fit Test for Benford's Law
Week 2: Data Integrity & Statistical Validation (Preparation Only)

This module evaluates whether a dataset's first significant digit (FSD)
distribution conforms to Benford's Law using a chi-squared test.

No real data is used at this stage.
No external dependencies are required.
"""

import math
from collections import Counter


def chi_squared_statistic(observed, expected):
    """
    Compute chi-squared statistic.
    """
    chi2 = 0.0
    for digit in range(1, 10):
        o = observed.get(digit, 0)
        e = expected.get(digit, 0)
        if e > 0:
            chi2 += (o - e) ** 2 / e
    return chi2


def chi_squared_p_value(chi2, degrees_of_freedom):
    """
    Approximate chi-squared p-value using survival function approximation.
    Avoids external libraries.
    """
    # Using incomplete gamma approximation
    k = degrees_of_freedom / 2.0
    x = chi2 / 2.0
    return math.exp(-x) * sum((x ** i) / math.factorial(i) for i in range(int(k)))


def expected_benford_distribution(n):
    """
    Expected Benford distribution for sample size n.
    """
    return {d: n * math.log10(1 + 1 / d) for d in range(1, 10)}


def benford_chi_squared_test(fsd_counts, alpha=0.05):
    """
    Perform chi-squared test and return statistic, p-value, and decision.
    """
    n = sum(fsd_counts.values())
    expected = expected_benford_distribution(n)

    chi2 = chi_squared_statistic(fsd_counts, expected)
    df = 8  # digits 1–9 → 9 categories → df = 8
    p_value = chi_squared_p_value(chi2, df)

    if p_value < alpha:
        decision = "Benford-inconsistent"
    else:
        decision = "Benford-consistent"

    return {
        "n": n,
        "chi_squared": chi2,
        "degrees_of_freedom": df,
        "p_value": p_value,
        "decision": decision,
    }


if __name__ == "__main__":
    # Example test using mock data (prep only)
    mock_fsd = Counter({1: 300, 2: 180, 3: 125, 4: 95, 5: 80, 6: 65, 7: 60, 8: 50, 9: 45})
    result = benford_chi_squared_test(mock_fsd)

    print("Chi-squared Test Result (Mock Data):")
    for key, value in result.items():
        print(f"{key}: {value}")
