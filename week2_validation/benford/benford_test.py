"""
Chi-squared Goodness-of-Fit Test for Benford's Law
Week 2: Data Integrity & Statistical Validation (Preparation Only)

This module evaluates whether a dataset's first significant digit (FSD)
distribution conforms to Benford's Law using a chi-squared test.

No real data is used at this stage.
No external dependencies are required.

WHAT IS A CHI-SQUARED TEST?
A chi-squared test is a statistical method that compares what we actually
observed in our data versus what we expected to see. It answers the question:
"Is the difference between observed and expected just random chance, or is
there a real discrepancy?"

A small chi-squared value means the data closely matches expectations.
A large chi-squared value means there's a significant difference.
"""

# 'math' provides mathematical functions like logarithms and factorials
import math
# 'Counter' is a tool that counts how many times each item appears in a list
from collections import Counter

# Overflow guard for p-value series: cap term magnitude to avoid OverflowError
FLOAT_SAFE_LIMIT = 1e300


def chi_squared_statistic(observed, expected):
    """
    Compute the chi-squared statistic.
    
    This measures HOW DIFFERENT our observed data is from what Benford's Law
    predicts. The formula looks at each digit and asks: "How far off are we?"
    
    A SMALL result (close to 0) means the data matches Benford's Law well.
    A LARGE result means the data deviates significantly from Benford's Law.
    
    Args:
        observed: Dictionary of actual digit counts from our data.
                  Example: {1: 290, 2: 185, 3: 120, ...}
        expected: Dictionary of expected counts based on Benford's Law.
                  Example: {1: 301, 2: 176, 3: 125, ...}
    
    Returns:
        A single number representing the chi-squared statistic.
    """
    # Start with zero and add up the differences for each digit
    chi2 = 0.0
    
    # Check each digit from 1 to 9
    for digit in range(1, 10):
        # Get how many times this digit actually appeared (observed count)
        # If the digit isn't in our data, assume 0 occurrences
        o = observed.get(digit, 0)
        
        # Get how many times Benford's Law says this digit SHOULD appear
        e = expected.get(digit, 0)
        
        # Only calculate if we expected some occurrences
        # (to avoid dividing by zero)
        if e > 0:
            # The chi-squared formula for each category:
            # (observed - expected)² / expected
            # This gives more weight to larger differences
            chi2 += (o - e) ** 2 / e
    
    # Return the total chi-squared value
    return chi2


def chi_squared_p_value(chi2, degrees_of_freedom):
    """
    Calculate an approximate p-value for the chi-squared test.
    
    WHAT IS A P-VALUE?
    A p-value tells us the probability of seeing our results (or more extreme)
    if the data truly followed Benford's Law. Think of it as a "surprise score":
    
    - HIGH p-value (like 0.5): "Not surprising - data looks normal"
    - LOW p-value (like 0.01): "Very surprising - something unusual is happening"
    
    By convention, if p-value < 0.05, we say the data does NOT follow Benford's Law.
    
    Args:
        chi2: The chi-squared statistic we calculated.
        degrees_of_freedom: A technical parameter (always 8 for 9 digit categories).
    
    Returns:
        A probability between 0 and 1.
    
    NOTE: This is an approximation used for diagnostic purposes only,
    not for making scientific claims or formal statistical inference.
    """
    k = degrees_of_freedom / 2.0
    x = chi2 / 2.0

    # Series sum with overflow guard: if term exceeds FLOAT_SAFE_LIMIT, stop to avoid OverflowError
    total = 0.0
    for i in range(int(k)):
        try:
            term = (x ** i) / math.factorial(i)
        except (OverflowError, ValueError):
            break
        if not math.isfinite(term) or term > FLOAT_SAFE_LIMIT:
            break
        total += term
    return math.exp(-x) * total


def expected_benford_distribution(n):
    """
    Calculate the expected Benford distribution for a given sample size.
    
    This tells us how many times EACH first digit (1-9) SHOULD appear
    if the data perfectly follows Benford's Law.
    
    For example, with 1000 numbers:
        - Digit 1 should appear about 301 times (30.1%)
        - Digit 2 should appear about 176 times (17.6%)
        - Digit 9 should appear about 46 times (4.6%)
    
    Args:
        n: The total number of data points in our sample.
    
    Returns:
        A dictionary showing expected count for each digit 1-9.
    """
    # The Benford's Law formula: probability of digit d = log10(1 + 1/d)
    # Multiply by n to get expected COUNT instead of probability.
    # We calculate this for each digit from 1 to 9.
    return {d: n * math.log10(1 + 1 / d) for d in range(1, 10)}


def benford_chi_squared_test(fsd_counts, alpha=0.05):
    """
    Perform the complete chi-squared test for Benford's Law compliance.
    
    This is the main function that puts everything together. It takes the
    digit counts from your data and tells you whether the data appears to
    follow Benford's Law or not.
    
    Args:
        fsd_counts: Dictionary of first significant digit counts from data.
                    Example: {1: 300, 2: 180, 3: 125, 4: 95, ...}
        alpha: The significance threshold (default 0.05 means 5%).
               If p-value is below this, we conclude data doesn't follow Benford.
    
    Returns:
        A dictionary containing:
        - n: Total number of data points analyzed
        - chi_squared: The test statistic (higher = more different from Benford)
        - degrees_of_freedom: Technical parameter (always 8)
        - p_value: Probability value (lower = more likely NOT Benford)
        - decision: "Benford-consistent" or "Benford-inconsistent"
    """
    # STEP 1: Count total data points by adding up all digit counts
    n = sum(fsd_counts.values())
    
    # STEP 2: Calculate what Benford's Law predicts for this sample size
    expected = expected_benford_distribution(n)

    # STEP 3: Calculate how different our data is from the prediction
    chi2 = chi_squared_statistic(fsd_counts, expected)
    
    # STEP 4: Set degrees of freedom
    # With 9 categories (digits 1-9), degrees of freedom = 9 - 1 = 8
    # This is a technical requirement of the chi-squared test
    df = 8  # digits 1–9 → 9 categories → df = 8
    
    # STEP 5: Calculate the p-value (probability of seeing this result by chance)
    p_value = chi_squared_p_value(chi2, df)

    # STEP 6: Make a decision based on the p-value
    # If p-value is very low (< alpha), the data likely doesn't follow Benford
    # If p-value is higher (>= alpha), the data is consistent with Benford
    if p_value < alpha:
        decision = "Benford-inconsistent"
    else:
        decision = "Benford-consistent"

    # Return all the test results in a dictionary for easy access
    return {
        "n": n,                      # How many data points we tested
        "chi_squared": chi2,         # The test statistic value
        "degrees_of_freedom": df,    # Technical parameter
        "p_value": p_value,          # Probability measure
        "decision": decision,        # Final verdict
    }


# This section only runs when the file is executed directly (not imported).
# It demonstrates the chi-squared test using example data.
if __name__ == "__main__":
    # Create example data that roughly follows Benford's Law pattern:
    # - Digit 1 appears most often (300 times)
    # - Digit 9 appears least often (45 times)
    # This pattern mimics what Benford's Law predicts
    mock_fsd = Counter({1: 300, 2: 180, 3: 125, 4: 95, 5: 80, 6: 65, 7: 60, 8: 50, 9: 45})
    
    # Run the chi-squared test on this mock data
    result = benford_chi_squared_test(mock_fsd)

    # Display all the test results
    print("Chi-squared Test Result (Mock Data):")
    for key, value in result.items():
        print(f"{key}: {value}")
