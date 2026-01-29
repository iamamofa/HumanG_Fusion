"""
Benford Control Data Generators
Week 2: Data Integrity & Statistical Validation (Preparation Only)

This script generates:
1. Benford-compliant synthetic data (positive control)
2. Non-Benford (uniform) synthetic data (negative control)

No real data is used.
No external connections are made.

WHAT IS BENFORD'S LAW?
Benford's Law is a mathematical observation that in many real-world datasets
(like financial records, scientific data, population numbers), the first digit
of numbers is not evenly distributed. Instead:
- The digit 1 appears as the first digit about 30% of the time
- The digit 2 appears about 18% of the time
- The digit 9 appears only about 5% of the time

This pattern is used to detect data anomalies or potential fraud.
"""

# 'math' provides mathematical functions like logarithms
import math
# 'random' allows us to generate random numbers for creating test data
import random
# 'Counter' is a tool that counts how many times each item appears in a list
from collections import Counter


def first_significant_digit(value):
    """
    Extract the first significant digit from a positive number.
    
    The "first significant digit" is the leftmost non-zero digit in a number.
    For example:
        - 1234 has first significant digit 1
        - 0.0045 has first significant digit 4
        - 987 has first significant digit 9
    """
    # Make the number positive (remove any negative sign)
    # This ensures we handle negative numbers correctly
    value = abs(value)
    
    # Zero has no significant digit because there's no non-zero digit
    # We return None to indicate "no valid digit found"
    if value == 0:
        return None
    
    # Convert the number to text, remove leading zeros and decimal points,
    # then take the first character. Finally, convert it back to a number.
    # Example: 0.0045 becomes "0.0045" -> "45" -> "4" -> 4
    return int(str(value).lstrip("0.")[0])


def generate_benford_compliant_data(n):
    """
    Generate Benford-compliant data using inverse transform sampling.
    
    This creates fake (synthetic) numbers that follow Benford's Law pattern.
    We use this as a "positive control" - data we KNOW should pass the
    Benford test - to verify our testing code works correctly.
    
    Think of it like testing a smoke detector with actual smoke to make
    sure it beeps when it should.
    """
    # Create an empty list to store our generated numbers
    data = []
    
    # Generate 'n' numbers (n is how many numbers we want)
    for _ in range(n):
        # Generate a random decimal between 0 and 1 (like 0.347 or 0.892)
        u = random.random()
        
        # This mathematical transformation (10 raised to the power of u)
        # creates numbers that naturally follow Benford's Law distribution.
        # The math behind this: when you take 10^u where u is between 0 and 1,
        # you get numbers between 1 and 10 with the right Benford pattern.
        value = 10 ** u
        
        # Add this number to our list
        data.append(value)
    
    # Return the complete list of generated numbers
    return data

# PRIMARY POSITIVE CONTROL:
# This function creates test data that we KNOW follows Benford's Law.
# We use it to verify our Benford analysis code is working correctly.
# This is never used on real patient or research data.


def generate_benford_compliant_data_multidecade(n):
    """
    Generate synthetic Benford-compliant data spanning different magnitudes.
    
    This creates fake numbers of various sizes (from small like 1.5 to large
    like 150,000) that all follow Benford's Law. Having numbers of different
    sizes makes the test data more realistic.
    
    POSITIVE CONTROL ONLY - used to verify our testing code works.
    
    Args:
        n: How many numbers to generate.
    
    Returns:
        A list of 'n' numbers that follow Benford's Law.
    """
    # Create an empty list to store our generated numbers
    result = []
    
    # Generate 'n' numbers
    for _ in range(n):
        # Generate a random decimal between 0 and 1
        u = random.random()
        
        # Create a number between 1 and 10 that follows Benford's pattern.
        # This is called the "significand" - the meaningful digits of a number.
        # Example: in 4,500 the significand would be 4.5
        significand = 10.0 ** u
        
        # Pick a random "magnitude" (how big the number should be).
        # 0 means ones (1-10), 1 means tens (10-100), 2 means hundreds, etc.
        # We go up to 5, which means numbers up to 100,000s.
        exponent = random.randint(0, 5)
        
        # Combine the significand with the magnitude to get the final number.
        # Example: significand 4.5 with exponent 2 gives 4.5 × 100 = 450
        value = significand * (10.0 ** exponent)
        
        # Add this number to our results
        result.append(value)
    
    # Return all the generated numbers
    return result

def generate_uniform_data(n):
    """
    Generate uniform random data that does NOT follow Benford's Law.
    
    This creates numbers where each first digit (1-9) appears equally often,
    which is NOT what Benford's Law predicts. We use this as a "negative control" -
    data we KNOW should FAIL the Benford test - to verify our testing code
    correctly identifies non-Benford data.
    
    Think of it like testing a smoke detector with fresh air to make sure
    it stays quiet when there's no smoke.
    """
    # Generate 'n' random numbers evenly spread between 1 and 10.
    # In this range, each first digit (1-9) has roughly equal chance,
    # which violates Benford's Law (where 1 should appear ~30% of the time).
    return [random.uniform(1, 10) for _ in range(n)]


def compute_fsd_distribution(data):
    """
    Compute first significant digit distribution.
    
    This function counts how many times each digit (1-9) appears as the
    first significant digit in our dataset. The result shows the pattern
    of first digits, which we can then compare to Benford's Law.
    
    Example output: {1: 301, 2: 175, 3: 125, ...} means digit 1 appeared
    301 times as the first digit, digit 2 appeared 175 times, etc.
    """
    # Create a counting tool to track how many times each digit appears
    fsd_counts = Counter()
    
    # Go through each number in our data
    for x in data:
        # Get the first significant digit of this number
        d = first_significant_digit(x)
        
        # Only count valid digits (skip zeros which return None)
        if d is not None:
            # Add 1 to the count for this digit
            fsd_counts[d] += 1
    
    # Return the final counts for each digit
    return fsd_counts


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
    # We do this for each digit from 1 to 9.
    return {d: n * math.log10(1 + 1 / d) for d in range(1, 10)}


# This section only runs when the file is executed directly (not imported).
# It demonstrates the control data generators and shows their distributions.
if __name__ == "__main__":
    # Number of test data points to generate.
    # 1000 is enough to see the Benford pattern clearly.
    N = 1000  # control size (safe, arbitrary)

    # STEP 1: Generate two types of test data
    # Positive control: data that SHOULD follow Benford's Law
    benford_data = generate_benford_compliant_data(N)
    # Negative control: data that should NOT follow Benford's Law
    uniform_data = generate_uniform_data(N)

    # STEP 2: Count the first digits in each dataset
    # This shows us the actual distribution of first digits
    benford_fsd = compute_fsd_distribution(benford_data)
    uniform_fsd = compute_fsd_distribution(uniform_data)

    # STEP 3: Calculate what Benford's Law predicts
    # This is what "perfect" Benford data would look like
    expected = expected_benford_distribution(N)

    # STEP 4: Display all results for comparison
    # The positive control should closely match the expected distribution
    print("Positive Control (Benford-compliant):")
    print(benford_fsd)
    # The negative control should look very different from expected
    print("\nNegative Control (Uniform):")
    print(uniform_fsd)
    # The theoretical expectation based on Benford's Law formula
    print("\nExpected Benford Distribution:")
    print(expected)
