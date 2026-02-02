# COSMIC Validation Statistical Methodology

*Week 2: Data Integrity & Statistical Validation*

---

## Overview

This document provides transparent documentation of the statistical methods used in COSMIC cross-validation for fusion recurrence data. All methods are chosen for their appropriateness to biological recurrence count data, which is typically non-normal, heavy-tailed, and rank-ordered.

---

## 1. Spearman Rank Correlation

### Why Spearman Correlation?

We use Spearman's rank correlation coefficient (ρ) rather than Pearson's correlation for the following reasons:

1. **Non-Parametric**: Spearman correlation does not assume data follows a normal distribution. Fusion recurrence counts are highly non-normal, typically exhibiting heavy right tails with most fusions having low counts and few having very high counts.

2. **Robust to Outliers**: Recurrence data often contains extreme values (e.g., BCR-ABL1 with very high recurrence). Spearman uses ranks rather than raw values, making it robust to these outliers.

3. **Appropriate for Ordinal Comparison**: We are comparing whether highly recurrent fusions in the dataset are also highly recurrent in COSMIC (rank-order agreement), not whether exact count values match. Spearman directly measures this rank-order agreement.

4. **Heavy-Tailed Distributions**: Biological recurrence data follows power-law/Zipf-like distributions. Spearman correlation is appropriate for such distributions as it captures monotonic relationships without being distorted by the long tail.

### Interpretation

- ρ ≈ 1.0: Strong positive rank agreement (high recurrence in dataset → high recurrence in COSMIC)
- ρ ≈ 0.0: No rank relationship
- ρ ≈ -1.0: Inverse rank relationship (rare in practice for recurrence data)

### Limitations

- Spearman correlation only measures monotonic relationships; it does not capture non-monotonic patterns.
- Statistical significance (p-value) assumes independent observations, which may not strictly hold for biologically related fusions.

---

## 2. Hypergeometric Enrichment Test

### Why Hypergeometric Test?

We use the hypergeometric test to assess whether top-ranked fusions in the dataset overlap with COSMIC's top fusions more than expected by chance.

1. **Enrichment vs Random Expectation**: The test compares observed overlap to what would be expected if fusions were randomly drawn from the population.

2. **Discrete Count Model**: The hypergeometric distribution is the correct model for sampling without replacement, which matches our scenario (a fusion cannot appear twice in the same list).

3. **Conservative for Enrichment**: The one-tailed p-value tests specifically whether enrichment is greater than chance, appropriate for validating biological consistency.

### Test Setup

- **Population size (N)**: Total number of unique fusion pairs across both datasets
- **Success states (K)**: Number of fusions in COSMIC's top-N
- **Draws (n)**: Number of fusions in dataset's top-N
- **Observed successes (k)**: Number of overlapping fusions in top-N

### Independence Limitation

**Important caveat**: The hypergeometric test assumes independence between fusion events. In biological reality, fusions involving genes in the same pathway or functional network are not independent. Pathway-level correlations and functional relationships introduce dependency between fusion events.

**Conservative interpretation recommended**: P-values should be treated as approximate rather than exact. Significant enrichment (p < 0.05) provides supportive evidence for biological consistency, but should not be over-interpreted as precise statistical proof.

---

## 3. Bootstrap Confidence Intervals

### Why Bootstrap CI for Spearman?

We provide bootstrap confidence intervals for the Spearman correlation coefficient because:

1. **Distribution-Free Uncertainty Quantification**: Bootstrap does not assume a specific sampling distribution for the correlation coefficient.

2. **Appropriate for Small Overlap**: When overlap between datasets is small, parametric methods for correlation CI can be unreliable.

3. **Transparent Uncertainty**: CI bounds communicate the precision of the correlation estimate, not just a point estimate.

### Method

1. Resample overlapping fusion pairs with replacement (B = 1000 iterations)
2. Compute Spearman ρ for each resample
3. Report the 2.5th and 97.5th percentiles as the 95% CI bounds

### Interpretation

If the CI contains 0, there is insufficient evidence for consistent rank agreement. Narrow CI indicates a precise estimate; wide CI indicates uncertainty.

---

## 4. Negative Control Validation

### Purpose

The negative control test validates that observed correlation is not an artifact of the analysis method or random chance.

### Method

1. Shuffle COSMIC recurrence counts (preserving distribution but breaking biological relationships)
2. Recompute Spearman correlation with shuffled data
3. Repeat 1000 times to establish null distribution
4. Compare real correlation to shuffled distribution

### Interpretation

If real ρ >> mean(shuffled ρ), this provides evidence that the observed agreement reflects genuine biological consistency rather than statistical artifact.

---

## 5. Composite Quality Score

### Components

The COSMIC validation score combines multiple independent signals:

| Component | Weight | Rationale |
|-----------|--------|-----------|
| Spearman strength | 0.30 | Primary measure of rank agreement |
| Spearman significance | 0.20 | Statistical confidence in correlation |
| Enrichment significance | 0.20 | Top-fusion overlap beyond chance |
| Negative control delta | 0.20 | Real vs shuffled correlation difference |
| Overlap enrichment ratio | 0.10 | Ratio of observed to expected overlap |

### Classification Thresholds

| Score Range | Classification |
|-------------|----------------|
| ≥ 0.70 | STRONG BIOLOGICAL AGREEMENT |
| 0.50 - 0.69 | MODERATE AGREEMENT |
| 0.30 - 0.49 | WEAK SIGNAL |
| < 0.30 | NO MEANINGFUL COSMIC AGREEMENT |

### Transparency

All component scores are reported individually, allowing researchers to understand which aspects of the validation contribute to the final classification.

---

## 6. Statistical Assumptions Summary

### Assumptions Made

1. Recurrence counts are comparable across datasets (same counting methodology)
2. Gene symbols can be meaningfully matched across datasets (alias normalization applied)
3. Overlapping fusions represent the same biological entity

### Assumptions NOT Made

1. ❌ Normal distribution of counts
2. ❌ Equal sample sizes
3. ❌ Independence of pathway-related fusions (acknowledged limitation)
4. ❌ Exact equivalence of COSMIC and dataset cohorts

---

## 7. Method Justification for Fusion Recurrence Comparison

### Why These Methods Are Appropriate

1. **Biological Context**: Fusion recurrence data inherently follows a heavy-tailed distribution where few fusions are highly recurrent and many are rare. Non-parametric rank-based methods are optimal for this distribution.

2. **Purpose Alignment**: Our goal is to assess whether recurrence rankings are consistent with COSMIC (plausibility check), not to prove exact equivalence. Rank correlation directly addresses this.

3. **Robustness**: The combination of Spearman correlation, bootstrap CI, hypergeometric enrichment, and negative control provides multiple independent lines of evidence, reducing the risk of false conclusions from any single metric.

4. **Transparency**: All metrics and their components are reported, allowing domain experts to apply their judgment to the results.

---

## References

1. Spearman, C. (1904). "The proof and measurement of association between two things". American Journal of Psychology.
2. Efron, B. (1979). "Bootstrap Methods: Another Look at the Jackknife". Annals of Statistics.
3. Fisher, R. A. (1935). "The logic of inductive inference". Journal of the Royal Statistical Society.

---

*Document Version: 1.0*
*Last Updated: Generated automatically during validation pipeline*
