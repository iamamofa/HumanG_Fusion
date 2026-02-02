# Week 2 Validation: Complete Run Commands

**Last Updated:** February 2, 2026

---

## 🚀 End-to-End Complete Pipeline Run

### **Full Pipeline with All Diagnostics & Tests**

```bash
# Run the complete Week 2 validation pipeline with all diagnostics
python -m week2_validation.run_week2 \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir path/to/output_directory \
    --run-all

# Example with actual file:
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --run-all
```

**What `--run-all` does:**
- ✅ Runs distribution diagnostics (skewness, kurtosis, quantiles)
- ✅ Runs Benford's Law analysis (first digit frequency)
- ✅ Runs log-normality tests (KS, Anderson-Darling)
- ✅ Runs Benford statistical controls
- ✅ Runs COSMIC cross-validation (with mock fallback if no COSMIC data provided)
- ✅ Generates all reports and visualizations
- ✅ Creates certification JSON
- ✅ Produces narrative Markdown report

---

## 📊 Generated Outputs

After running `--run-all`, you'll get these files in your output directory:

### **JSON Outputs:**
```
output_dir/
├── week2_status_{stem}.json                    # Pipeline execution status
├── week2_dataset_certification_{stem}.json     # Dataset approval certification
├── week2_diagnostic_results_{stem}.json        # All diagnostic metrics
└── mock_cosmic_distribution_validation.json    # Mock COSMIC validation (if mock used)
```

### **CSV Outputs:**
```
output_dir/
└── week2_cleaned_dataset_{stem}.csv            # Cleaned dataset
```

### **Visualizations (PNG):**
```
output_dir/
├── protein_distribution.png                    # Histogram of protein lengths
└── skewness_diagnostic.png                     # Skewness gauge visualization
```

### **Reports (Markdown):**
```
output_dir/
└── Statistical_Integrity_Report_{stem}.md      # Human-readable narrative report
```

### **Logs:**
```
output_dir/
└── run_log_{stem}.txt                          # Detailed execution log
```

---

## 🧪 Run Complete Test Suite

### **Run All Tests (Including New Scientific Defensibility Tests)**

```bash
# Run all Week 2 validation tests
python -m pytest week2_validation/tests/ -v

# Run with coverage report
python -m pytest week2_validation/tests/ -v --cov=week2_validation --cov-report=html

# Run only COSMIC validation tests (including new features)
python -m pytest week2_validation/tests/test_cosmic_validation.py -v

# Run specific new feature tests
python -m pytest week2_validation/tests/test_cosmic_validation.py::test_bootstrap_ci_contains_rho \
    week2_validation/tests/test_cosmic_validation.py::test_reproducibility_lock_metadata \
    week2_validation/tests/test_cosmic_validation.py::test_mock_distribution_metrics_computed \
    week2_validation/tests/test_cosmic_validation.py::test_score_breakdown_sums_correctly \
    week2_validation/tests/test_cosmic_validation.py::test_statistical_methodology_docs_exist \
    -v
```

---

## 🎯 Complete End-to-End Workflow (Pipeline + Tests)

### **Option 1: Sequential Execution**

```bash
# Step 1: Run full validation pipeline
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --run-all

# Step 2: Run test suite
python -m pytest week2_validation/tests/ -v

# Step 3: Verify outputs exist
ls demo_output/
```

### **Option 2: PowerShell Script (Windows)**

```powershell
# Create a PowerShell script: run_complete_validation.ps1

# Run pipeline
python -m week2_validation.run_week2 `
    --fusion-data week2_validation/demo_fusion.csv `
    --output-dir demo_output `
    --run-all

# Check exit code
if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Pipeline completed successfully" -ForegroundColor Green
} else {
    Write-Host "❌ Pipeline failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

# Run tests
python -m pytest week2_validation/tests/ -v

# Check test results
if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ All tests passed" -ForegroundColor Green
} else {
    Write-Host "❌ Tests failed" -ForegroundColor Red
    exit $LASTEXITCODE
}

# List outputs
Write-Host "`n📁 Generated outputs:" -ForegroundColor Cyan
Get-ChildItem demo_output/ | Format-Table Name, Length, LastWriteTime
```

### **Option 3: Bash Script (Linux/Mac)**

```bash
#!/bin/bash
# Create a bash script: run_complete_validation.sh

set -e  # Exit on error

echo "🚀 Running Week 2 Complete Validation Pipeline..."

# Run pipeline
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --run-all

echo "✅ Pipeline completed successfully"

# Run tests
echo "🧪 Running test suite..."
python -m pytest week2_validation/tests/ -v

echo "✅ All tests passed"

# List outputs
echo ""
echo "📁 Generated outputs:"
ls -lh demo_output/

echo ""
echo "🎉 Complete validation workflow finished successfully!"
```

---

## 🔬 Advanced Options

### **With Custom COSMIC Data**

```bash
# Use your own COSMIC reference data instead of mock
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --cosmic-data path/to/cosmic_census.tsv \
    --run-all
```

### **Individual Diagnostic Components**

```bash
# Run only specific diagnostics (instead of --run-all)
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --run-diagnostics \
    --run-benford \
    --run-cosmic

# Just COSMIC validation
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --run-cosmic \
    --cosmic-data path/to/cosmic_census.tsv
```

### **Dry Run (Validation Only)**

```bash
# Check inputs without running diagnostics
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --dry-run
```

---

## 📋 Supported Input Formats

The pipeline supports multiple file formats:

```bash
# CSV
python -m week2_validation.run_week2 --fusion-data data.csv --output-dir output --run-all

# TSV
python -m week2_validation.run_week2 --fusion-data data.tsv --output-dir output --run-all

# JSON
python -m week2_validation.run_week2 --fusion-data data.json --output-dir output --run-all

# Parquet
python -m week2_validation.run_week2 --fusion-data data.parquet --output-dir output --run-all

# Excel
python -m week2_validation.run_week2 --fusion-data data.xlsx --output-dir output --run-all
```

---

## 🧬 Generate Mock COSMIC Validation Report

### **Standalone Mock COSMIC Realism Validation**

```bash
# Generate mock COSMIC distribution validation report
python week2_validation/cosmic/mock_cosmic_validation.py \
    path/to/mock_cosmic_census.csv \
    output_directory/

# Example:
python week2_validation/cosmic/mock_cosmic_validation.py \
    week2_validation/cosmic/mock_cosmic_census.csv \
    demo_output/
```

**Output:** `output_directory/mock_cosmic_distribution_validation.json`

**Contains:**
- Skewness, kurtosis, tail heaviness
- Gini coefficient (inequality measure)
- Zero inflation rate
- Realism assessment (REALISTIC/ACCEPTABLE/DEGRADED)
- Compatibility statement

---

## 📊 Example: Complete Workflow with Multiple Files

```bash
# Process multiple fusion datasets
for file in data/*.csv; do
    stem=$(basename "$file" .csv)
    echo "Processing $stem..."
    python -m week2_validation.run_week2 \
        --fusion-data "$file" \
        --output-dir "output_$stem" \
        --run-all
done

# Run tests after all processing
python -m pytest week2_validation/tests/ -v
```

---

## 🔍 Verify Implementation Features

### **Check New Scientific Defensibility Features**

```bash
# Verify statistical methodology documentation exists
cat week2_validation/cosmic/statistical_methodology.md | head -50

# Verify mock validation module exists
python -c "import week2_validation.cosmic.mock_cosmic_validation as m; print('✅ Mock validation module loaded')"

# Verify bootstrap CI is available
python -c "from week2_validation.cosmic.diagnostics import compute_bootstrap_spearman_ci; print('✅ Bootstrap CI available')"

# Check version
python -c "import week2_validation.cosmic.diagnostics as d; print(f'✅ Code version: {d.COSMIC_VALIDATION_CODE_VERSION}')"
```

---

## 📈 Performance Notes

**Typical Execution Times:**
- Small dataset (< 100 rows): ~5-10 seconds
- Medium dataset (100-1000 rows): ~10-30 seconds
- Large dataset (1000+ rows): ~30-60 seconds
- Bootstrap CI computation: +2-5 seconds per dataset
- Full test suite: ~5-10 seconds

**Resource Usage:**
- Memory: < 500 MB for typical datasets
- CPU: Single-threaded (bootstrap uses NumPy vectorization)
- Disk: Minimal (< 10 MB per run for outputs)

---

## 🛠️ Troubleshooting

### **Common Issues:**

1. **Module not found:**
   ```bash
   # Ensure you're in the project root
   cd /path/to/HumanG_Fusion
   python -m week2_validation.run_week2 --help
   ```

2. **Missing dependencies:**
   ```bash
   pip install pandas numpy scipy matplotlib openpyxl pyarrow
   ```

3. **Permission denied on output:**
   ```bash
   # Ensure output directory is writable
   mkdir -p demo_output
   chmod 755 demo_output
   ```

4. **COSMIC data not found:**
   ```bash
   # Pipeline will use mock COSMIC automatically
   # Or provide explicit COSMIC file:
   python -m week2_validation.run_week2 \
       --fusion-data data.csv \
       --output-dir output \
       --cosmic-data cosmic.tsv \
       --run-all
   ```

---

## 📚 Related Documentation

- **Implementation Summary:** `week2_validation/IMPLEMENTATION_SUMMARY.md`
- **Statistical Methodology:** `week2_validation/cosmic/statistical_methodology.md`
- **Main README:** `week2_validation/README.md`
- **COSMIC Verification:** `week2_validation/COSMIC_OPERATIONAL_VERIFICATION_FINAL_REPORT.md`

---

## ✅ Quick Checklist

Before running end-to-end validation:

- [ ] Python 3.8+ installed
- [ ] All dependencies installed (`pip install -r requirements.txt`)
- [ ] Input fusion data file exists
- [ ] Output directory is writable
- [ ] (Optional) COSMIC reference data available

To run complete validation:

- [ ] Run pipeline with `--run-all` flag
- [ ] Run test suite with `pytest`
- [ ] Verify all outputs generated
- [ ] Review narrative report
- [ ] Check certification status

---

## 🎓 Example Output Inspection

```bash
# After running pipeline, inspect key outputs:

# 1. Check pipeline status
cat demo_output/week2_status_demo_fusion.json | grep -E '"status"|"exit_code"'

# 2. Check certification
cat demo_output/week2_dataset_certification_demo_fusion.json | grep "approved_for_powerlaw_modeling"

# 3. View narrative report
cat demo_output/Statistical_Integrity_Report_demo_fusion.md | less

# 4. Check COSMIC validation score
cat demo_output/week2_diagnostic_results_demo_fusion.json | grep -A 5 "cosmic_validation"

# 5. Verify reproducibility metadata
cat demo_output/week2_diagnostic_results_demo_fusion.json | grep -A 10 "reproducibility_lock"
```

---

**For the absolute simplest end-to-end run:**

```bash
python -m week2_validation.run_week2 \
    --fusion-data week2_validation/demo_fusion.csv \
    --output-dir demo_output \
    --run-all && \
python -m pytest week2_validation/tests/ -v
```

This single command runs the entire pipeline with all diagnostics and then runs all tests. ✅
