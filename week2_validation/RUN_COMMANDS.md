# Week 2 Validation: Complete Run Commands

**Last Updated:** February 2026  
**COSMIC Migration:** Real COSMIC Fusion v103 GRCh38 now used as default (no mock fallback)

---

## ⚡ Quick Reference

**Most Common Command (from project root):**
```bash
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all
```

**Key Points:**
- ✅ Uses **real COSMIC Fusion v103 GRCh38** by default (no mock fallback)
- ✅ Supports CSV, TSV, JSON, Parquet, XLSX input formats
- ✅ All outputs written to `--output-dir`
- ✅ Run from project root (e.g. `HumanG_Fusion/`)
- ✅ **No demo datasets included** — provide your own fusion data file. For sample data, run `python week2_validation/scripts/generate_demo_files.py`

---

## 🚀 End-to-End Complete Pipeline Run

### **Full Pipeline with All Diagnostics & Tests**

```bash
# Run the complete Week 2 validation pipeline with all diagnostics
# From project root (e.g. HumanG_Fusion/)
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all

# Alternative using module syntax (requires PYTHONPATH from project root):
# Windows PowerShell:
$env:PYTHONPATH = (Get-Location).Path
python -m week2_validation.run_week2 \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all

# Linux/macOS:
export PYTHONPATH=$(pwd)
python -m week2_validation.run_week2 \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all
```

**What `--run-all` does:**
- ✅ Runs distribution diagnostics (skewness, kurtosis, quantiles)
- ✅ Runs Benford's Law analysis (first digit frequency)
- ✅ Runs log-normality tests (KS, Anderson-Darling)
- ✅ Runs Benford statistical controls
- ✅ Runs COSMIC cross-validation (uses real COSMIC Fusion v103 GRCh38 by default)
- ✅ Generates all reports and visualizations
- ✅ Creates certification JSON
- ✅ Produces narrative Markdown report

**COSMIC Reference Data:**
- **Default:** Automatically loads `week2_validation/cosmic/Cosmic_Fusion_v103_GRCh38.tsv`
- **Custom:** Use `--cosmic-data` flag to provide your own COSMIC file
- **No Fallback:** Pipeline will fail if COSMIC cannot be loaded (ensures data quality)

---

## 📊 Generated Outputs

After running `--run-all`, you'll get these files in your output directory:

### **JSON Outputs:**
```
output_dir/
├── week2_status_{stem}.json                    # Pipeline execution status
├── week2_dataset_certification_{stem}.json     # Dataset approval certification
├── week2_diagnostic_results_{stem}.json        # All diagnostic metrics (includes COSMIC results)
└── (COSMIC validation metrics included in diagnostic_results JSON)
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
# From project root (HumanG_Fusion/)
python -m pytest week2_validation/tests/ -v

# Run with coverage report
python -m pytest week2_validation/tests/ -v --cov=week2_validation --cov-report=html

# Run only COSMIC validation tests (including real COSMIC loading)
python -m pytest week2_validation/tests/test_cosmic_validation.py -v

# Run specific feature tests
python -m pytest week2_validation/tests/test_cosmic_validation.py::test_real_cosmic_loading \
    week2_validation/tests/test_cosmic_validation.py::test_bootstrap_ci_contains_rho \
    week2_validation/tests/test_cosmic_validation.py::test_reproducibility_lock_metadata \
    week2_validation/tests/test_cosmic_validation.py::test_score_breakdown_sums_correctly \
    week2_validation/tests/test_cosmic_validation.py::test_statistical_methodology_docs_exist \
    -v
```

---

## 🎯 Complete End-to-End Workflow (Pipeline + Tests)

### **Option 1: Sequential Execution**

```bash
# Step 1: Run full validation pipeline
# From project root (e.g. HumanG_Fusion/)
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all

# Step 2: Run test suite
python -m pytest week2_validation/tests/ -v

# Step 3: Verify outputs exist
# Windows PowerShell:
Get-ChildItem output_name/
# Linux/Mac:
ls output_name/
```

### **Option 2: PowerShell Script (Windows)**

```powershell
# Create a PowerShell script: run_complete_validation.ps1
# Run from project root (HumanG_Fusion/)

# Set PYTHONPATH if needed
$env:PYTHONPATH = (Get-Location).Path

# Run pipeline
python week2_validation/run_week2.py `
    --fusion-data path/to/your_fusion_data.csv `
    --output-dir output_name `
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
Get-ChildItem output_name/ | Format-Table Name, Length, LastWriteTime
```

### **Option 3: Bash Script (Linux/Mac)**

```bash
#!/bin/bash
# Create a bash script: run_complete_validation.sh
# Run from project root (HumanG_Fusion/)

set -e  # Exit on error

echo "🚀 Running Week 2 Complete Validation Pipeline..."

# Run pipeline
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all

echo "✅ Pipeline completed successfully"

# Run tests
echo "🧪 Running test suite..."
python -m pytest week2_validation/tests/ -v

echo "✅ All tests passed"

# List outputs
echo ""
echo "📁 Generated outputs:"
ls -lh output_name/

echo ""
echo "🎉 Complete validation workflow finished successfully!"
```

---

## 🔬 Advanced Options

### **With Custom COSMIC Data**

```bash
# Use your own COSMIC reference data instead of default
# Default is: week2_validation/cosmic/Cosmic_Fusion_v103_GRCh38.tsv
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --cosmic-data path/to/your_cosmic_data.tsv \
    --run-all

# COSMIC file can be in COSMIC Fusion format (FIVE_PRIME_GENE_SYMBOL, THREE_PRIME_GENE_SYMBOL)
# or standard format (gene_1, gene_2, recurrence_count)
# Pipeline automatically transforms COSMIC Fusion format to standard format
```

### **Individual Diagnostic Components**

```bash
# Run only specific diagnostics (instead of --run-all)
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-diagnostics \
    --run-benford \
    --run-cosmic

# Just COSMIC validation (uses default real COSMIC)
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-cosmic

# COSMIC validation with custom COSMIC file
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-cosmic \
    --cosmic-data path/to/your_cosmic_data.tsv
```

### **Dry Run (Validation Only)**

```bash
# Check inputs without running diagnostics
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --dry-run
```

---

## 📋 Supported Input Formats

The pipeline supports multiple file formats:

```bash
# CSV
python week2_validation/run_week2.py --fusion-data data.csv --output-dir output --run-all

# TSV
python week2_validation/run_week2.py --fusion-data data.tsv --output-dir output --run-all

# JSON
python week2_validation/run_week2.py --fusion-data data.json --output-dir output --run-all

# Parquet
python week2_validation/run_week2.py --fusion-data data.parquet --output-dir output --run-all

# Excel
python week2_validation/run_week2.py --fusion-data data.xlsx --output-dir output --run-all
```

---

## 📊 Example: Complete Workflow with Multiple Files

```bash
# Process multiple fusion datasets
for file in data/*.csv; do
    stem=$(basename "$file" .csv)
    echo "Processing $stem..."
    python week2_validation/run_week2.py \
        --fusion-data "$file" \
        --output-dir "output_$stem" \
        --run-all
done

# Run tests after all processing
python -m pytest week2_validation/tests/ -v
```

**PowerShell version:**
```powershell
# Process multiple fusion datasets
Get-ChildItem data/*.csv | ForEach-Object {
    $stem = $_.BaseName
    Write-Host "Processing $stem..."
    python week2_validation/run_week2.py `
        --fusion-data $_.FullName `
        --output-dir "output_$stem" `
        --run-all
}

# Run tests after all processing
python -m pytest week2_validation/tests/ -v
```

---

## 🔍 Verify Implementation Features

### **Check Scientific Defensibility Features**

```bash
# Verify statistical methodology documentation exists
cat week2_validation/cosmic/statistical_methodology.md | head -50

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
   mkdir -p output_name
   chmod 755 output_name
   ```

4. **COSMIC data not found:**
   ```bash
   # Pipeline requires COSMIC data (default: Cosmic_Fusion_v103_GRCh38.tsv)
   # Ensure file exists at: week2_validation/cosmic/Cosmic_Fusion_v103_GRCh38.tsv
   # Or provide explicit COSMIC file:
   python week2_validation/run_week2.py \
       --fusion-data data.csv \
       --output-dir output \
       --cosmic-data path/to/cosmic.tsv \
       --run-all
   
   # If COSMIC file is missing, pipeline will raise RuntimeError
   # This ensures data quality (no silent fallback to mock data)
   ```

---

## 📚 Related Documentation

- **Statistical Methodology:** `week2_validation/cosmic/statistical_methodology.md`
- **Main README:** `week2_validation/README.md`

---

## ✅ Quick Checklist

Before running end-to-end validation:

- [ ] Python 3.8+ installed
- [ ] All dependencies installed (`pip install -r week2_validation/requirements.txt`)
- [ ] Input fusion data file exists
- [ ] Output directory is writable
- [ ] Real COSMIC file exists: `week2_validation/cosmic/Cosmic_Fusion_v103_GRCh38.tsv` (or use `--cosmic-data` flag)

To run complete validation:

- [ ] Run pipeline with `--run-all` flag
- [ ] Run test suite with `pytest`
- [ ] Verify all outputs generated
- [ ] Review narrative report
- [ ] Check certification status

---

## 🎓 Example Output Inspection

Output filenames use `{stem}` from your input filename (e.g. `your_fusion_data.csv` → stem is `your_fusion_data`).

```bash
# After running pipeline, inspect key outputs (replace {stem} with your input filename stem):

# 1. Check pipeline status
cat output_name/validation_status_{stem}.json | grep -E '"status"|"exit_code"'

# 2. Check certification
cat output_name/validation_dataset_certification_{stem}.json | grep "approved_for_powerlaw_modeling"

# 3. View narrative report
cat output_name/Statistical_Integrity_Report_{stem}.md | less

# 4. Check COSMIC validation score
cat output_name/validation_diagnostic_results_{stem}.json | grep -A 5 "cosmic_validation"

# 5. Verify reproducibility metadata
cat output_name/validation_diagnostic_results_{stem}.json | grep -A 10 "reproducibility_lock"
```

---

**For the absolute simplest end-to-end run:**

```bash
# From project root (e.g. HumanG_Fusion/)
# Replace path/to/your_fusion_data.csv with your actual fusion data file
python week2_validation/run_week2.py \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all && \
python -m pytest week2_validation/tests/ -v
```

**PowerShell version:**
```powershell
# From project root (e.g. HumanG_Fusion/)
python week2_validation/run_week2.py `
    --fusion-data path/to/your_fusion_data.csv `
    --output-dir output_name `
    --run-all; if ($LASTEXITCODE -eq 0) { python -m pytest week2_validation/tests/ -v }
```

This command runs the entire pipeline with all diagnostics (including real COSMIC validation) and then runs all tests. ✅

**What happens:**
1. Pipeline loads your fusion data file
2. Automatically loads real COSMIC Fusion v103 GRCh38 reference data
3. Runs all diagnostics (distribution, Benford, log-normality, COSMIC)
4. Generates all outputs (JSON, CSV, Markdown reports, PDF/HTML)
5. Runs full test suite
