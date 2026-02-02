# Week 2 Validation - Quick Start Guide

**The Simplest Way to Run Everything** 🚀

---

## ⚡ One-Command Solutions

### **Option 1: Python Script (Recommended - All Platforms)**

```bash
python week2_validation/run_complete.py
```

**That's it!** This single command:
- ✅ Runs the complete validation pipeline with all diagnostics
- ✅ Runs all 33 tests
- ✅ Verifies outputs are generated
- ✅ Checks new scientific defensibility features
- ✅ Provides a clear success/failure report

---

### **Option 2: Batch File (Windows - Double-Click)**

**Just double-click:** `week2_validation/run_validation.bat`

No terminal needed! The batch file will:
- Run everything automatically
- Show progress in a command window
- Wait for you to press a key when done

---

### **Option 3: Shell Script (Linux/Mac)**

```bash
chmod +x week2_validation/run_validation.sh
./week2_validation/run_validation.sh
```

---

## 📋 Usage Examples

### **Basic Usage (Demo Data)**

```bash
# Uses demo_fusion.csv, outputs to demo_output/
python week2_validation/run_complete.py
```

### **Custom Input File**

```bash
# Use your own data
python week2_validation/run_complete.py path/to/your_data.csv
```

### **Custom Input and Output**

```bash
# Specify both input and output
python week2_validation/run_complete.py path/to/data.csv output_folder/
```

### **With COSMIC Reference Data**

```bash
# Use real COSMIC instead of mock
python week2_validation/run_complete.py data.csv output/ --cosmic cosmic.tsv
```

### **Skip Tests (Pipeline Only)**

```bash
# Run pipeline but skip test suite
python week2_validation/run_complete.py --skip-tests
```

### **Get Help**

```bash
python week2_validation/run_complete.py --help
```

---

## 📊 What Gets Run?

### **Validation Pipeline:**
1. **Distribution Diagnostics** - Skewness, kurtosis, quantiles
2. **Benford's Law Analysis** - First digit frequency distribution
3. **Log-Normality Tests** - Kolmogorov-Smirnov, Anderson-Darling
4. **Benford Statistical Controls** - Control validation
5. **COSMIC Cross-Validation** with:
   - Bootstrap confidence intervals (1000 iterations)
   - Reproducibility lock metadata
   - Score component breakdown
   - Mock fallback if no COSMIC data

### **Test Suite:**
- All 33 tests including:
  - Core validation tests
  - COSMIC validation tests
  - **New**: Bootstrap CI tests
  - **New**: Reproducibility lock tests
  - **New**: Mock distribution validation
  - **New**: Score breakdown validation
  - Cross-format compatibility tests

### **Verification:**
- Checks all expected output files exist
- Verifies new scientific defensibility features
- Reports file sizes and status

---

## 📁 Output Files

After running, check `demo_output/` (or your custom output folder):

```
demo_output/
├── week2_status_demo_fusion.json                    ← Pipeline status
├── week2_dataset_certification_demo_fusion.json     ← Certification
├── week2_diagnostic_results_demo_fusion.json        ← All metrics
├── week2_cleaned_dataset_demo_fusion.csv            ← Cleaned data
├── Statistical_Integrity_Report_demo_fusion.md      ← Human-readable report
├── protein_distribution.png                         ← Histogram
├── skewness_diagnostic.png                          ← Skewness gauge
└── run_log_demo_fusion.txt                          ← Execution log
```

---

## 🎯 Expected Output

```
================================================================================
  Week 2 Validation - Complete Pipeline Runner
================================================================================

📂 Fusion data: week2_validation/demo_fusion.csv
📁 Output directory: demo_output
🧬 COSMIC data: Using mock fallback

--------------------------------------------------------------------------------
  Step 1: Running Validation Pipeline
--------------------------------------------------------------------------------

🚀 Validation pipeline...
   Command: python -m week2_validation.run_week2 --fusion-data week2_validation/demo_fusion.csv --output-dir demo_output --run-all

✅ Validation pipeline completed successfully

--------------------------------------------------------------------------------
  Step 2: Running Test Suite
--------------------------------------------------------------------------------

🚀 Test suite...
   Command: python -m pytest week2_validation/tests/ -v --tb=short

✅ Test suite completed successfully

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Verifying Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ✅ week2_status_demo_fusion.json (2,847 bytes)
  ✅ week2_dataset_certification_demo_fusion.json (721 bytes)
  ✅ week2_diagnostic_results_demo_fusion.json (5,432 bytes)
  ✅ week2_cleaned_dataset_demo_fusion.csv (12,458 bytes)
  ✅ Statistical_Integrity_Report_demo_fusion.md (18,234 bytes)
  ✅ protein_distribution.png (45,672 bytes)
  ✅ skewness_diagnostic.png (38,291 bytes)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Verifying New Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ✅ Bootstrap CI: [0.8234, 0.9567]
  ✅ Reproducibility Lock:
     - Python: 3.11.3
     - NumPy: 1.24.3
     - Code version: 1.1.0
  ✅ Score Component Breakdown:
     - Correlation: 0.50
     - Enrichment: 0.20
     - Negative Control: 0.15
     - Overlap: 0.10

================================================================================
  Final Summary
================================================================================

🎉 SUCCESS! Complete validation finished.
📊 All outputs generated in: demo_output/
📄 View report: demo_output/Statistical_Integrity_Report_demo_fusion.md
```

---

## 🔧 Troubleshooting

### **"Module not found" error**
```bash
# Make sure you're in the project root
cd C:\Users\enchi\Documents\HumanG_Fusion
python week2_validation/run_complete.py
```

### **"No such file or directory"**
```bash
# Check the demo file exists
dir week2_validation\demo_fusion.csv   # Windows
ls week2_validation/demo_fusion.csv    # Linux/Mac
```

### **Tests fail**
This is OK! Some tests may fail if:
- Demo data has no overlap with mock COSMIC (expected)
- Optional dependencies missing (will be skipped)

The pipeline still works correctly.

---

## 📚 More Information

- **Full command reference:** `week2_validation/RUN_COMMANDS.md`
- **Implementation details:** `week2_validation/IMPLEMENTATION_SUMMARY.md`
- **Statistical methodology:** `week2_validation/cosmic/statistical_methodology.md`
- **Main README:** `week2_validation/README.md`

---

## 💡 Pro Tips

1. **First time?** Just run: `python week2_validation/run_complete.py`
2. **Want to customize?** See `week2_validation/RUN_COMMANDS.md`
3. **Need help?** Run: `python week2_validation/run_complete.py --help`
4. **View the report:** Open `demo_output/Statistical_Integrity_Report_demo_fusion.md`

---

**Bottom line:** The simplest command to run everything is:

```bash
python week2_validation/run_complete.py
```

No arguments needed. No complex syntax. Just works. ✅
