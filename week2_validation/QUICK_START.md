# Week 2 Validation - Quick Start Guide

**Run the pipeline on your fusion data.** No demo datasets included — provide your own data file.

---

## ⚡ Run the Pipeline

### **Recommended: run_week2.py (full control)**

```bash
python -m week2_validation.run_week2 \
    --fusion-data path/to/your_fusion_data.csv \
    --output-dir output_name \
    --run-all \
    --generate-report \
    --generate-pdf
```

**From project root** (e.g. `HumanG_Fusion/`). Replace `path/to/your_fusion_data.csv` with your file and `output_name` with your output folder name.

---

### **Simplified: run_complete.py (pipeline + tests)**

```bash
python week2_validation/run_complete.py path/to/your_fusion_data.csv output_name
```

This runs the validation pipeline and the test suite, then verifies outputs.

---

### **Batch / shell scripts**

**Windows:**
```batch
run_validation.bat path\to\your_fusion_data.csv output_name
```

**Linux/Mac:**
```bash
chmod +x week2_validation/run_validation.sh
./week2_validation/run_validation.sh path/to/your_fusion_data.csv output_name
```

---

## 📋 Input Requirements

- **Formats:** CSV, TSV, JSON, Parquet, or XLSX
- **Required columns:** `fusion_id`, `gene_1`, `gene_2`, `protein_length`, `recurrence_count`
- **Constraints:** `protein_length` > 0, `recurrence_count` >= 0, no nulls in required columns

---

## 📁 Output Files

Outputs are written to `week2_validation/results/output_name/` (or the name you chose):

```
output_name/
├── validation_status_{stem}.json              ← Pipeline status
├── validation_dataset_certification_{stem}.json  ← Certification
├── validation_diagnostic_results_{stem}.json    ← All metrics
├── validation_cleaned_dataset_{stem}.csv        ← Cleaned data
├── Statistical_Integrity_Report_{stem}.md       ← Narrative report
├── data_integrity_validation_report.pdf         ← PDF report (with --generate-pdf)
├── data_integrity_validation_report.html        ← HTML report
├── protein_distribution.png                     ← Histogram
├── skewness_diagnostic.png                      ← Skewness gauge
├── benford_analysis.png                         ← Benford plot
└── run_log_{stem}.txt                           ← Execution log
```

`{stem}` is derived from your input filename (e.g. `my_fusion.parquet` → `my_fusion`).

---

## 🧬 COSMIC Reference

The pipeline uses **real COSMIC Fusion v103 GRCh38** by default. The file must exist at:

```
week2_validation/cosmic/Cosmic_Fusion_v103_GRCh38.tsv
```

To use a different COSMIC file:

```bash
python -m week2_validation.run_week2 \
    --fusion-data path/to/your_data.csv \
    --output-dir output_name \
    --cosmic-data path/to/your_cosmic.tsv \
    --run-all
```

---

## 🔧 Troubleshooting

### **"Module not found"**
```bash
# Run from project root
cd /path/to/HumanG_Fusion
python -m week2_validation.run_week2 --help
```

### **"No such file or directory"**
Ensure your fusion data file path is correct and the file exists.

### **"COSMIC reference file not found"**
Ensure `Cosmic_Fusion_v103_GRCh38.tsv` exists in `week2_validation/cosmic/`, or provide `--cosmic-data`.

---

## 📚 More Information

- **Full command reference:** `week2_validation/RUN_COMMANDS.md`
- **Statistical methodology:** `week2_validation/cosmic/statistical_methodology.md`
- **Main README:** `week2_validation/README.md`
