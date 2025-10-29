#  RNA-Seq & Fusion Detection Reference Setup (GRCh38)

This repository provides **automated setup scripts** to prepare reference resources for RNA-seq alignment, fusion detection, and contamination screening tools such as **STAR**, **STAR-Fusion**, and **Kraken2**.

Two setup options are provided:
1. **Full Setup** – for high-performance systems or production workflows.
2. **Lightweight Setup** – for testing, limited storage, or learning environments.

---

## Repository Contents

| File | Description |
|------|--------------|
| `setup_references.sh` | Full setup: downloads complete genome, annotations, CTAT library, and standard Kraken2 DB (~150 GB). |
| `setup_references_light.sh` | Lightweight setup: uses a subset of GRCh38 (chr21 & chr22) and small Kraken2 DB (~10–15 GB). |
| `README.md` | Documentation and setup instructions. |

---

## ⚙️ Installation Requirements

Before running either script, ensure the following software and tools are installed on your system:

| Tool | Purpose | Installation Command (Ubuntu/Debian) |
|------|----------|--------------------------------------|
| **STAR** | RNA-seq alignment and index building | `sudo apt install star` |
| **wget** | File download utility | `sudo apt install wget` |
| **gunzip** | File decompression | `sudo apt install gzip` |
| **kraken2** | Taxonomic classification | `sudo apt install kraken2` |
| **awk** | Text file processing | `sudo apt install gawk` |
| **bash** | Shell interpreter (default) | *(pre-installed on most systems)* |

> 💡 Tip: For cluster or HPC environments, load these via your module system, e.g. `module load star kraken2`.

---

## 🚀 How to Use the Scripts

### **Option 1 — Full Setup**

**Recommended for production pipelines, servers, or large-scale analyses.**

```bash
chmod +x setup_references.sh
./setup_references.sh
```

### **Option 2 - Minimal Setup**
```bash
chmod +x setup_references_light.sh
./setup_references_light.sh
```