## Build and Run

### 1. Build the image
```bash
docker build -t bioinf:latest .
```

### 2. Run interactively
```bash
docker run -it --rm bioinf:latest
```

**Inside the container:**
```bash
(bioinf) bash-4.4$ fastqc --version
(bioinf) bash-4.4$ samtools view --help
```

---

### 🧱 Optional Improvements

If you want a lighter base image:
- Replace `ubuntu:22.04` with `condaforge/miniforge3`.
- Then you can skip the Miniconda install section.
