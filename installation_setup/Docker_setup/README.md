## Pull, Build and Run

### 1. Pull from Docker Hub (Recommended)
```bash
docker pull your-dockerhub-username/bioinf:latest
```

### 2. Run interactively
```bash
docker run -it --rm your-dockerhub-username/bioinf:latest
```

### 3. Build the image locally (Alternative)
```bash
docker build -t your-dockerhub-username/bioinf:latest .
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
