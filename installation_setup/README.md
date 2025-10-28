## How to use it

```bash
chmod +x install.py
./install.py
```

### What it does:
- Installs required apt packages (if you have sudo)
- Installs Miniconda under `~/miniconda3` (if missing)
- Creates a `bioinf` environment via mamba
- Installs all conda + pip packages
- If any conda step fails → auto-creates a Python venv as fallback

---

### ✅ After installation

**Activate your environment:**

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bioinf
```

**Or (if fallback was used):**

```bash
source ~/pyenv-bioinf/bin/activate
```

