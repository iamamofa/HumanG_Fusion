```markdown
## How to run

- **Make the script executable:**
  ```bash
  chmod +x install.sh
  ```

- **Run it:**
  ```bash
  ./install.sh
  ```

- **After the script finishes, activate the conda env (preferred path) with:**
  ```bash
  # if you installed Miniconda at $HOME/miniconda3
  source $HOME/miniconda3/etc/profile.d/conda.sh
  conda activate bioinf
  ```

- **Or (fallback venv):**
  ```bash
  source ~/pyenv-bioinf/bin/activate
  ```

---

### Important notes & caveats (please read)

- **STAR-Fusion** and some other heavy tools may require additional C/C++ libraries or specific versions and can be large. The conda install usually works but sometimes needs extra dependencies (or longer disk space). If a particular package fails in conda, re-run the conda install command manually to see the error and decide whether to pin a different version.
- **star-fusion** often requires **CTAT resource libraries** (large) that are not installed by default. The script does **not** fetch CTAT libs — you should download those separately ([STAR-Fusion docs](https://github.com/STAR-Fusion/STAR-Fusion/wiki)).
- If your system has limited RAM/CPU or disk, the conda install of many bio packages may fail. In that case, using a **Docker image** (your Dockerfile snippet) is more reproducible.
```
