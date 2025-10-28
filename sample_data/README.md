# Get your Sample Data Set for Test

<b>First Install</b>  SRA toolkit
### Installation Guide

---

## 1. Linux (Ubuntu/Debian/CentOS/etc.)

### Option 1 — via conda (recommended if you use Conda)
```bash
conda install -c bioconda sra-tools
```
That’s it — now you can run:
```bash
fastq-dump --split-files <your SRA Accession Number>
```

### Option 2 — manual installation
1. Go to the [SRA Toolkit download page](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/)
2. Download the latest Linux tarball:
   ```bash
   wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
   ```
3. Extract and move to `/usr/local`:
   ```bash
   tar -xvzf sratoolkit.current-ubuntu64.tar.gz
   sudo mv sratoolkit.* /usr/local/sratoolkit
   ```
4. Add it to your PATH:
   ```bash
   echo 'export PATH=$PATH:/usr/local/sratoolkit/bin' >> ~/.bashrc
   source ~/.bashrc
   ```
5. Test:
   ```bash
   fastq-dump --version
   ```

---

### 🍏 2. macOS

#### Option 1 — via Homebrew
```bash
brew install sratoolkit
```

#### Option 2 — via conda
```bash
conda install -c bioconda sra-tools
```

#### Option 3 — manual
Download from the [official NCBI site](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/), extract, and add to PATH just like in Linux.

---

### 🪟 3. Windows

1. Go to:
   [🔗 SRA Toolkit Download](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
2. Download the Windows 64-bit installer (.zip)
3. Unzip it somewhere (e.g. `C:\sratoolkit`)
4. Add the path `C:\sratoolkit\bin` to your System PATH
   *(Control Panel → System → Advanced System Settings → Environment Variables)*
5. Open a new terminal (cmd or PowerShell) and test:
   ```bash
   fastq-dump --version
   ```

---

### ✅ Test Human Samples

To download paired-end FASTQ files from SRA:
```bash
fastq-dump --split-files DRR791410
```
If you want to store output in a specific directory:
```bash
fastq-dump --split-files --outdir ./data DRR791410
```


### Other Test Accession Numbers
- DRR791410
- DRR791411
- DRR791412
- DRR791413
- DRR791414