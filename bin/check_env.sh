#!/usr/bin/env bash
# Simple helper script to print environment versions for debugging
set -euo pipefail

echo "STAR: $(STAR --version 2>/dev/null || echo 'not found')"
echo "Arriba: $(arriba -h 2>&1 | head -n1 || echo 'not found')"
echo "STAR-Fusion: $(STAR-Fusion --version 2>/dev/null || echo 'not found')"
echo "Kraken2: $(kraken2 --version 2>/dev/null || echo 'not found')"
