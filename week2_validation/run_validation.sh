#!/bin/bash
# Week 2 Validation - Simple Runner for Linux/Mac
# Usage: ./run_validation.sh path/to/fusion_data.csv output_name
# Example: ./run_validation.sh path/to/my_fusion.parquet my_output

echo "================================================================================"
echo "  Week 2 Validation - Complete Pipeline Runner"
echo "================================================================================"
echo ""

python week2_validation/run_complete.py "$@"

echo ""
echo "Done!"
