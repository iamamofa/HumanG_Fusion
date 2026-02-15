@echo off
REM Week 2 Validation - Simple Runner for Windows
REM Usage: run_validation.bat path\to\fusion_data.csv output_name
REM Example: run_validation.bat path\to\my_fusion.parquet my_output

echo ================================================================================
echo   Week 2 Validation - Complete Pipeline Runner
echo ================================================================================
echo.

python week2_validation\run_complete.py %*

echo.
echo Press any key to exit...
pause > nul
