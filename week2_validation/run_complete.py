#!/usr/bin/env python3
"""
Week 2 Validation: Complete End-to-End Runner

Single-command execution of the entire Week 2 validation pipeline including:
- Full diagnostic pipeline (distribution, Benford, log-normality, COSMIC)
- Complete test suite
- Output verification

Usage:
    python week2_validation/run_complete.py [fusion_data] [output_dir]
    
Examples:
    # Use demo data
    python week2_validation/run_complete.py
    
    # Use custom data
    python week2_validation/run_complete.py path/to/data.csv output/
    
    # With COSMIC data
    python week2_validation/run_complete.py data.csv output/ --cosmic cosmic.tsv
"""

import sys
import subprocess
import argparse
from pathlib import Path


def print_banner(text, char="="):
    """Print a formatted banner."""
    print(f"\n{char * 80}")
    print(f"  {text}")
    print(f"{char * 80}\n")


def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"🚀 {description}...")
    print(f"   Command: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,
            text=True,
        )
        print(f"✅ {description} completed successfully\n")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with exit code {e.returncode}\n")
        return False
    except Exception as e:
        print(f"❌ Error during {description}: {e}\n")
        return False


def verify_outputs(output_dir, dataset_stem):
    """Verify expected output files exist."""
    print_banner("Verifying Outputs", "~")
    
    output_path = Path(output_dir)
    expected_files = [
        f"week2_status_{dataset_stem}.json",
        f"week2_dataset_certification_{dataset_stem}.json",
        f"week2_diagnostic_results_{dataset_stem}.json",
        f"week2_cleaned_dataset_{dataset_stem}.csv",
        f"Statistical_Integrity_Report_{dataset_stem}.md",
        "protein_distribution.png",
        "skewness_diagnostic.png",
    ]
    
    all_exist = True
    for filename in expected_files:
        filepath = output_path / filename
        if filepath.exists():
            size = filepath.stat().st_size
            print(f"  ✅ {filename} ({size:,} bytes)")
        else:
            print(f"  ❌ {filename} (missing)")
            all_exist = False
    
    print()
    return all_exist


def verify_new_features(output_dir, dataset_stem):
    """Verify new scientific defensibility features are present."""
    print_banner("Verifying New Features", "~")
    
    import json
    
    diagnostic_file = Path(output_dir) / f"week2_diagnostic_results_{dataset_stem}.json"
    
    if not diagnostic_file.exists():
        print("  ❌ Diagnostic results file not found")
        return False
    
    try:
        with open(diagnostic_file, 'r') as f:
            data = json.load(f)
        
        cosmic = data.get('cosmic', {})
        
        # Check bootstrap CI
        if 'rho_ci_lower' in cosmic and 'rho_ci_upper' in cosmic:
            print(f"  ✅ Bootstrap CI: [{cosmic.get('rho_ci_lower')}, {cosmic.get('rho_ci_upper')}]")
        else:
            print("  ⚠️  Bootstrap CI not found")
        
        # Check reproducibility lock
        if 'reproducibility_lock' in cosmic:
            lock = cosmic['reproducibility_lock']
            print(f"  ✅ Reproducibility Lock:")
            print(f"     - Python: {lock.get('python_version', 'N/A')}")
            print(f"     - NumPy: {lock.get('numpy_version', 'N/A')}")
            print(f"     - Code version: {lock.get('cosmic_validation_code_version', 'N/A')}")
        else:
            print("  ⚠️  Reproducibility lock not found")
        
        # Check score breakdown
        if 'score_component_breakdown' in cosmic:
            breakdown = cosmic['score_component_breakdown']
            print(f"  ✅ Score Component Breakdown:")
            print(f"     - Correlation: {breakdown.get('correlation_component', 'N/A')}")
            print(f"     - Enrichment: {breakdown.get('enrichment_component', 'N/A')}")
            print(f"     - Negative Control: {breakdown.get('negative_control_component', 'N/A')}")
            print(f"     - Overlap: {breakdown.get('overlap_component', 'N/A')}")
        else:
            print("  ⚠️  Score breakdown not found")
        
        print()
        return True
        
    except Exception as e:
        print(f"  ❌ Error reading diagnostic file: {e}\n")
        return False


def main():
    # Fix encoding issues on Windows
    if sys.platform == 'win32':
        try:
            import io
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')
        except Exception:
            pass  # Fallback to default encoding
    
    parser = argparse.ArgumentParser(
        description="Run complete Week 2 validation pipeline with all diagnostics and tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use demo data (default)
  python week2_validation/run_complete.py
  
  # Use custom data
  python week2_validation/run_complete.py path/to/data.csv output/
  
  # With COSMIC reference data
  python week2_validation/run_complete.py data.csv output/ --cosmic cosmic.tsv
  
  # Skip tests
  python week2_validation/run_complete.py --skip-tests
        """
    )
    
    parser.add_argument(
        'fusion_data',
        nargs='?',
        default='week2_validation/demo_fusion.csv',
        help='Path to fusion data file (default: week2_validation/demo_fusion.csv)'
    )
    
    parser.add_argument(
        'output_dir',
        nargs='?',
        default='demo_output',
        help='Output directory (default: demo_output, will be placed in week2_validation/results/)'
    )
    
    parser.add_argument(
        '--cosmic',
        dest='cosmic_data',
        help='Path to COSMIC reference data (optional, uses mock if not provided)'
    )
    
    parser.add_argument(
        '--skip-tests',
        action='store_true',
        help='Skip running the test suite'
    )
    
    parser.add_argument(
        '--skip-verification',
        action='store_true',
        help='Skip output verification'
    )
    
    args = parser.parse_args()
    
    # Place output inside week2_validation/results/
    week2_validation_dir = Path(__file__).parent  # week2_validation directory
    results_dir = week2_validation_dir / "results"
    
    # Create results directory if it doesn't exist
    results_dir.mkdir(exist_ok=True)
    
    # Resolve the final output path
    output_path = results_dir / args.output_dir
    output_path_str = str(output_path)
    
    # Print welcome banner
    print_banner("Week 2 Validation - Complete Pipeline Runner", "=")
    print(f"📂 Fusion data: {args.fusion_data}")
    print(f"📁 Output directory: {output_path_str}")
    if args.cosmic_data:
        print(f"🧬 COSMIC data: {args.cosmic_data}")
    else:
        print(f"🧬 COSMIC data: Using mock fallback")
    print()
    
    # Determine dataset stem for output verification
    dataset_stem = Path(args.fusion_data).stem
    
    # Step 1: Run validation pipeline
    print_banner("Step 1: Running Validation Pipeline", "-")
    
    pipeline_cmd = [
        sys.executable,
        '-m', 'week2_validation.run_week2',
        '--fusion-data', args.fusion_data,
        '--output-dir', output_path_str,
        '--run-all',
    ]
    
    if args.cosmic_data:
        pipeline_cmd.extend(['--cosmic-data', args.cosmic_data])
    
    if not run_command(pipeline_cmd, "Validation pipeline"):
        print("⚠️  Pipeline failed. Continuing to tests anyway...\n")
    
    # Step 2: Run tests (unless skipped)
    if not args.skip_tests:
        print_banner("Step 2: Running Test Suite", "-")
        
        test_cmd = [
            sys.executable,
            '-m', 'pytest',
            'week2_validation/tests/',
            '-v',
            '--tb=short',
        ]
        
        if not run_command(test_cmd, "Test suite"):
            print("⚠️  Some tests failed. Check output above.\n")
    else:
        print_banner("Step 2: Skipping Test Suite (--skip-tests)", "-")
    
    # Step 3: Verify outputs (unless skipped)
    if not args.skip_verification:
        outputs_ok = verify_outputs(output_path_str, dataset_stem)
        features_ok = verify_new_features(output_path_str, dataset_stem)
        
        # Final summary
        print_banner("Final Summary", "=")
        
        if outputs_ok and features_ok:
            print("🎉 SUCCESS! Complete validation finished.")
            print(f"📊 All outputs generated in: {output_path_str}/")
            print(f"📄 View report: {output_path_str}/Statistical_Integrity_Report_{dataset_stem}.md")
            print()
            return 0
        else:
            print("⚠️  Validation completed with warnings.")
            print("   Some outputs or features may be missing.")
            print(f"   Check: {output_path_str}/")
            print()
            return 1
    else:
        print_banner("Verification Skipped (--skip-verification)", "=")
        print("✅ Pipeline and tests completed.")
        print(f"📊 Outputs should be in: {args.output_dir}/")
        print()
        return 0


if __name__ == '__main__':
    sys.exit(main())
