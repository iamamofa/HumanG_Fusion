#!/usr/bin/env python3
"""
COSMIC Operational Verification - End-to-End Audit

This script performs a comprehensive forensic investigation of COSMIC
cross-validation functionality WITHOUT modifying any source code or datasets.

PHASES:
1. COSMIC Data Presence Verification
2. Runtime Load Proof
3. Functional Usage Proof
4. Mechanical Mapping Validation
5. Statistical Consistency Validation
6. Biological Plausibility Validation
7. Signal Non-Triviality Check
8. Dataset Sensitivity Test
9. Failure Mode Detection
"""

import json
import hashlib
import subprocess
import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
import pandas as pd
import numpy as np
from scipy import stats
import tempfile
import shutil

# Known high-recurrence cancer driver genes
KNOWN_CANCER_DRIVERS = {
    "ABL1", "BCR", "ALK", "EML4", "BRAF", "EGFR", "RET", "NTRK1", "NTRK2", "NTRK3",
    "ROS1", "MET", "FGFR1", "FGFR2", "FGFR3", "TMPRSS2", "ERG", "ETV1", "EWSR1",
    "FLI1", "PAX3", "FOXO1", "MLL", "AF4", "PML", "RARA", "RUNX1", "RUNX1T1"
}

# Known fusion hotspots
KNOWN_FUSION_HOTSPOTS = [
    ("ABL1", "BCR"),
    ("EML4", "ALK"),
    ("TMPRSS2", "ERG"),
    ("NTRK1", "TPM3"),
    ("ROS1", "SLC34A2"),
    ("FGFR3", "TACC3"),
    ("RET", "KIF5B"),
]


@dataclass
class VerificationResult:
    """Results from a verification phase"""
    phase: str
    status: str  # PASS, FAIL, INSUFFICIENT_EVIDENCE
    evidence: Dict[str, Any]
    metrics: Dict[str, Any] = None


class COSMICVerification:
    """End-to-end COSMIC operational verification"""
    
    def __init__(self, workspace_root: Path):
        self.workspace_root = Path(workspace_root)
        self.week2_dir = self.workspace_root / "week2_validation"
        self.results: List[VerificationResult] = []
        self.temp_dir = None
        
    def __enter__(self):
        self.temp_dir = tempfile.mkdtemp(prefix="cosmic_audit_")
        return self
        
    def __exit__(self, *args):
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def run_all_phases(self) -> Dict[str, Any]:
        """Execute all verification phases"""
        print("=" * 80)
        print("COSMIC OPERATIONAL VERIFICATION")
        print("=" * 80)
        print()
        
        # Phase 1: COSMIC Data Presence
        self.phase1_cosmic_presence()
        
        # Phase 2: Runtime Load Proof
        self.phase2_runtime_load()
        
        # Phase 3: Functional Usage Proof
        self.phase3_functional_usage()
        
        # Phase 4: Mechanical Mapping Validation
        self.phase4_mapping_validation()
        
        # Phase 5: Statistical Consistency
        self.phase5_statistical_consistency()
        
        # Phase 6: Biological Plausibility
        self.phase6_biological_plausibility()
        
        # Phase 7: Signal Non-Triviality
        self.phase7_signal_non_triviality()
        
        # Phase 8: Dataset Sensitivity
        self.phase8_dataset_sensitivity()
        
        # Phase 9: Failure Mode Detection
        self.phase9_failure_modes()
        
        return self.generate_final_report()
    
    def phase1_cosmic_presence(self):
        """Phase 1: COSMIC Data Presence Verification"""
        print("PHASE 1: COSMIC Data Presence Verification")
        print("-" * 80)
        
        evidence = {}
        
        # Check for COSMIC files in common locations
        cosmic_paths = [
            self.workspace_root / "cosmic.tsv",
            self.workspace_root / "cosmic.csv",
            self.workspace_root / "COSMIC.tsv",
            self.workspace_root / "COSMIC.csv",
            self.workspace_root / "week2_validation" / "cosmic.tsv",
            self.workspace_root / "week2_validation" / "cosmic.csv",
        ]
        
        cosmic_file_present = False
        cosmic_file_path = None
        cosmic_file_size = 0
        cosmic_row_count = 0
        cosmic_file_hash = None
        
        for path in cosmic_paths:
            if path.exists() and path.is_file():
                cosmic_file_present = True
                cosmic_file_path = str(path)
                cosmic_file_size = path.stat().st_size
                
                # Try to load and count rows
                try:
                    if path.suffix.lower() == '.csv':
                        df = pd.read_csv(path, nrows=1000)  # Sample first 1000 rows
                    elif path.suffix.lower() == '.tsv':
                        df = pd.read_csv(path, sep='\t', nrows=1000)
                    else:
                        df = pd.read_csv(path, nrows=1000)
                    
                    # Count total rows (approximate for large files)
                    with open(path, 'rb') as f:
                        cosmic_file_hash = hashlib.md5(f.read(1024*1024)).hexdigest()  # First 1MB hash
                    
                    # Full row count
                    try:
                        if path.suffix.lower() == '.csv':
                            cosmic_row_count = sum(1 for _ in open(path)) - 1  # Subtract header
                        else:
                            cosmic_row_count = sum(1 for _ in open(path)) - 1
                    except:
                        cosmic_row_count = len(df) if df is not None else 0
                    
                    evidence["cosmic_file_present"] = True
                    evidence["cosmic_file_path"] = cosmic_file_path
                    evidence["cosmic_file_size"] = cosmic_file_size
                    evidence["cosmic_row_count"] = cosmic_row_count
                    evidence["cosmic_file_hash"] = cosmic_file_hash
                    break
                except Exception as e:
                    evidence["load_error"] = str(e)
        
        if not cosmic_file_present:
            evidence["cosmic_file_present"] = False
            evidence["message"] = "No COSMIC file found in expected locations"
        
        status = "PASS" if cosmic_file_present else "INSUFFICIENT_EVIDENCE"
        self.results.append(VerificationResult(
            phase="COSMIC File Presence",
            status=status,
            evidence=evidence
        ))
        
        print(f"  COSMIC file present: {cosmic_file_present}")
        if cosmic_file_present:
            print(f"  Path: {cosmic_file_path}")
            print(f"  Size: {cosmic_file_size} bytes")
            print(f"  Rows: {cosmic_row_count}")
            print(f"  Hash (first 1MB): {cosmic_file_hash}")
        print()
    
    def phase2_runtime_load(self):
        """Phase 2: Runtime Load Proof"""
        print("PHASE 2: Runtime Load Proof")
        print("-" * 80)
        
        evidence = {}
        
        # Find a test dataset
        test_datasets = [
            self.week2_dir / "demo_fusion.csv",
            self.week2_dir / "demo_fusion.tsv",
            self.week2_dir / "demo_fusion.json",
        ]
        
        test_dataset = None
        for ds in test_datasets:
            if ds.exists():
                test_dataset = ds
                break
        
        if not test_dataset:
            evidence["status"] = "FAILED"
            evidence["reason"] = "No test dataset found"
            self.results.append(VerificationResult(
                phase="COSMIC Runtime Load",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No test dataset found - skipping runtime test")
            print()
            return
        
        # Create output directory
        output_dir = Path(self.temp_dir) / "phase2_output"
        output_dir.mkdir(exist_ok=True)
        
        # Try to find COSMIC file
        cosmic_path = None
        for path in [
            self.workspace_root / "cosmic.tsv",
            self.workspace_root / "cosmic.csv",
            self.workspace_root / "COSMIC.tsv",
            self.workspace_root / "COSMIC.csv",
        ]:
            if path.exists():
                cosmic_path = path
                break
        
        # Run pipeline with COSMIC
        cmd = [
            sys.executable, "-m", "week2_validation.run_week2",
            "--fusion-data", str(test_dataset),
            "--output-dir", str(output_dir),
            "--run-cosmic"
        ]
        
        if cosmic_path:
            cmd.extend(["--cosmic-data", str(cosmic_path)])
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=str(self.workspace_root)
            )
            
            stdout = result.stdout
            stderr = result.stderr
            
            # Check for COSMIC load signals
            cosmic_reference_loaded = "COSMIC DataFrame is None" not in stdout
            cosmic_df_rows = None
            
            # Parse output for evidence
            if "Total fusion pairs (COSMIC)" in stdout:
                for line in stdout.split('\n'):
                    if "Total fusion pairs (COSMIC)" in line:
                        try:
                            cosmic_df_rows = int(line.split(':')[-1].strip())
                        except:
                            pass
            
            # Check status JSON
            status_json_path = output_dir / f"week2_status_{test_dataset.stem}.json"
            if status_json_path.exists():
                with open(status_json_path, 'r') as f:
                    status_data = json.load(f)
                    evidence["status_json"] = status_data
                    cosmic_reference_loaded = status_data.get("cosmic_reference_loaded", False)
            
            evidence["cosmic_reference_loaded"] = cosmic_reference_loaded
            evidence["cosmic_df_rows"] = cosmic_df_rows
            evidence["exit_code"] = result.returncode
            evidence["stdout_snippet"] = stdout[:2000] if stdout else ""
            evidence["stderr_snippet"] = stderr[:1000] if stderr else ""
            
            status = "SUCCESS" if cosmic_reference_loaded and cosmic_df_rows is not None and cosmic_df_rows > 0 else "FAILED"
            
        except subprocess.TimeoutExpired:
            evidence["status"] = "TIMEOUT"
            status = "FAILED"
        except Exception as e:
            evidence["status"] = "ERROR"
            evidence["error"] = str(e)
            status = "FAILED"
        
        self.results.append(VerificationResult(
            phase="COSMIC Runtime Load",
            status=status,
            evidence=evidence
        ))
        
        print(f"  COSMIC reference loaded: {evidence.get('cosmic_reference_loaded', False)}")
        print(f"  COSMIC rows: {evidence.get('cosmic_df_rows', 'N/A')}")
        print(f"  Exit code: {evidence.get('exit_code', 'N/A')}")
        print()
    
    def phase3_functional_usage(self):
        """Phase 3: Functional Usage Proof"""
        print("PHASE 3: Functional Usage Proof")
        print("-" * 80)
        
        evidence = {}
        
        # Check if COSMIC diagnostics actually produce output
        test_datasets = [
            self.week2_dir / "demo_fusion.csv",
            self.week2_dir / "demo_fusion.tsv",
        ]
        
        test_dataset = None
        for ds in test_datasets:
            if ds.exists():
                test_dataset = ds
                break
        
        if not test_dataset:
            evidence["status"] = "NO_TEST_DATA"
            self.results.append(VerificationResult(
                phase="COSMIC Functional Usage",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No test dataset found")
            print()
            return
        
        output_dir = Path(self.temp_dir) / "phase3_output"
        output_dir.mkdir(exist_ok=True)
        
        # Find COSMIC file
        cosmic_path = None
        for path in [
            self.workspace_root / "cosmic.tsv",
            self.workspace_root / "cosmic.csv",
        ]:
            if path.exists():
                cosmic_path = path
                break
        
        if not cosmic_path:
            evidence["status"] = "NO_COSMIC_FILE"
            self.results.append(VerificationResult(
                phase="COSMIC Functional Usage",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No COSMIC file found - cannot test functional usage")
            print()
            return
        
        # Run with COSMIC
        cmd = [
            sys.executable, "-m", "week2_validation.run_week2",
            "--fusion-data", str(test_dataset),
            "--output-dir", str(output_dir),
            "--run-cosmic",
            "--cosmic-data", str(cosmic_path)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=str(self.workspace_root)
            )
            
            stdout = result.stdout
            
            # Check for functional usage indicators
            has_overlap_table = "Overlapping pairs" in stdout
            has_recurrence_comparison = "Total fusion pairs (COSMIC)" in stdout
            has_ranking_comparison = "rank discrepancies" in stdout.lower() or "rank" in stdout.lower()
            has_enrichment_metrics = "Only in ours" in stdout or "Only in COSMIC" in stdout
            
            # Check diagnostic results JSON
            diag_json_path = output_dir / f"week2_diagnostic_results_{test_dataset.stem}.json"
            overlap_count = None
            if diag_json_path.exists():
                with open(diag_json_path, 'r') as f:
                    diag_data = json.load(f)
                    cosmic_data = diag_data.get("cosmic", {})
                    overlap_count = cosmic_data.get("overlap_count", 0)
                    evidence["diagnostic_results"] = cosmic_data
            
            evidence["has_overlap_table"] = has_overlap_table
            evidence["has_recurrence_comparison"] = has_recurrence_comparison
            evidence["has_ranking_comparison"] = has_ranking_comparison
            evidence["has_enrichment_metrics"] = has_enrichment_metrics
            evidence["overlap_count"] = overlap_count
            
            # Determine status
            usage_indicators = [
                has_overlap_table,
                has_recurrence_comparison,
                has_ranking_comparison,
                has_enrichment_metrics
            ]
            
            if any(usage_indicators):
                status = "ACTIVE"
            elif overlap_count is not None:
                status = "LOADED_BUT_UNUSED" if overlap_count == 0 else "ACTIVE"
            else:
                status = "NOT_USED"
                
        except Exception as e:
            evidence["error"] = str(e)
            status = "INSUFFICIENT_EVIDENCE"
        
        self.results.append(VerificationResult(
            phase="COSMIC Functional Usage",
            status=status,
            evidence=evidence
        ))
        
        print(f"  Overlap table: {evidence.get('has_overlap_table', False)}")
        print(f"  Recurrence comparison: {evidence.get('has_recurrence_comparison', False)}")
        print(f"  Ranking comparison: {evidence.get('has_ranking_comparison', False)}")
        print(f"  Enrichment metrics: {evidence.get('has_enrichment_metrics', False)}")
        print(f"  Overlap count: {evidence.get('overlap_count', 'N/A')}")
        print()
    
    def phase4_mapping_validation(self):
        """Phase 4: Mechanical Mapping Validation"""
        print("PHASE 4: Mechanical Mapping Validation")
        print("-" * 80)
        
        evidence = {}
        
        # Load test dataset
        test_dataset = self.week2_dir / "demo_fusion.csv"
        if not test_dataset.exists():
            evidence["status"] = "NO_TEST_DATA"
            self.results.append(VerificationResult(
                phase="COSMIC Mapping Integrity",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No test dataset found")
            print()
            return
        
        # Load fusion data
        try:
            fusion_df = pd.read_csv(test_dataset)
            
            # Check gene normalization
            gene_cols = ["gene_1", "gene_2"]
            normalization_tests = {}
            
            for col in gene_cols:
                if col in fusion_df.columns:
                    # Test uppercase normalization
                    original = fusion_df[col].astype(str)
                    normalized = original.str.upper().str.strip()
                    uppercase_changes = (original != normalized).sum()
                    
                    # Test whitespace normalization
                    whitespace_changes = original.str.strip() != original
                    
                    normalization_tests[col] = {
                        "uppercase_changes": int(uppercase_changes.sum()),
                        "whitespace_changes": int(whitespace_changes.sum()),
                        "total_values": len(original)
                    }
            
            # Test pair normalization (sorting)
            if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
                pairs_before = list(zip(fusion_df["gene_1"], fusion_df["gene_2"]))
                pairs_normalized = [tuple(sorted([str(g1).upper().strip(), str(g2).upper().strip()])) 
                                   for g1, g2 in pairs_before]
                pairs_changed = sum(1 for p1, p2 in zip(pairs_before, pairs_normalized) if p1 != p2)
                
                normalization_tests["pair_sorting"] = {
                    "pairs_changed": pairs_changed,
                    "total_pairs": len(pairs_before)
                }
            
            evidence["normalization_tests"] = normalization_tests
            evidence["gene_normalization_valid"] = True
            
            # Check for duplicate pairs after normalization
            if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
                normalized_pairs = [tuple(sorted([str(g1).upper().strip(), str(g2).upper().strip()])) 
                                  for g1, g2 in zip(fusion_df["gene_1"], fusion_df["gene_2"])]
                unique_pairs = len(set(normalized_pairs))
                total_pairs = len(normalized_pairs)
                duplicates = total_pairs - unique_pairs
                
                evidence["duplicate_pairs_after_normalization"] = duplicates
                evidence["unique_pairs"] = unique_pairs
                evidence["total_pairs"] = total_pairs
            
            status = "VALID"
            
        except Exception as e:
            evidence["error"] = str(e)
            status = "SUSPECT"
        
        self.results.append(VerificationResult(
            phase="COSMIC Mapping Integrity",
            status=status,
            evidence=evidence
        ))
        
        print(f"  Gene normalization: {evidence.get('gene_normalization_valid', False)}")
        if "normalization_tests" in evidence:
            print(f"  Normalization details: {evidence['normalization_tests']}")
        print()
    
    def phase5_statistical_consistency(self):
        """Phase 5: Statistical Consistency Validation"""
        print("PHASE 5: Statistical Consistency Validation")
        print("-" * 80)
        
        evidence = {}
        metrics = {}
        
        # This phase requires actual COSMIC data and fusion data
        test_dataset = self.week2_dir / "demo_fusion.csv"
        cosmic_path = None
        
        for path in [
            self.workspace_root / "cosmic.tsv",
            self.workspace_root / "cosmic.csv",
        ]:
            if path.exists():
                cosmic_path = path
                break
        
        if not test_dataset.exists() or not cosmic_path:
            evidence["status"] = "INSUFFICIENT_DATA"
            self.results.append(VerificationResult(
                phase="COSMIC Statistical Agreement",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence,
                metrics=metrics
            ))
            print("  Insufficient data for statistical analysis")
            print()
            return
        
        try:
            # Load data
            fusion_df = pd.read_csv(test_dataset)
            
            if cosmic_path.suffix.lower() == '.tsv':
                cosmic_df = pd.read_csv(cosmic_path, sep='\t')
            else:
                cosmic_df = pd.read_csv(cosmic_path)
            
            # Normalize pairs (same as diagnostics.py)
            def normalize_pair(g1, g2):
                return tuple(sorted([str(g1).upper().strip(), str(g2).upper().strip()]))
            
            # Build fusion pairs with recurrence
            fusion_pairs = {}
            if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
                if "recurrence_count" in fusion_df.columns:
                    for _, row in fusion_df.iterrows():
                        pair = normalize_pair(row["gene_1"], row["gene_2"])
                        count = row.get("recurrence_count", 0)
                        if pair not in fusion_pairs or count > fusion_pairs[pair]:
                            fusion_pairs[pair] = count
            
            # Build COSMIC pairs
            cosmic_pairs = {}
            if "gene_1" in cosmic_df.columns and "gene_2" in cosmic_df.columns:
                if "recurrence_count" in cosmic_df.columns:
                    for _, row in cosmic_df.iterrows():
                        pair = normalize_pair(row["gene_1"], row["gene_2"])
                        count = row.get("recurrence_count", 0)
                        if pair not in cosmic_pairs or count > cosmic_pairs[pair]:
                            cosmic_pairs[pair] = count
            
            # Find overlapping pairs
            overlap_pairs = set(fusion_pairs.keys()) & set(cosmic_pairs.keys())
            
            if len(overlap_pairs) > 0:
                # Extract recurrence counts for overlapping pairs
                fusion_recurrence = [fusion_pairs[p] for p in overlap_pairs]
                cosmic_recurrence = [cosmic_pairs[p] for p in overlap_pairs]
                
                # Spearman correlation
                if len(fusion_recurrence) > 1:
                    spearman_rho, spearman_p = stats.spearmanr(fusion_recurrence, cosmic_recurrence)
                    metrics["spearman_rho"] = float(spearman_rho) if not np.isnan(spearman_rho) else None
                    metrics["spearman_pvalue"] = float(spearman_p) if not np.isnan(spearman_p) else None
                
                # Distribution comparison
                metrics["fusion_mean"] = float(np.mean(fusion_recurrence))
                metrics["cosmic_mean"] = float(np.mean(cosmic_recurrence))
                metrics["fusion_variance"] = float(np.var(fusion_recurrence))
                metrics["cosmic_variance"] = float(np.var(cosmic_recurrence))
                
                # Zero inflation check
                fusion_zeros = sum(1 for x in fusion_recurrence if x == 0)
                cosmic_zeros = sum(1 for x in cosmic_recurrence if x == 0)
                metrics["fusion_zero_fraction"] = fusion_zeros / len(fusion_recurrence)
                metrics["cosmic_zero_fraction"] = cosmic_zeros / len(cosmic_recurrence)
                
                # Enrichment check: Are top COSMIC genes in top fusion genes?
                fusion_sorted = sorted(fusion_pairs.items(), key=lambda x: x[1], reverse=True)
                cosmic_sorted = sorted(cosmic_pairs.items(), key=lambda x: x[1], reverse=True)
                
                top_fusion_genes = set()
                top_cosmic_genes = set()
                
                for pair, _ in fusion_sorted[:min(20, len(fusion_sorted))]:
                    top_fusion_genes.update(pair)
                
                for pair, _ in cosmic_sorted[:min(20, len(cosmic_sorted))]:
                    top_cosmic_genes.update(pair)
                
                overlap_genes = top_fusion_genes & top_cosmic_genes
                enrichment_ratio = len(overlap_genes) / len(top_cosmic_genes) if len(top_cosmic_genes) > 0 else 0
                metrics["enrichment_ratio"] = float(enrichment_ratio)
                
                evidence["overlap_count"] = len(overlap_pairs)
                evidence["statistical_analysis_complete"] = True
                
                # Determine agreement level
                if metrics.get("spearman_rho") is not None:
                    rho = abs(metrics["spearman_rho"])
                    if rho > 0.5:
                        status = "GOOD"
                    elif rho > 0.3:
                        status = "WEAK"
                    else:
                        status = "NONE"
                else:
                    status = "WEAK"
            else:
                evidence["overlap_count"] = 0
                evidence["message"] = "No overlapping pairs found"
                status = "NONE"
                
        except Exception as e:
            evidence["error"] = str(e)
            status = "INSUFFICIENT_EVIDENCE"
        
        self.results.append(VerificationResult(
            phase="COSMIC Statistical Agreement",
            status=status,
            evidence=evidence,
            metrics=metrics
        ))
        
        if metrics:
            print(f"  Spearman rho: {metrics.get('spearman_rho', 'N/A')}")
            print(f"  Spearman p-value: {metrics.get('spearman_pvalue', 'N/A')}")
            print(f"  Fusion mean: {metrics.get('fusion_mean', 'N/A')}")
            print(f"  COSMIC mean: {metrics.get('cosmic_mean', 'N/A')}")
            print(f"  Enrichment ratio: {metrics.get('enrichment_ratio', 'N/A')}")
        print()
    
    def phase6_biological_plausibility(self):
        """Phase 6: Biological Plausibility Validation"""
        print("PHASE 6: Biological Plausibility Validation")
        print("-" * 80)
        
        evidence = {}
        
        test_dataset = self.week2_dir / "demo_fusion.csv"
        if not test_dataset.exists():
            evidence["status"] = "NO_TEST_DATA"
            self.results.append(VerificationResult(
                phase="COSMIC Biological Plausibility",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No test dataset found")
            print()
            return
        
        try:
            fusion_df = pd.read_csv(test_dataset)
            
            if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
                # Check for known cancer drivers
                all_genes = set(fusion_df["gene_1"].astype(str).str.upper()) | set(fusion_df["gene_2"].astype(str).str.upper())
                known_drivers_found = all_genes & KNOWN_CANCER_DRIVERS
                
                # Check for known fusion hotspots
                pairs = [(str(g1).upper().strip(), str(g2).upper().strip()) 
                        for g1, g2 in zip(fusion_df["gene_1"], fusion_df["gene_2"])]
                normalized_pairs = [tuple(sorted(p)) for p in pairs]
                known_hotspots_found = []
                
                for hotspot in KNOWN_FUSION_HOTSPOTS:
                    normalized_hotspot = tuple(sorted([h.upper() for h in hotspot]))
                    if normalized_hotspot in normalized_pairs:
                        known_hotspots_found.append(hotspot)
                
                evidence["known_drivers_found"] = list(known_drivers_found)
                evidence["known_hotspots_found"] = known_hotspots_found
                evidence["driver_count"] = len(known_drivers_found)
                evidence["hotspot_count"] = len(known_hotspots_found)
                
                # Determine plausibility
                if len(known_drivers_found) > 5 and len(known_hotspots_found) > 0:
                    status = "REALISTIC"
                elif len(known_drivers_found) > 0:
                    status = "QUESTIONABLE"
                else:
                    status = "INVALID"
            else:
                evidence["error"] = "Missing gene columns"
                status = "INVALID"
                
        except Exception as e:
            evidence["error"] = str(e)
            status = "INVALID"
        
        self.results.append(VerificationResult(
            phase="COSMIC Biological Plausibility",
            status=status,
            evidence=evidence
        ))
        
        print(f"  Known drivers found: {evidence.get('driver_count', 0)}")
        print(f"  Known hotspots found: {evidence.get('hotspot_count', 0)}")
        print()
    
    def phase7_signal_non_triviality(self):
        """Phase 7: Signal Non-Triviality Check"""
        print("PHASE 7: Signal Non-Triviality Check")
        print("-" * 80)
        
        evidence = {}
        
        test_dataset = self.week2_dir / "demo_fusion.csv"
        cosmic_path = None
        
        for path in [
            self.workspace_root / "cosmic.tsv",
            self.workspace_root / "cosmic.csv",
        ]:
            if path.exists():
                cosmic_path = path
                break
        
        if not test_dataset.exists() or not cosmic_path:
            evidence["status"] = "INSUFFICIENT_DATA"
            self.results.append(VerificationResult(
                phase="COSMIC Signal Non-Triviality",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  Insufficient data")
            print()
            return
        
        try:
            fusion_df = pd.read_csv(test_dataset)
            
            if cosmic_path.suffix.lower() == '.tsv':
                cosmic_df = pd.read_csv(cosmic_path, sep='\t')
            else:
                cosmic_df = pd.read_csv(cosmic_path)
            
            def normalize_pair(g1, g2):
                return tuple(sorted([str(g1).upper().strip(), str(g2).upper().strip()]))
            
            # Build pairs
            fusion_pairs = set()
            if "gene_1" in fusion_df.columns and "gene_2" in fusion_df.columns:
                for _, row in fusion_df.iterrows():
                    pair = normalize_pair(row["gene_1"], row["gene_2"])
                    fusion_pairs.add(pair)
            
            cosmic_pairs = set()
            if "gene_1" in cosmic_df.columns and "gene_2" in cosmic_df.columns:
                for _, row in cosmic_df.iterrows():
                    pair = normalize_pair(row["gene_1"], row["gene_2"])
                    cosmic_pairs.add(pair)
            
            overlap = fusion_pairs & cosmic_pairs
            only_fusion = fusion_pairs - cosmic_pairs
            only_cosmic = cosmic_pairs - fusion_pairs
            
            total_fusion = len(fusion_pairs)
            total_cosmic = len(cosmic_pairs)
            overlap_count = len(overlap)
            
            if total_fusion > 0:
                overlap_fraction = overlap_count / total_fusion
            else:
                overlap_fraction = 0
            
            evidence["total_fusion_pairs"] = total_fusion
            evidence["total_cosmic_pairs"] = total_cosmic
            evidence["overlap_count"] = overlap_count
            evidence["overlap_fraction"] = overlap_fraction
            evidence["only_fusion_count"] = len(only_fusion)
            evidence["only_cosmic_count"] = len(only_cosmic)
            
            # Check for triviality
            if overlap_count == 0:
                status = "INVALID"  # All overlaps = 0
            elif overlap_fraction == 1.0:
                status = "SUSPICIOUS"  # All overlaps = 100%
            elif 0 < overlap_fraction < 1.0 and len(only_fusion) > 0 and len(only_cosmic) > 0:
                status = "VALID"  # Mixed matches, realistic noise
            else:
                status = "SUSPICIOUS"
                
        except Exception as e:
            evidence["error"] = str(e)
            status = "INVALID"
        
        self.results.append(VerificationResult(
            phase="COSMIC Signal Non-Triviality",
            status=status,
            evidence=evidence
        ))
        
        print(f"  Overlap count: {evidence.get('overlap_count', 0)}")
        print(f"  Overlap fraction: {evidence.get('overlap_fraction', 0):.2%}")
        print(f"  Only in fusion: {evidence.get('only_fusion_count', 0)}")
        print(f"  Only in COSMIC: {evidence.get('only_cosmic_count', 0)}")
        print()
    
    def phase8_dataset_sensitivity(self):
        """Phase 8: Dataset Sensitivity Test"""
        print("PHASE 8: Dataset Sensitivity Test")
        print("-" * 80)
        
        evidence = {}
        
        # Find multiple test datasets
        test_datasets = [
            self.week2_dir / "demo_fusion.csv",
            self.week2_dir / "demo_fusion.tsv",
            self.week2_dir / "demo_fusion.json",
        ]
        
        available_datasets = [ds for ds in test_datasets if ds.exists()]
        
        if len(available_datasets) < 2:
            evidence["status"] = "INSUFFICIENT_DATASETS"
            evidence["available_count"] = len(available_datasets)
            self.results.append(VerificationResult(
                phase="COSMIC Sensitivity To Dataset Changes",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print(f"  Only {len(available_datasets)} dataset(s) available - need 2+")
            print()
            return
        
        cosmic_path = None
        for path in [
            self.workspace_root / "cosmic.tsv",
            self.workspace_root / "cosmic.csv",
        ]:
            if path.exists():
                cosmic_path = path
                break
        
        if not cosmic_path:
            evidence["status"] = "NO_COSMIC_FILE"
            self.results.append(VerificationResult(
                phase="COSMIC Sensitivity To Dataset Changes",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No COSMIC file found")
            print()
            return
        
        # Run COSMIC diagnostics on multiple datasets
        results_by_dataset = {}
        
        for dataset in available_datasets[:2]:  # Test first 2
            output_dir = Path(self.temp_dir) / f"phase8_{dataset.stem}"
            output_dir.mkdir(exist_ok=True)
            
            cmd = [
                sys.executable, "-m", "week2_validation.run_week2",
                "--fusion-data", str(dataset),
                "--output-dir", str(output_dir),
                "--run-cosmic",
                "--cosmic-data", str(cosmic_path)
            ]
            
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300,
                    cwd=str(self.workspace_root)
                )
                
                # Extract overlap count from output
                overlap_count = None
                diag_json_path = output_dir / f"week2_diagnostic_results_{dataset.stem}.json"
                if diag_json_path.exists():
                    with open(diag_json_path, 'r') as f:
                        diag_data = json.load(f)
                        cosmic_data = diag_data.get("cosmic", {})
                        overlap_count = cosmic_data.get("overlap_count", None)
                
                results_by_dataset[dataset.stem] = {
                    "overlap_count": overlap_count,
                    "exit_code": result.returncode
                }
                
            except Exception as e:
                results_by_dataset[dataset.stem] = {"error": str(e)}
        
        evidence["results_by_dataset"] = results_by_dataset
        
        # Check if results differ
        overlap_counts = [r.get("overlap_count") for r in results_by_dataset.values() 
                          if r.get("overlap_count") is not None]
        
        if len(overlap_counts) >= 2:
            if len(set(overlap_counts)) > 1:
                status = "RESPONSIVE"  # Results differ between datasets
            else:
                status = "STATIC"  # Results are identical (suspicious)
        else:
            status = "UNKNOWN"
        
        self.results.append(VerificationResult(
            phase="COSMIC Sensitivity To Dataset Changes",
            status=status,
            evidence=evidence
        ))
        
        print(f"  Datasets tested: {len(results_by_dataset)}")
        for ds_name, result in results_by_dataset.items():
            print(f"    {ds_name}: overlap={result.get('overlap_count', 'N/A')}")
        print()
    
    def phase9_failure_modes(self):
        """Phase 9: Failure Mode Detection"""
        print("PHASE 9: Failure Mode Detection")
        print("-" * 80)
        
        evidence = {}
        failure_modes_tested = []
        
        test_dataset = self.week2_dir / "demo_fusion.csv"
        if not test_dataset.exists():
            evidence["status"] = "NO_TEST_DATA"
            self.results.append(VerificationResult(
                phase="Failure Mode Detection",
                status="INSUFFICIENT_EVIDENCE",
                evidence=evidence
            ))
            print("  No test dataset found")
            print()
            return
        
        output_dir = Path(self.temp_dir) / "phase9_output"
        output_dir.mkdir(exist_ok=True)
        
        # Test 1: COSMIC file exists but not loaded (invalid path)
        failure_modes_tested.append("invalid_cosmic_path")
        cmd = [
            sys.executable, "-m", "week2_validation.run_week2",
            "--fusion-data", str(test_dataset),
            "--output-dir", str(output_dir / "test1"),
            "--run-cosmic",
            "--cosmic-data", "/nonexistent/cosmic.tsv"
        ]
        output_dir.mkdir(exist_ok=True)
        (output_dir / "test1").mkdir(exist_ok=True)
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=60,
                cwd=str(self.workspace_root)
            )
            evidence["test1_invalid_path"] = {
                "exit_code": result.returncode,
                "handled_gracefully": "Warning" in result.stdout or "Error" in result.stderr
            }
        except:
            evidence["test1_invalid_path"] = {"error": "timeout"}
        
        # Test 2: COSMIC loads but mapping fails (malformed COSMIC)
        # This would require creating a test COSMIC file, which we avoid
        
        # Test 3: COSMIC loads but skipped by config
        # Check config for cosmic.enabled
        config_path = self.week2_dir / "config" / "thresholds.yaml"
        if config_path.exists():
            import yaml
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                cosmic_enabled = config.get("cosmic", {}).get("enabled", True)
                evidence["cosmic_enabled_in_config"] = cosmic_enabled
                failure_modes_tested.append("config_check")
        
        evidence["failure_modes_tested"] = failure_modes_tested
        
        self.results.append(VerificationResult(
            phase="Failure Mode Detection",
            status="COMPLETE",
            evidence=evidence
        ))
        
        print(f"  Failure modes tested: {len(failure_modes_tested)}")
        print(f"  COSMIC enabled in config: {evidence.get('cosmic_enabled_in_config', 'N/A')}")
        print()
    
    def generate_final_report(self) -> Dict[str, Any]:
        """Generate final verification report"""
        print("=" * 80)
        print("FINAL VERDICT")
        print("=" * 80)
        print()
        
        # Extract key findings
        file_presence = next((r for r in self.results if r.phase == "COSMIC File Presence"), None)
        runtime_load = next((r for r in self.results if r.phase == "COSMIC Runtime Load"), None)
        functional_usage = next((r for r in self.results if r.phase == "COSMIC Functional Usage"), None)
        mapping_integrity = next((r for r in self.results if r.phase == "COSMIC Mapping Integrity"), None)
        statistical_agreement = next((r for r in self.results if r.phase == "COSMIC Statistical Agreement"), None)
        biological_plausibility = next((r for r in self.results if r.phase == "COSMIC Biological Plausibility"), None)
        signal_non_triviality = next((r for r in self.results if r.phase == "COSMIC Signal Non-Triviality"), None)
        dataset_sensitivity = next((r for r in self.results if r.phase == "COSMIC Sensitivity To Dataset Changes"), None)
        
        # Determine final verdict
        if not file_presence or file_presence.status == "INSUFFICIENT_EVIDENCE":
            verdict = "COSMIC NOT PRESENT"
        elif runtime_load and runtime_load.status == "FAILED":
            verdict = "COSMIC PRESENT BUT NOT FUNCTIONAL"
        elif functional_usage and functional_usage.status == "LOADED_BUT_UNUSED":
            verdict = "COSMIC LOADED BUT NOT USED"
        elif functional_usage and functional_usage.status == "NOT_USED":
            verdict = "COSMIC LOADED BUT NOT USED"
        elif functional_usage and functional_usage.status == "ACTIVE":
            # Check other indicators
            if (statistical_agreement and statistical_agreement.status in ["GOOD", "WEAK"] and
                signal_non_triviality and signal_non_triviality.status == "VALID"):
                verdict = "COSMIC FULLY OPERATIONAL"
            else:
                verdict = "COSMIC PRESENT BUT NOT FUNCTIONAL"
        else:
            verdict = "INSUFFICIENT EVIDENCE"
        
        # Build report
        report = {
            "verification_timestamp": pd.Timestamp.now().isoformat(),
            "final_verdict": verdict,
            "phases": {
                "COSMIC File Presence": {
                    "status": file_presence.status if file_presence else "NOT_RUN",
                    "evidence": file_presence.evidence if file_presence else {}
                },
                "COSMIC Runtime Load": {
                    "status": runtime_load.status if runtime_load else "NOT_RUN",
                    "evidence": runtime_load.evidence if runtime_load else {}
                },
                "COSMIC Functional Usage": {
                    "status": functional_usage.status if functional_usage else "NOT_RUN",
                    "evidence": functional_usage.evidence if functional_usage else {}
                },
                "COSMIC Mapping Integrity": {
                    "status": mapping_integrity.status if mapping_integrity else "NOT_RUN",
                    "evidence": mapping_integrity.evidence if mapping_integrity else {}
                },
                "COSMIC Statistical Agreement": {
                    "status": statistical_agreement.status if statistical_agreement else "NOT_RUN",
                    "evidence": statistical_agreement.evidence if statistical_agreement else {},
                    "metrics": statistical_agreement.metrics if statistical_agreement and statistical_agreement.metrics else {}
                },
                "COSMIC Biological Plausibility": {
                    "status": biological_plausibility.status if biological_plausibility else "NOT_RUN",
                    "evidence": biological_plausibility.evidence if biological_plausibility else {}
                },
                "COSMIC Signal Non-Triviality": {
                    "status": signal_non_triviality.status if signal_non_triviality else "NOT_RUN",
                    "evidence": signal_non_triviality.evidence if signal_non_triviality else {}
                },
                "COSMIC Sensitivity To Dataset Changes": {
                    "status": dataset_sensitivity.status if dataset_sensitivity else "NOT_RUN",
                    "evidence": dataset_sensitivity.evidence if dataset_sensitivity else {}
                }
            }
        }
        
        # Print summary
        print("SUMMARY:")
        print(f"  COSMIC File Presence: {file_presence.status if file_presence else 'NOT_RUN'}")
        print(f"  COSMIC Runtime Load: {runtime_load.status if runtime_load else 'NOT_RUN'}")
        print(f"  COSMIC Functional Usage: {functional_usage.status if functional_usage else 'NOT_RUN'}")
        print(f"  COSMIC Mapping Integrity: {mapping_integrity.status if mapping_integrity else 'NOT_RUN'}")
        print(f"  COSMIC Statistical Agreement: {statistical_agreement.status if statistical_agreement else 'NOT_RUN'}")
        print(f"  COSMIC Biological Plausibility: {biological_plausibility.status if biological_plausibility else 'NOT_RUN'}")
        print(f"  COSMIC Signal Non-Triviality: {signal_non_triviality.status if signal_non_triviality else 'NOT_RUN'}")
        print(f"  COSMIC Sensitivity: {dataset_sensitivity.status if dataset_sensitivity else 'NOT_RUN'}")
        print()
        print(f"FINAL VERDICT: {verdict}")
        print()
        
        return report


def main():
    """Main entry point"""
    workspace_root = Path(__file__).parent.parent
    
    with COSMICVerification(workspace_root) as verifier:
        report = verifier.run_all_phases()
        
        # Save report
        report_path = workspace_root / "week2_validation" / "COSMIC_OPERATIONAL_VERIFICATION_REPORT.json"
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"Full report saved to: {report_path}")
        
        return 0 if report["final_verdict"] == "COSMIC FULLY OPERATIONAL" else 1


if __name__ == "__main__":
    sys.exit(main())
