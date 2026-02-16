# Week 2 Integration Contract

## Input Contract (what Week 2 expects from Week 1)
- File: CSV, TSV, JSON, Parquet, or XLSX
- Required columns: fusion_id, gene_1, gene_2, protein_length, recurrence_count
- Constraints: protein_length > 0, recurrence_count >= 0, no nulls in required columns
- Alternative: geneA/geneB/samples_detected format (auto-adapted)

## Output Contract (what Week 2 provides to Week 3)
- validation_status_{stem}.json  -- exit code, status, quality gates
- validation_cleaned_dataset_{stem}.csv  -- frozen validated data (unchanged)
- validation_dataset_certification_{stem}.json  -- certification for power-law modeling
- validation_diagnostic_results_{stem}.json  -- all diagnostic metrics
- run_manifest_{stem}.json  -- complete run metadata

## Exit Codes
- 0: Success (APPROVED or CONDITIONAL)
- 10: Input schema error (bad columns or types)
- 20: Freeze error
- 25: Quality gate rejected
- 30: Diagnostic runtime error
- 40: Config error
- 99: Unknown error

## Week 3 should check:
- Exit code == 0 before proceeding
- validation_status_{stem}.json -> status == "SUCCESS"
- validation_dataset_certification_{stem}.json -> approved_for_powerlaw_modeling == true
