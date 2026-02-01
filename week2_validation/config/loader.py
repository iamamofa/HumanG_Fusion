"""
config/loader.py - Safe configuration loading and validation for Week 2: Data Integrity & Statistical Validation pipeline.

Loads and validates thresholds.yaml with strict validation.
Does not infer defaults. Treats configuration as read-only.

WHAT DOES THIS FILE DO?
This module reads the thresholds.yaml configuration file and converts it
into Python objects that the pipeline can use. It validates every field
to ensure the configuration is correct before the pipeline runs.

The configuration controls:
- Whether real data analysis is allowed (freeze state)
- Benford's Law test settings
- Distribution test settings
- Visualization rules
- COSMIC cross-validation rules
- Reporting rules
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict

import yaml


# =============================================================================
# CUSTOM ERROR TYPE
# =============================================================================

class ConfigError(Exception):
    """
    Raised when configuration loading or validation fails.
    
    This error indicates that something is wrong with the thresholds.yaml
    file - either it's missing, malformed, or contains invalid values.
    """
    pass


# =============================================================================
# NESTED DATACLASSES FOR CONFIGURATION SECTIONS
# =============================================================================

@dataclass(frozen=True)
class MinimumSampleSize:
    """
    Minimum sample size requirements for different tests.
    
    These thresholds ensure statistical tests are only run when there
    is enough data for meaningful results.
    """
    benford: int                # Minimum N for Benford's Law applicability
    distribution_tests: int     # Minimum N for KS / Anderson-Darling tests


@dataclass(frozen=True)
class BenfordDecisionLabels:
    """
    Labels used to describe Benford's Law test outcomes.
    
    These provide human-readable descriptions of test results.
    """
    pass_label: str         # Label when data is Benford-consistent
    fail_label: str         # Label when data is Benford-inconsistent
    inconclusive: str       # Label when test is not applicable


@dataclass(frozen=True)
class BenfordConfig:
    """
    Configuration for Benford's Law analysis.
    
    Controls how Benford's Law tests are performed and interpreted.
    """
    enabled: bool                           # Whether Benford testing is enabled
    first_significant_digit_only: bool      # Only test first digit (not second)
    goodness_of_fit_test: str               # Test type (e.g., "chi_squared")
    significance_level: float               # Alpha level for hypothesis test
    interpretation_note: str                # Important context for results
    decision_labels: BenfordDecisionLabels  # Labels for test outcomes


@dataclass(frozen=True)
class SimulationConfig:
    """
    Configuration for simulation and reproducibility.
    
    Controls random seed behavior to ensure reproducible results.
    """
    random_seed_policy: str     # Policy for random seeds ("fixed" or "random")
    default_seed: int           # Default seed value when policy is "fixed"
    note: str                   # Explanation of seed policy


@dataclass(frozen=True)
class StatisticalTestConfig:
    """
    Configuration for a single statistical test (KS or Anderson-Darling).
    """
    enabled: bool               # Whether this test is enabled
    significance_level: float   # Alpha level for hypothesis test


@dataclass(frozen=True)
class DistributionTestsConfig:
    """
    Configuration for distribution testing (KS and Anderson-Darling tests).
    
    These tests assess whether data follows expected distributions.
    """
    ks_test: StatisticalTestConfig          # Kolmogorov-Smirnov test settings
    anderson_darling: StatisticalTestConfig # Anderson-Darling test settings
    interpretation_note: str                # Important context for results


@dataclass(frozen=True)
class VisualizationConfig:
    """
    Configuration for visualization/plotting rules.
    
    Controls what types of plots can be generated and saved.
    """
    allow_log_scale: bool           # Allow logarithmic scale plots
    allow_linear_scale: bool        # Allow linear scale plots
    save_real_data_plots: bool      # Allow saving plots of real data


@dataclass(frozen=True)
class CosmicConfig:
    """
    Configuration for COSMIC database cross-validation.
    
    Controls how fusion genes are validated against the COSMIC database.
    """
    enabled: bool                       # Whether COSMIC validation is enabled
    require_exact_gene_match: bool      # Require exact gene name matches
    allow_partial_matches: bool         # Allow partial/fuzzy matches
    save_results_before_freeze: bool    # Allow saving results before freeze


@dataclass(frozen=True)
class ReportingConfig:
    """
    Configuration for reporting rules.
    
    Controls what interpretations and approvals are allowed before freeze.
    """
    allow_interpretation_before_freeze: bool    # Allow interpreting results
    allow_dataset_approval_before_freeze: bool  # Allow approving datasets


@dataclass(frozen=True)
class Week2Config:
    """
    Complete, immutable configuration for the Data Integrity & Statistical Validation pipeline.
    
    This is the main configuration object that contains all settings
    needed to run the pipeline. All fields are validated upon loading.
    """
    # Top-level freeze state controls
    dataset_frozen_required: bool       # Whether dataset must be frozen
    allow_real_data_analysis: bool      # Whether real data analysis is allowed
    
    # Section configurations
    minimum_sample_size: MinimumSampleSize
    benford: BenfordConfig
    simulation: SimulationConfig
    distribution_tests: DistributionTestsConfig
    visualization: VisualizationConfig
    cosmic: CosmicConfig
    reporting: ReportingConfig


# =============================================================================
# REQUIRED KEYS FOR VALIDATION
# =============================================================================

REQUIRED_TOP_LEVEL_KEYS = frozenset({
    "dataset_frozen_required",
    "allow_real_data_analysis",
    "minimum_sample_size",
    "benford",
    "simulation",
    "distribution_tests",
    "visualization",
    "cosmic",
    "reporting",
})


# =============================================================================
# VALIDATION HELPER FUNCTIONS
# =============================================================================

def _validate_yaml_structure(data: Any) -> Dict[str, Any]:
    """
    Validate that loaded YAML is a dictionary.
    
    Args:
        data: Parsed YAML data.
    
    Returns:
        Validated dictionary.
    
    Raises:
        ConfigError: If data is not a dictionary.
    """
    if data is None:
        raise ConfigError("Configuration file is empty")
    
    if not isinstance(data, dict):
        raise ConfigError(
            f"Configuration must be a mapping, got {type(data).__name__}"
        )
    
    return data


def _validate_top_level_keys(config: Dict[str, Any]) -> None:
    """
    Validate that all required top-level keys are present.
    
    Args:
        config: Configuration dictionary.
    
    Raises:
        ConfigError: If required keys are missing.
    """
    present_keys = set(config.keys())
    missing_keys = REQUIRED_TOP_LEVEL_KEYS - present_keys
    
    if missing_keys:
        raise ConfigError(
            f"Missing required configuration keys: {sorted(missing_keys)}"
        )


def _validate_bool(value: Any, field_name: str) -> bool:
    """
    Validate and return a boolean value.
    
    Args:
        value: Value to validate.
        field_name: Name of the field for error messages.
    
    Returns:
        Validated boolean.
    
    Raises:
        ConfigError: If value is not a boolean.
    """
    if not isinstance(value, bool):
        raise ConfigError(
            f"'{field_name}' must be a boolean, got {type(value).__name__}"
        )
    return value


def _validate_int(value: Any, field_name: str, min_value: int = None) -> int:
    """
    Validate and return an integer value.
    
    Args:
        value: Value to validate.
        field_name: Name of the field for error messages.
        min_value: Optional minimum allowed value.
    
    Returns:
        Validated integer.
    
    Raises:
        ConfigError: If value is not a valid integer.
    """
    # Check for boolean first (bool is subclass of int in Python)
    if isinstance(value, bool):
        raise ConfigError(
            f"'{field_name}' must be an integer, got boolean"
        )
    
    if not isinstance(value, int):
        raise ConfigError(
            f"'{field_name}' must be an integer, got {type(value).__name__}"
        )
    
    if min_value is not None and value < min_value:
        raise ConfigError(
            f"'{field_name}' must be at least {min_value}, got {value}"
        )
    
    return value


def _validate_float(value: Any, field_name: str, 
                    min_value: float = None, max_value: float = None) -> float:
    """
    Validate and return a float value.
    
    Args:
        value: Value to validate.
        field_name: Name of the field for error messages.
        min_value: Optional minimum allowed value.
        max_value: Optional maximum allowed value.
    
    Returns:
        Validated float.
    
    Raises:
        ConfigError: If value is not a valid number.
    """
    if isinstance(value, bool):
        raise ConfigError(
            f"'{field_name}' must be a number, got boolean"
        )
    
    if not isinstance(value, (int, float)):
        raise ConfigError(
            f"'{field_name}' must be a number, got {type(value).__name__}"
        )
    
    float_value = float(value)
    
    if min_value is not None and float_value < min_value:
        raise ConfigError(
            f"'{field_name}' must be at least {min_value}, got {float_value}"
        )
    
    if max_value is not None and float_value > max_value:
        raise ConfigError(
            f"'{field_name}' must be at most {max_value}, got {float_value}"
        )
    
    return float_value


def _validate_string(value: Any, field_name: str) -> str:
    """
    Validate and return a string value.
    
    Args:
        value: Value to validate.
        field_name: Name of the field for error messages.
    
    Returns:
        Validated string.
    
    Raises:
        ConfigError: If value is not a string.
    """
    if not isinstance(value, str):
        raise ConfigError(
            f"'{field_name}' must be a string, got {type(value).__name__}"
        )
    return value


def _validate_dict(value: Any, field_name: str) -> Dict[str, Any]:
    """
    Validate and return a dictionary value.
    
    Args:
        value: Value to validate.
        field_name: Name of the field for error messages.
    
    Returns:
        Validated dictionary.
    
    Raises:
        ConfigError: If value is not a dictionary.
    """
    if not isinstance(value, dict):
        raise ConfigError(
            f"'{field_name}' must be a mapping, got {type(value).__name__}"
        )
    return value


# =============================================================================
# SECTION PARSING FUNCTIONS
# =============================================================================

def _parse_minimum_sample_size(data: Dict[str, Any]) -> MinimumSampleSize:
    """Parse and validate the minimum_sample_size section."""
    section = _validate_dict(data.get("minimum_sample_size"), "minimum_sample_size")
    
    return MinimumSampleSize(
        benford=_validate_int(section.get("benford"), "minimum_sample_size.benford", min_value=1),
        distribution_tests=_validate_int(section.get("distribution_tests"), "minimum_sample_size.distribution_tests", min_value=1),
    )


def _parse_benford_decision_labels(data: Dict[str, Any]) -> BenfordDecisionLabels:
    """Parse and validate the benford.decision_labels section."""
    section = _validate_dict(data, "benford.decision_labels")
    
    return BenfordDecisionLabels(
        pass_label=_validate_string(section.get("pass"), "benford.decision_labels.pass"),
        fail_label=_validate_string(section.get("fail"), "benford.decision_labels.fail"),
        inconclusive=_validate_string(section.get("inconclusive"), "benford.decision_labels.inconclusive"),
    )


def _parse_benford(data: Dict[str, Any]) -> BenfordConfig:
    """Parse and validate the benford section."""
    section = _validate_dict(data.get("benford"), "benford")
    
    return BenfordConfig(
        enabled=_validate_bool(section.get("enabled"), "benford.enabled"),
        first_significant_digit_only=_validate_bool(section.get("first_significant_digit_only"), "benford.first_significant_digit_only"),
        goodness_of_fit_test=_validate_string(section.get("goodness_of_fit_test"), "benford.goodness_of_fit_test"),
        significance_level=_validate_float(section.get("significance_level"), "benford.significance_level", min_value=0.0, max_value=1.0),
        interpretation_note=_validate_string(section.get("interpretation_note"), "benford.interpretation_note"),
        decision_labels=_parse_benford_decision_labels(section.get("decision_labels")),
    )


def _parse_simulation(data: Dict[str, Any]) -> SimulationConfig:
    """Parse and validate the simulation section."""
    section = _validate_dict(data.get("simulation"), "simulation")
    
    return SimulationConfig(
        random_seed_policy=_validate_string(section.get("random_seed_policy"), "simulation.random_seed_policy"),
        default_seed=_validate_int(section.get("default_seed"), "simulation.default_seed", min_value=0),
        note=_validate_string(section.get("note"), "simulation.note"),
    )


def _parse_statistical_test(data: Dict[str, Any], prefix: str) -> StatisticalTestConfig:
    """Parse and validate a statistical test configuration."""
    section = _validate_dict(data, prefix)
    
    return StatisticalTestConfig(
        enabled=_validate_bool(section.get("enabled"), f"{prefix}.enabled"),
        significance_level=_validate_float(section.get("significance_level"), f"{prefix}.significance_level", min_value=0.0, max_value=1.0),
    )


def _parse_distribution_tests(data: Dict[str, Any]) -> DistributionTestsConfig:
    """Parse and validate the distribution_tests section."""
    section = _validate_dict(data.get("distribution_tests"), "distribution_tests")
    
    return DistributionTestsConfig(
        ks_test=_parse_statistical_test(section.get("ks_test"), "distribution_tests.ks_test"),
        anderson_darling=_parse_statistical_test(section.get("anderson_darling"), "distribution_tests.anderson_darling"),
        interpretation_note=_validate_string(section.get("interpretation_note"), "distribution_tests.interpretation_note"),
    )


def _parse_visualization(data: Dict[str, Any]) -> VisualizationConfig:
    """Parse and validate the visualization section."""
    section = _validate_dict(data.get("visualization"), "visualization")
    
    return VisualizationConfig(
        allow_log_scale=_validate_bool(section.get("allow_log_scale"), "visualization.allow_log_scale"),
        allow_linear_scale=_validate_bool(section.get("allow_linear_scale"), "visualization.allow_linear_scale"),
        save_real_data_plots=_validate_bool(section.get("save_real_data_plots"), "visualization.save_real_data_plots"),
    )


def _parse_cosmic(data: Dict[str, Any]) -> CosmicConfig:
    """Parse and validate the cosmic section."""
    section = _validate_dict(data.get("cosmic"), "cosmic")
    
    return CosmicConfig(
        enabled=_validate_bool(section.get("enabled"), "cosmic.enabled"),
        require_exact_gene_match=_validate_bool(section.get("require_exact_gene_match"), "cosmic.require_exact_gene_match"),
        allow_partial_matches=_validate_bool(section.get("allow_partial_matches"), "cosmic.allow_partial_matches"),
        save_results_before_freeze=_validate_bool(section.get("save_results_before_freeze"), "cosmic.save_results_before_freeze"),
    )


def _parse_reporting(data: Dict[str, Any]) -> ReportingConfig:
    """Parse and validate the reporting section."""
    section = _validate_dict(data.get("reporting"), "reporting")
    
    return ReportingConfig(
        allow_interpretation_before_freeze=_validate_bool(section.get("allow_interpretation_before_freeze"), "reporting.allow_interpretation_before_freeze"),
        allow_dataset_approval_before_freeze=_validate_bool(section.get("allow_dataset_approval_before_freeze"), "reporting.allow_dataset_approval_before_freeze"),
    )


# =============================================================================
# MAIN LOADING FUNCTION
# =============================================================================

def load_config(config_path: Path) -> Week2Config:
    """
    Load and validate the Data Integrity & Statistical Validation configuration from a YAML file.
    
    This is the main entry point for loading configuration. It reads the
    thresholds.yaml file, validates all fields, and returns an immutable
    configuration object.
    
    Args:
        config_path: Path to the thresholds.yaml file.
    
    Returns:
        Validated, immutable Week2Config object containing all settings.
    
    Raises:
        ConfigError: If the file cannot be read or validation fails.
    
    Example:
        config = load_config(Path("config/thresholds.yaml"))
        if config.benford.enabled:
            # Run Benford analysis...
    """
    # Validate the path argument
    if config_path is None:
        raise ConfigError("Config path cannot be None")
    
    if not isinstance(config_path, Path):
        raise ConfigError("Config path must be a Path object")
    
    if not config_path.exists():
        raise ConfigError(f"Configuration file does not exist: {config_path}")
    
    if not config_path.is_file():
        raise ConfigError(f"Configuration path is not a regular file: {config_path}")
    
    # Read and parse the YAML file
    try:
        with config_path.open("r", encoding="utf-8") as f:
            raw_data = yaml.safe_load(f)
    except PermissionError as exc:
        raise ConfigError("Permission denied reading configuration file") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML syntax: {exc}") from exc
    except OSError as exc:
        raise ConfigError(f"Error reading configuration file: {exc}") from exc
    
    # Validate structure and keys
    config = _validate_yaml_structure(raw_data)
    _validate_top_level_keys(config)
    
    # Parse and validate all sections, building the final config object
    return Week2Config(
        dataset_frozen_required=_validate_bool(config.get("dataset_frozen_required"), "dataset_frozen_required"),
        allow_real_data_analysis=_validate_bool(config.get("allow_real_data_analysis"), "allow_real_data_analysis"),
        minimum_sample_size=_parse_minimum_sample_size(config),
        benford=_parse_benford(config),
        simulation=_parse_simulation(config),
        distribution_tests=_parse_distribution_tests(config),
        visualization=_parse_visualization(config),
        cosmic=_parse_cosmic(config),
        reporting=_parse_reporting(config),
    )


# =============================================================================
# CONVENIENCE FUNCTION FOR DEFAULT PATH
# =============================================================================

def load_default_config() -> Week2Config:
    """
    Load configuration from the default location (config/thresholds.yaml).
    
    This is a convenience function that loads from the standard config
    location relative to this file.
    
    Returns:
        Validated, immutable Week2Config object.
    
    Raises:
        ConfigError: If the file cannot be read or validation fails.
    """
    # Get the directory where this loader.py file is located
    config_dir = Path(__file__).parent
    default_path = config_dir / "thresholds.yaml"
    
    return load_config(default_path)
