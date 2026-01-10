# FeaVar

[![PyPI version](https://badge.fury.io/py/feavar.svg)](https://badge.fury.io/py/feavar)
[![Python versions](https://img.shields.io/pypi/pyversions/feavar.svg)](https://pypi.org/project/feavar/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/burkesquires/feavar/workflows/Tests/badge.svg)](https://github.com/burkesquires/feavar/actions)

**FeaVar** (Feature Variant) is a Python package for computing clusters of unique substrings or "Sequence Feature Variant Types" (SFVTs) for user-selected positions in a set of aligned sequences.

## Features

- **SFVT Computation**: Identifies unique variant types based on user-defined positions in aligned sequences
- **Reference-Aware**: Automatically adjusts positions to account for gaps in the reference sequence
- **Metadata Integration**: Merge variant results with external metadata for downstream analysis
- **Visualization**: Generate bar charts and stacked plots for variant distributions
- **Flexible Input**: Supports Clustal, FASTA, and other alignment formats via Biopython
- **Modern Python**: Type hints, dataclasses, and clean API design

## Installation

### From PyPI (recommended)

```bash
pip install feavar
```

### With plotting support

```bash
pip install feavar[plot]
```

### From source

```bash
git clone https://github.com/burkesquires/feavar.git
cd feavar
pip install -e ".[dev]"
```

## Quick Start

### Command Line

```bash
# Basic usage (uses delta naming by default)
feavar -a alignment.clw -r CY021716 -p "124-142"

# Use traditional ranked naming (VT-001, VT-002, etc.)
feavar -a alignment.clw -r CY021716 -p "124-142" -n ranked

# With metadata
feavar -a alignment.clw -r CY021716 -p "124-142" -m metadata.tsv -t 15

# Save to specific directory
feavar -a sequences.fasta -r REF001 -p "50-100" -o results/
```

### Python API

```python
from feavar import FeaVarAnalysis

# Create and run analysis (delta naming by default)
analysis = FeaVarAnalysis(
    alignment_path="sequences.clw",
    reference_id="CY021716",
    positions="124-142",
)

result = analysis.run()

# View results
print(f"Found {result.num_variant_types} unique variants")
print(result.summary_df.head(10))

# Get top variants
top_variants = result.get_top_variants(n=10)
print(top_variants)

# Save results
analysis.save_results(prefix="my_analysis")

# Merge with metadata
merged_df = analysis.merge_metadata("metadata.tsv")
```

### Using Traditional Ranked Naming

```python
# Use ranked naming scheme (VT-001, VT-002, etc.)
analysis = FeaVarAnalysis(
    alignment_path="sequences.clw",
    reference_id="CY021716",
    positions="124-142",
    naming_scheme="ranked",  # Traditional frequency-based naming
)
```
```

### Visualization

```python
from feavar import FeaVarAnalysis, VariantPlotter

# Run analysis
analysis = FeaVarAnalysis(
    alignment_path="sequences.clw",
    reference_id="CY021716", 
    positions="124-142",
)
result = analysis.run()

# Create plots
plotter = VariantPlotter(output_dir="plots")

# Bar chart of variant counts
plotter.plot_variant_counts(result.summary_df, top_n=10)
plotter.save("variant_counts.png")

# Horizontal distribution
plotter.plot_variant_distribution(result.summary_df, top_n=20)
plotter.save("variant_distribution.png")
```

## Position Syntax

Positions can be specified in several formats:

| Format | Example | Description |
|--------|---------|-------------|
| Single | `100` | Single position |
| List | `100,110,120` | Multiple positions |
| Range | `100-110` | Inclusive range |
| Mixed | `100-110,120,130-140` | Combination |

All positions are **1-based** (first position is 1, not 0).

### Supported Alignment Formats

FeaVar automatically detects the alignment format from the file extension:

| Extension | Format | Description |
|-----------|--------|-------------|
| `.clw`, `.clustal`, `.aln` | clustal | Clustal W/X format |
| `.fasta`, `.fa`, `.fas`, `.fna`, `.faa` | fasta | FASTA format |
| `.phy`, `.phylip`, `.ph` | phylip | PHYLIP format |
| `.nex`, `.nexus`, `.nxs` | nexus | NEXUS format |
| `.sto`, `.stockholm`, `.stk` | stockholm | Stockholm format |

You can also explicitly specify the format using the `-f` option.

## Variant Naming Schemes

FeaVar supports two naming schemes for variant types:

### Delta Naming (Default)

The **delta** scheme generates stable, reproducible names based on differences from the reference sequence. Names remain consistent regardless of variant frequency or the order sequences appear.

**Format:** `VT-{position}{residue}.{position}{residue}...`

| Variant Sequence | Reference | Delta Name |
|-----------------|-----------|------------|
| `GAAGACAGG` | `GAAGACAGG` | `VT-REF` |
| `GAAGACAGT` | `GAAGACAGG` | `VT-9T` |
| `GAAGGCAGT` | `GAAGACAGG` | `VT-5G.9T` |
| `TAAGACAGG` | `GAAGACAGG` | `VT-1T` |

**Advantages:**
- Names are stable across different runs
- Names are meaningful (show what changed)
- Adding new sequences doesn't change existing names
- Easy to compare variants

### Ranked Naming (Traditional)

The **ranked** scheme assigns names based on frequency (most common = VT-001).

| Variant | Count | Ranked Name |
|---------|-------|-------------|
| Most common | 150 | `VT-001` |
| Second most | 89 | `VT-002` |
| Third most | 45 | `VT-003` |

**Note:** Ranked names may change when new sequences are added, as frequencies can shift.

### Usage

```bash
# Use delta naming (default)
feavar -a sequences.clw -r REF001 -p "100-120"

# Use ranked naming
feavar -a sequences.clw -r REF001 -p "100-120" -n ranked
```

```python
# Delta naming (default)
analysis = FeaVarAnalysis(
    alignment_path="sequences.clw",
    reference_id="REF001",
    positions="100-120",
)

# Ranked naming
analysis = FeaVarAnalysis(
    alignment_path="sequences.clw",
    reference_id="REF001",
    positions="100-120",
    naming_scheme="ranked",
)
```

## API Reference

### FeaVarAnalysis

The main analysis class.

```python
FeaVarAnalysis(
    alignment_path: str,      # Path to alignment file
    reference_id: str,        # Reference sequence identifier
    positions: str,           # Position string (e.g., "100-110")
    alignment_format: str = None,  # Alignment format (None = auto-detect)
    output_dir: str = ".",    # Output directory
    naming_scheme: str = "delta",  # Naming scheme ("delta" or "ranked")
)
```

**Methods:**

- `run() -> AnalysisResult`: Run the analysis
- `validate() -> bool`: Validate inputs before running
- `merge_metadata(path) -> DataFrame`: Merge results with metadata
- `save_results(prefix) -> Dict[str, Path]`: Save results to CSV files

### AnalysisResult

Container for analysis results.

**Attributes:**

- `variants_df`: DataFrame with (accession, variant_type) columns
- `summary_df`: DataFrame with (variant_type, count, VT) columns
- `num_sequences`: Total sequences analyzed
- `num_variant_types`: Unique variant types found

**Methods:**

- `get_top_variants(n) -> DataFrame`: Get top N variants by count
- `get_variant_for_accession(id) -> str`: Get variant for specific accession

### PositionParser

Utility for parsing position strings.

```python
from feavar import PositionParser

parser = PositionParser()
positions = parser.parse("100-110,120,130")
# Returns: [100, 101, 102, ..., 110, 120, 130]
```

### AlignmentHandler

Low-level alignment handling.

```python
from feavar import AlignmentHandler

handler = AlignmentHandler("sequences.clw", "clustal")
ref = handler.get_reference("SEQ001")
print(ref.length)  # Sequence length without gaps
```

## Development

### Setup

```bash
git clone https://github.com/burkesquires/feavar.git
cd feavar
pip install -e ".[all]"
```

### Running Tests

```bash
# Run all tests
pytest

# With coverage
pytest --cov=feavar --cov-report=html

# Run specific test file
pytest tests/test_parser.py -v
```

### Code Quality

```bash
# Format code
black feavar tests

# Lint
ruff check feavar tests

# Type checking
mypy feavar
```

## References

- Noronha, J. M., et al. (2012). Influenza virus sequence feature variant type analysis. *Journal of Virology*, 86(10), 5857–5866. [doi:10.1128/JVI.06901-11](https://doi.org/10.1128/JVI.06901-11)

- Karp, D. R., et al. (2009). Novel sequence feature variant type analysis of the HLA genetic association in systemic sclerosis. *Human Molecular Genetics*, 19(4), 707–719. [doi:10.1093/hmg/ddp521](https://doi.org/10.1093/hmg/ddp521)

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Changelog

### v2.0.0 (2024)

- Complete refactoring with modern Python practices
- New modular architecture with separate components
- Type hints throughout
- Comprehensive test suite
- Click-based CLI
- Improved error handling with custom exceptions
- Better documentation
