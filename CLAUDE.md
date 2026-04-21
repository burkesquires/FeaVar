# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FeaVar (Feature Variant) is a Python bioinformatics package for computing **Sequence Feature Variant Types (SFVTs)** in multiple sequence alignments (MSAs). Given an alignment and a set of positions, it extracts the residues at those positions for every sequence, groups sequences by their unique substring patterns, and assigns variant type names.

## Commands

```bash
# Install for development (includes all optional deps)
pip install ".[all]"

# Run all tests
pytest

# Run a single test
pytest tests/test_parser.py::TestPositionParser::test_parse_range -v

# Skip slow/integration tests
pytest -m "not slow and not integration"

# With coverage
pytest --cov=feavar --cov-report=html

# Format
black feavar tests

# Lint
ruff check feavar tests

# Type check
mypy feavar

# CLI
feavar --help
feavar analyze alignment.clw --reference REF001 --positions "100-110" --naming delta
```

## Architecture

The package lives entirely in `feavar/` (ignore the legacy `FeaVar/` directory at the root which is superseded). Entry point for CLI is `feavar.cli:run_cli` (registered as the `feavar` command).

### Data flow

```
alignment file
    → AlignmentHandler  (loads MSA via BioPython, locates reference, maps gaps)
    → PositionParser    (parses "100-110", "5,10,15", or mixed into sorted int list)
    → gap-adjusted positions
    → AlignmentHandler.extract_all_variants()  (generator over sequences)
    → FeaVarAnalysis    (orchestrates above, builds DataFrames)
    → AnalysisResult    (variants_df + summary_df, optional metadata merge)
```

### Key modules

| Module | Responsibility |
|---|---|
| `core.py` | `FeaVarAnalysis` orchestrator + `AnalysisResult` dataclass |
| `alignment.py` | `AlignmentHandler` (MSA I/O, gap-aware position mapping) + `ReferenceSequence` dataclass |
| `parser.py` | `PositionParser` — parses position strings, validates bounds |
| `naming.py` | `NamingScheme` enum + Delta/Ranked naming logic |
| `visualization.py` | `VariantPlotter` — matplotlib bar/stacked charts (optional dep) |
| `exceptions.py` | Exception hierarchy rooted at `FeaVarError` |

### Naming schemes

Two mutually exclusive schemes (pass `naming_scheme=` to `FeaVarAnalysis`):

- **Delta** (default, stable): name encodes position + residue change, e.g. `VT-5G.12C`. Survives resequencing campaigns without renaming.
- **Ranked**: traditional `VT-001`, `VT-002` ordered by frequency — simpler but changes when new sequences are added.

### Gap-aware positions

Positions are always **1-based relative to the reference sequence** (excluding gaps). `AlignmentHandler` translates them to alignment column indices by walking the reference and skipping gap characters (`-`). This is the core domain concept and the reason for the `ReferenceSequence` dataclass.

### `AnalysisResult`

The return value of `FeaVarAnalysis.run()`. Its two DataFrames:

- `variants_df` — one row per sequence: `accession`, `variant_type` (raw substring)
- `summary_df` — one row per unique variant type: `variant_type`, `count`, `vt_name`

Useful methods: `get_top_variants(n)`, `get_variant_for_accession(id)`.

## Conventions

- **Type hints everywhere**; mypy runs in strict mode (`disallow_untyped_defs = true`).
- **Black** (line length 88) + **Ruff** handle all formatting/linting. `ruff` rules: E, W, F, I, B, C4, UP.
- Test markers: `@pytest.mark.slow` and `@pytest.mark.integration` — must be declared in `pyproject.toml` (`--strict-markers` is set).
- Alignment format is inferred from file extension by `AlignmentHandler`; supported: Clustal (`.clw`, `.clustal`), FASTA, Phylip, Nexus, Stockholm.
- `VariantPlotter` methods return `self` for chaining; matplotlib import is deferred and raises a clear error if the `[plot]` extra is missing.
- Custom exceptions from `exceptions.py` should be raised (not bare `ValueError`/`RuntimeError`) so callers can catch domain-specific errors.
