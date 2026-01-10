"""
Core FeaVar analysis module.

Provides the main FeaVarAnalysis class for computing sequence feature variant types.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any

import pandas as pd

from feavar.alignment import AlignmentHandler, ReferenceSequence
from feavar.parser import PositionParser
from feavar.exceptions import FeaVarError, MetadataError
from feavar.naming import (
    NamingScheme,
    compute_variant_names,
    get_reference_variant_sequence,
    generate_ranked_name,
)

logger = logging.getLogger(__name__)


@dataclass
class AnalysisResult:
    """
    Container for FeaVar analysis results.
    
    Attributes
    ----------
    variants_df : pd.DataFrame
        DataFrame with columns ['accession', 'variant_type'].
    summary_df : pd.DataFrame
        DataFrame with columns ['variant_type', 'count', 'VT'].
    reference_id : str
        The reference sequence identifier used.
    positions : List[int]
        Original (unadjusted) positions analyzed.
    adjusted_positions : List[int]
        Positions adjusted for gaps in reference.
    num_sequences : int
        Total number of sequences analyzed.
    num_variant_types : int
        Number of unique variant types found.
    """
    variants_df: pd.DataFrame
    summary_df: pd.DataFrame
    reference_id: str
    positions: List[int]
    adjusted_positions: List[int]
    num_sequences: int = field(init=False)
    num_variant_types: int = field(init=False)
    
    def __post_init__(self):
        self.num_sequences = len(self.variants_df)
        self.num_variant_types = len(self.summary_df)
    
    def get_top_variants(self, n: int = 10) -> pd.DataFrame:
        """Get the top N most common variants."""
        return self.summary_df.head(n)
    
    def get_variant_for_accession(self, accession: str) -> Optional[str]:
        """Get the variant type for a specific accession."""
        matches = self.variants_df[
            self.variants_df['accession'].str.contains(accession)
        ]
        if len(matches) > 0:
            return matches.iloc[0]['variant_type']
        return None


class FeaVarAnalysis:
    """
    Main class for Sequence Feature Variant Type analysis.
    
    This class orchestrates the complete FeaVar workflow:
    1. Load and validate alignment
    2. Parse and validate positions
    3. Adjust positions for reference gaps
    4. Extract variant types
    5. Compute statistics and optionally merge with metadata
    
    Parameters
    ----------
    alignment_path : str
        Path to the multiple sequence alignment file.
    reference_id : str
        Identifier of the reference sequence.
    positions : str
        Position string (e.g., "100-110" or "100,110,120").
    alignment_format : str, optional
        Alignment file format. If None, auto-detects from extension.
    output_dir : str or Path, optional
        Directory for output files (default: current directory).
    naming_scheme : str, optional
        Variant naming scheme: "delta" (default) or "ranked".
        - "delta": Stable names based on differences from reference (e.g., VT-5G.12C)
        - "ranked": Traditional rank-based names (e.g., VT-001 for most common)
        
    Examples
    --------
    >>> analysis = FeaVarAnalysis(
    ...     alignment_path="sequences.clw",
    ...     reference_id="CY021716",
    ...     positions="124-142"
    ... )
    >>> result = analysis.run()
    >>> print(result.summary_df.head())
    
    # Use ranked naming (traditional)
    >>> analysis = FeaVarAnalysis(
    ...     alignment_path="sequences.clw",
    ...     reference_id="CY021716",
    ...     positions="124-142",
    ...     naming_scheme="ranked"
    ... )
    """
    
    def __init__(
        self,
        alignment_path: str,
        reference_id: str,
        positions: str,
        alignment_format: Optional[str] = None,
        output_dir: Optional[str] = None,
        naming_scheme: str = "delta",
    ):
        """Initialize the FeaVar analysis."""
        self.alignment_path = alignment_path
        self.reference_id = reference_id
        self.positions_str = positions
        self.alignment_format = alignment_format  # None means auto-detect
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        
        # Set naming scheme
        try:
            self.naming_scheme = NamingScheme(naming_scheme.lower())
        except ValueError:
            valid = [s.value for s in NamingScheme]
            raise ValueError(
                f"Invalid naming scheme '{naming_scheme}'. Must be one of: {valid}"
            )
        
        # Components
        self._parser = PositionParser()
        self._handler: Optional[AlignmentHandler] = None
        self._reference: Optional[ReferenceSequence] = None
        
        # Parsed data
        self._positions: Optional[List[int]] = None
        self._adjusted_positions: Optional[List[int]] = None
        
        # Results
        self._result: Optional[AnalysisResult] = None
    
    @property
    def handler(self) -> AlignmentHandler:
        """Get the alignment handler, loading if necessary."""
        if self._handler is None:
            # Pass None or "auto" to enable format inference
            fmt = None if self.alignment_format in (None, "auto") else self.alignment_format
            self._handler = AlignmentHandler(
                self.alignment_path,
                fmt,
            )
        return self._handler
    
    @property
    def reference(self) -> ReferenceSequence:
        """Get the reference sequence."""
        if self._reference is None:
            self._reference = self.handler.get_reference(self.reference_id)
        return self._reference
    
    @property
    def positions(self) -> List[int]:
        """Get the parsed positions."""
        if self._positions is None:
            self._positions = self._parser.parse(self.positions_str)
        return self._positions
    
    @property
    def adjusted_positions(self) -> List[int]:
        """Get the gap-adjusted positions."""
        if self._adjusted_positions is None:
            self._adjusted_positions = self.handler.adjust_positions(
                self.positions,
                self.reference,
            )
        return self._adjusted_positions
    
    def validate(self) -> bool:
        """
        Validate all inputs before running analysis.
        
        Returns
        -------
        bool
            True if all validations pass.
            
        Raises
        ------
        FeaVarError
            If any validation fails.
        """
        logger.info("Running pre-flight validation...")
        
        # This will trigger loading and validation of each component
        _ = self.handler  # Load alignment
        _ = self.reference  # Find reference
        _ = self.positions  # Parse positions
        
        # Validate positions are within reference bounds
        self._parser.validate_positions(self.positions, self.reference.length)
        
        # Get adjusted positions
        _ = self.adjusted_positions
        
        # Validate no gaps at target positions
        if not self.handler.validate_positions_no_gaps(
            self.adjusted_positions, 
            self.reference
        ):
            logger.warning(
                "Some adjusted positions fall on gap characters in reference"
            )
        
        logger.info("Validation complete - all checks passed")
        return True
    
    def run(self) -> AnalysisResult:
        """
        Run the complete FeaVar analysis.
        
        Returns
        -------
        AnalysisResult
            The analysis results.
        """
        # Validate inputs
        self.validate()
        
        # Extract variants for all sequences
        logger.info("Extracting variant types...")
        variants = list(
            self.handler.extract_all_variants(self.adjusted_positions)
        )
        
        # Create variants DataFrame
        variants_df = pd.DataFrame(variants, columns=['accession', 'variant_type'])
        
        # Compute summary statistics
        summary_df = self._compute_summary(variants_df)
        
        # Store result
        self._result = AnalysisResult(
            variants_df=variants_df,
            summary_df=summary_df,
            reference_id=self.reference.identifier,
            positions=self.positions,
            adjusted_positions=self.adjusted_positions,
        )
        
        logger.info(
            f"Analysis complete: {self._result.num_sequences} sequences, "
            f"{self._result.num_variant_types} unique variant types"
        )
        
        return self._result
    
    def _compute_summary(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute variant type summary statistics.
        
        Parameters
        ----------
        variants_df : pd.DataFrame
            DataFrame with 'accession' and 'variant_type' columns.
            
        Returns
        -------
        pd.DataFrame
            Summary with 'variant_type', 'count', and 'VT' columns.
        """
        # Count occurrences of each variant type
        summary = (
            variants_df
            .groupby('variant_type')
            .size()
            .reset_index(name='count')
            .sort_values('count', ascending=False)
            .reset_index(drop=True)
        )
        
        # Generate VT names based on naming scheme
        if self.naming_scheme == NamingScheme.DELTA:
            # Get reference variant sequence
            ref_variant = get_reference_variant_sequence(
                variants_df, 
                self.reference_id
            )
            
            if ref_variant is None:
                # Fallback: use most common variant as reference
                ref_variant = summary.iloc[0]['variant_type']
                logger.warning(
                    f"Reference '{self.reference_id}' not found in variants. "
                    f"Using most common variant as reference for naming."
                )
            
            # Generate delta names for each variant
            name_map = compute_variant_names(
                variants_df,
                ref_variant,
                self.positions,  # Use original 1-based positions
                NamingScheme.DELTA,
            )
            summary['VT'] = summary['variant_type'].map(name_map)
            
            logger.info(f"Using delta naming scheme (reference variant: {ref_variant[:20]}...)")
            
        else:  # NamingScheme.RANKED
            # Traditional rank-based naming (VT-001, VT-002, etc.)
            summary['VT'] = [generate_ranked_name(i + 1) for i in range(len(summary))]
            logger.info("Using ranked naming scheme")
        
        return summary
    
    def merge_metadata(
        self,
        metadata_path: str,
        merge_column: str = 'accession',
        delimiter: str = '\t',
    ) -> pd.DataFrame:
        """
        Merge analysis results with metadata file.
        
        Parameters
        ----------
        metadata_path : str
            Path to metadata file (CSV or TSV).
        merge_column : str, optional
            Column to merge on (default: 'accession').
        delimiter : str, optional
            Delimiter for metadata file (default: tab).
            
        Returns
        -------
        pd.DataFrame
            Merged DataFrame with variants and metadata.
            
        Raises
        ------
        MetadataError
            If metadata file cannot be loaded or merged.
        """
        if self._result is None:
            raise FeaVarError("Must run analysis before merging metadata")
        
        try:
            metadata = pd.read_csv(metadata_path, delimiter=delimiter)
        except Exception as e:
            raise MetadataError(
                f"Failed to load metadata file '{metadata_path}': {e}"
            ) from e
        
        if merge_column not in metadata.columns:
            raise MetadataError(
                f"Merge column '{merge_column}' not found in metadata. "
                f"Available columns: {list(metadata.columns)}"
            )
        
        # Merge variants with metadata
        merged = self._result.variants_df.merge(
            metadata,
            on=merge_column,
            how='outer',
        )
        
        # Add VT labels from summary
        vt_map = dict(zip(
            self._result.summary_df['variant_type'],
            self._result.summary_df['VT'],
        ))
        merged['VT'] = merged['variant_type'].map(vt_map)
        
        logger.info(f"Merged {len(merged)} records with metadata")
        
        return merged
    
    def save_results(
        self,
        prefix: str = "feavar",
        save_variants: bool = True,
        save_summary: bool = True,
    ) -> Dict[str, Path]:
        """
        Save analysis results to CSV files.
        
        Parameters
        ----------
        prefix : str, optional
            Prefix for output filenames (default: "feavar").
        save_variants : bool, optional
            Whether to save the variants DataFrame (default: True).
        save_summary : bool, optional
            Whether to save the summary DataFrame (default: True).
            
        Returns
        -------
        Dict[str, Path]
            Dictionary mapping result type to file path.
        """
        if self._result is None:
            raise FeaVarError("Must run analysis before saving results")
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        saved_files = {}
        
        if save_variants:
            path = self.output_dir / f"{prefix}_variants.csv"
            self._result.variants_df.to_csv(path, index=False)
            saved_files['variants'] = path
            logger.info(f"Saved variants to {path}")
        
        if save_summary:
            path = self.output_dir / f"{prefix}_summary.csv"
            self._result.summary_df.to_csv(path, index=False)
            saved_files['summary'] = path
            logger.info(f"Saved summary to {path}")
        
        return saved_files
    
    def get_result(self) -> Optional[AnalysisResult]:
        """Get the analysis result, or None if not yet run."""
        return self._result
