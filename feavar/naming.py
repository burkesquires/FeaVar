"""
Variant naming utilities for FeaVar.

Provides different naming schemes for variant types.
"""

import logging
from enum import Enum

import pandas as pd

logger = logging.getLogger(__name__)


class NamingScheme(Enum):
    """Available variant naming schemes."""

    RANKED = "ranked"  # Traditional: VT-001, VT-002 based on frequency
    DELTA = "delta"  # Delta notation: VT-5G.12C based on differences from reference


def generate_delta_name(
    variant_sequence: str,
    reference_sequence: str,
    positions: list[int] | None = None,
) -> str:
    """
    Generate a stable variant name based on differences from the reference.

    Uses delta notation showing only positions and new residues that differ
    from the reference sequence.

    Parameters
    ----------
    variant_sequence : str
        The variant sequence to name.
    reference_sequence : str
        The reference sequence to compare against.
    positions : List[int], optional
        The 1-based positions being analyzed. If provided, position numbers
        in the name will correspond to these positions. If None, positions
        are numbered 1 to len(sequence).

    Returns
    -------
    str
        The delta-notation variant name (e.g., "VT-5G.12C.15A" or "VT-REF").

    Examples
    --------
    >>> generate_delta_name("GAAGACAGT", "GAAGACAGG")
    'VT-9T'
    >>> generate_delta_name("GAAGGCAGT", "GAAGACAGG")
    'VT-5G.9T'
    >>> generate_delta_name("GAAGACAGG", "GAAGACAGG")
    'VT-REF'
    """
    if len(variant_sequence) != len(reference_sequence):
        raise ValueError(
            f"Sequence lengths must match: variant={len(variant_sequence)}, "
            f"reference={len(reference_sequence)}"
        )

    # If positions not provided, use 1-based indexing
    if positions is None:
        positions = list(range(1, len(variant_sequence) + 1))

    if len(positions) != len(variant_sequence):
        raise ValueError(
            f"Positions length ({len(positions)}) must match sequence length "
            f"({len(variant_sequence)})"
        )

    # Find differences
    differences = []
    for i, (var_res, ref_res) in enumerate(
        zip(variant_sequence, reference_sequence, strict=True)
    ):
        if var_res != ref_res:
            # Use the actual position from the positions list (1-based)
            pos = positions[i]
            differences.append((pos, var_res))

    # Generate name
    if not differences:
        return "VT-REF"

    # Sort by position (should already be sorted, but ensure it)
    differences.sort(key=lambda x: x[0])

    # Build name: position + residue, dot-separated
    delta_parts = [f"{pos}{residue}" for pos, residue in differences]
    return f"VT-{'.'.join(delta_parts)}"


def generate_ranked_name(rank: int) -> str:
    """
    Generate a traditional ranked variant name.

    Parameters
    ----------
    rank : int
        The rank of the variant (1 = most common).

    Returns
    -------
    str
        The ranked variant name (e.g., "VT-001").

    Examples
    --------
    >>> generate_ranked_name(1)
    'VT-001'
    >>> generate_ranked_name(42)
    'VT-042'
    """
    return f"VT-{rank:03d}"


def compute_variant_names(
    variants_df: pd.DataFrame,
    reference_variant: str,
    positions: list[int],
    scheme: NamingScheme = NamingScheme.DELTA,
) -> dict[str, str]:
    """
    Compute names for all unique variants.

    Parameters
    ----------
    variants_df : pd.DataFrame
        DataFrame with 'variant_type' column containing variant sequences.
    reference_variant : str
        The reference variant sequence.
    positions : List[int]
        The 1-based positions being analyzed.
    scheme : NamingScheme
        The naming scheme to use.

    Returns
    -------
    Dict[str, str]
        Mapping from variant sequence to variant name.

    Examples
    --------
    >>> names = compute_variant_names(df, "GAAGACAGG", [1,2,3,4,5,6,7,8,9])
    >>> names["GAAGACAGG"]
    'VT-REF'
    >>> names["GAAGACAGT"]
    'VT-9T'
    """
    # Get unique variants with counts for ranking
    variant_counts = variants_df["variant_type"].value_counts()
    unique_variants = sorted(
        (str(variant) for variant in variant_counts.index.tolist()),
        key=lambda variant: (-int(variant_counts[variant]), variant),
    )

    name_map = {}

    if scheme == NamingScheme.DELTA:
        for variant in unique_variants:
            name_map[variant] = generate_delta_name(
                variant, reference_variant, positions
            )

    elif scheme == NamingScheme.RANKED:
        for rank, variant in enumerate(unique_variants, start=1):
            name_map[variant] = generate_ranked_name(rank)

    else:
        raise ValueError(f"Unknown naming scheme: {scheme}")

    logger.info(f"Generated {len(name_map)} variant names using {scheme.value} scheme")

    return name_map


def get_reference_variant_sequence(
    variants_df: pd.DataFrame,
    reference_accession: str,
) -> str | None:
    """
    Get the variant sequence for the reference accession.

    Parameters
    ----------
    variants_df : pd.DataFrame
        DataFrame with 'accession' and 'variant_type' columns.
    reference_accession : str
        The reference accession identifier.

    Returns
    -------
    str or None
        The variant sequence for the reference, or None if not found.
    """
    # Try exact match first
    match = variants_df[variants_df["accession"] == reference_accession]

    if match.empty:
        # Try partial match (reference_accession might be substring)
        match = variants_df[
            variants_df["accession"].str.contains(
                reference_accession,
                na=False,
                regex=False,
            )
        ]

    if match.empty:
        logger.warning(
            f"Reference accession '{reference_accession}' not found in variants"
        )
        return None

    return str(match.iloc[0]["variant_type"])
