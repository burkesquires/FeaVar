"""
Alignment handling utilities for FeaVar.

Provides classes for loading, validating, and processing multiple sequence alignments.
"""

import logging
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from feavar.exceptions import (
    AlignmentError,
    InvalidPositionError,
    ReferenceNotFoundError,
)

logger = logging.getLogger(__name__)


# Mapping of file extensions to BioPython format names
EXTENSION_FORMAT_MAP = {
    # Clustal formats
    ".clw": "clustal",
    ".clustal": "clustal",
    ".aln": "clustal",
    # FASTA formats
    ".fasta": "fasta",
    ".fa": "fasta",
    ".fas": "fasta",
    ".fna": "fasta",
    ".faa": "fasta",
    ".ffn": "fasta",
    # Phylip formats
    ".phy": "phylip",
    ".phylip": "phylip",
    ".ph": "phylip",
    # Nexus formats
    ".nex": "nexus",
    ".nexus": "nexus",
    ".nxs": "nexus",
    # Stockholm format
    ".sto": "stockholm",
    ".stockholm": "stockholm",
    ".stk": "stockholm",
    # EMBOSS formats
    ".msf": "msf",
    # MAF format
    ".maf": "maf",
}


def infer_format_from_extension(filepath: str | Path) -> str | None:
    """
    Infer alignment format from file extension.

    Parameters
    ----------
    filepath : str or Path
        Path to the alignment file.

    Returns
    -------
    str or None
        The inferred format name, or None if unknown.

    Examples
    --------
    >>> infer_format_from_extension("sequences.clw")
    'clustal'
    >>> infer_format_from_extension("sequences.fasta")
    'fasta'
    >>> infer_format_from_extension("sequences.xyz")
    None
    """
    path = Path(filepath)
    ext = path.suffix.lower()
    return EXTENSION_FORMAT_MAP.get(ext)


@dataclass
class ReferenceSequence:
    """
    Container for reference sequence information.

    Attributes
    ----------
    identifier : str
        The sequence identifier.
    sequence : str
        The aligned sequence (may contain gaps).
    raw_sequence : str
        The sequence with gaps removed.
    gap_positions : List[int]
        0-based indices where gaps occur in the aligned sequence.
    """

    identifier: str
    sequence: str
    raw_sequence: str
    gap_positions: list[int]

    @property
    def length(self) -> int:
        """Length of the raw (ungapped) sequence."""
        return len(self.raw_sequence)

    @property
    def aligned_length(self) -> int:
        """Length of the aligned sequence (including gaps)."""
        return len(self.sequence)


class AlignmentHandler:
    """
    Handler for multiple sequence alignment files.

    Provides methods for loading alignments, finding reference sequences,
    and extracting sequence features.

    Parameters
    ----------
    alignment_path : str or Path
        Path to the alignment file.
    alignment_format : str, optional
        Format of the alignment file. If None or "auto", will attempt to
        infer from file extension (default: None).
        Supports any format recognized by BioPython's AlignIO.

    Examples
    --------
    >>> handler = AlignmentHandler("sequences.clw")  # Auto-detects clustal
    >>> handler = AlignmentHandler("sequences.fasta")  # Auto-detects fasta
    >>> handler = AlignmentHandler("sequences.txt", "clustal")  # Explicit format
    >>> ref = handler.get_reference("CY021716")
    >>> print(ref.length)
    338
    """

    SUPPORTED_FORMATS = ["clustal", "fasta", "phylip", "nexus", "stockholm"]

    def __init__(self, alignment_path: str, alignment_format: str | None = None):
        """Initialize the alignment handler."""
        self.alignment_path = Path(alignment_path)

        # Infer format from extension if not provided or if "auto"
        if alignment_format is None or alignment_format.lower() == "auto":
            inferred = infer_format_from_extension(self.alignment_path)
            if inferred:
                self.alignment_format = inferred
                logger.info(
                    f"Inferred format '{inferred}' from extension '{self.alignment_path.suffix}'"
                )
            else:
                # Default to clustal if can't infer
                self.alignment_format = "clustal"
                logger.warning(
                    f"Could not infer format from extension '{self.alignment_path.suffix}'. "
                    f"Defaulting to 'clustal'. Use -f option to specify format."
                )
        else:
            self.alignment_format = alignment_format.lower()

        self._alignment: MultipleSeqAlignment | None = None
        self._reference: ReferenceSequence | None = None

        self._validate_format()
        self._load_alignment()

    def _validate_format(self) -> None:
        """Validate the alignment format is supported."""
        if self.alignment_format not in self.SUPPORTED_FORMATS:
            logger.warning(
                f"Format '{self.alignment_format}' not in known formats. "
                f"Supported: {self.SUPPORTED_FORMATS}. Attempting anyway..."
            )

    def _load_alignment(self) -> None:
        """Load the alignment file."""
        if not self.alignment_path.exists():
            raise AlignmentError(f"Alignment file not found: {self.alignment_path}")

        try:
            self._alignment = AlignIO.read(
                str(self.alignment_path), self.alignment_format
            )
            identifier_counts = Counter(record.id for record in self._alignment)
            duplicate_ids = [
                identifier
                for identifier, count in identifier_counts.items()
                if count > 1
            ]
            if duplicate_ids:
                examples = ", ".join(sorted(duplicate_ids)[:5])
                raise AlignmentError(
                    "Alignment sequence identifiers must be unique; "
                    f"duplicates include: {examples}"
                )
            logger.info(
                f"Loaded alignment with {len(self._alignment)} sequences "
                f"of length {self._alignment.get_alignment_length()}"
            )
        except FileNotFoundError as e:
            raise AlignmentError(
                f"Alignment file not found: {self.alignment_path}"
            ) from e
        except ValueError as e:
            error_msg = str(e) if str(e) else "Unknown parsing error"
            raise AlignmentError(
                f"Failed to parse alignment file '{self.alignment_path}' "
                f"as {self.alignment_format} format: {error_msg}\n"
                f"Hint: Check that the file format matches the -f/--format option. "
                f"Try -f fasta if this is a FASTA file."
            ) from e
        except Exception as e:
            error_msg = str(e) if str(e) else f"Unknown error ({type(e).__name__})"
            # Log the full exception for debugging
            logger.exception(f"Error loading alignment: {error_msg}")
            raise AlignmentError(
                f"Failed to load alignment file '{self.alignment_path}': {error_msg}"
            ) from e

    @property
    def alignment(self) -> MultipleSeqAlignment:
        """Get the loaded alignment."""
        if self._alignment is None:
            raise AlignmentError("Alignment not loaded")
        return self._alignment

    @property
    def num_sequences(self) -> int:
        """Number of sequences in the alignment."""
        return len(self.alignment)

    @property
    def alignment_length(self) -> int:
        """Length of the alignment (including gaps)."""
        return int(self.alignment.get_alignment_length())

    def get_sequence_ids(self) -> list[str]:
        """Get list of all sequence identifiers in the alignment."""
        return [record.id for record in self.alignment]

    def find_reference(self, identifier: str) -> SeqRecord | None:
        """
        Find a sequence by identifier (exact match preferred).

        Parameters
        ----------
        identifier : str
            Full or partial sequence identifier to search for.

        Returns
        -------
        SeqRecord or None
            The matching sequence record, or None if not found.

        Raises
        ------
        ReferenceNotFoundError
            If the identifier is ambiguous and matches multiple sequences.
        """
        records = list(self.alignment)

        # Prefer exact match
        exact_matches = [record for record in records if record.id == identifier]
        if len(exact_matches) == 1:
            return cast(SeqRecord, exact_matches[0])
        if len(exact_matches) > 1:
            available = [record.id for record in exact_matches[:5]]
            raise ReferenceNotFoundError(
                f"Reference sequence '{identifier}' is ambiguous; "
                f"found {len(exact_matches)} exact matches (e.g., {available}...)"
            )

        # Fall back to partial match
        partial_matches = [record for record in records if identifier in record.id]
        if len(partial_matches) == 1:
            return cast(SeqRecord, partial_matches[0])
        if len(partial_matches) > 1:
            available = [record.id for record in partial_matches[:5]]
            raise ReferenceNotFoundError(
                f"Reference sequence '{identifier}' is ambiguous; "
                f"found {len(partial_matches)} partial matches (e.g., {available}...)"
            )

        return None

    def get_reference(self, identifier: str) -> ReferenceSequence:
        """
        Get reference sequence information.

        Parameters
        ----------
        identifier : str
            The reference sequence identifier.

        Returns
        -------
        ReferenceSequence
            Reference sequence data container.

        Raises
        ------
        ReferenceNotFoundError
            If the reference sequence is not found in the alignment.
        """
        record = self.find_reference(identifier)

        if record is None:
            available = self.get_sequence_ids()[:5]
            raise ReferenceNotFoundError(
                f"Reference sequence '{identifier}' not found in alignment. "
                f"Available sequences include: {available}..."
            )

        sequence = str(record.seq)
        gap_positions = [i for i, char in enumerate(sequence) if char == "-"]
        raw_sequence = sequence.replace("-", "")

        self._reference = ReferenceSequence(
            identifier=str(record.id),
            sequence=sequence,
            raw_sequence=raw_sequence,
            gap_positions=gap_positions,
        )

        logger.info(
            f"Found reference '{record.id}' with {len(raw_sequence)} residues "
            f"({len(gap_positions)} gaps in alignment)"
        )

        return self._reference

    def build_position_map(self, reference: ReferenceSequence) -> dict[int, int]:
        """
        Build a mapping from sequence positions to alignment positions.

        This accounts for gaps in the reference sequence.

        Parameters
        ----------
        reference : ReferenceSequence
            The reference sequence.

        Returns
        -------
        Dict[int, int]
            Mapping of 1-based sequence position to 1-based alignment position.

        Examples
        --------
        For sequence "--TTA--GGG" (gaps at positions 0,1,5,6):
        - Sequence position 1 (T) -> Alignment position 3
        - Sequence position 2 (T) -> Alignment position 4
        - Sequence position 3 (A) -> Alignment position 5
        """
        position_map = {}
        seq_pos = 0  # 0-based sequence position (excluding gaps)

        for align_pos, char in enumerate(reference.sequence):
            if char != "-":
                seq_pos += 1
                # Store as 1-based positions
                position_map[seq_pos] = align_pos + 1

        return position_map

    def adjust_positions(
        self,
        positions: list[int],
        reference: ReferenceSequence,
    ) -> list[int]:
        """
        Adjust sequence positions to alignment positions.

        Parameters
        ----------
        positions : List[int]
            1-based positions in the ungapped reference sequence.
        reference : ReferenceSequence
            The reference sequence.

        Returns
        -------
        List[int]
            1-based positions in the aligned sequence.

        Raises
        ------
        InvalidPositionError
            If any position is out of bounds.
        """
        # Validate positions
        max_pos = reference.length
        invalid = [p for p in positions if p > max_pos or p < 1]
        if invalid:
            raise InvalidPositionError(
                f"Positions {invalid} are out of range [1, {max_pos}]"
            )

        position_map = self.build_position_map(reference)
        adjusted = [position_map[p] for p in positions]

        logger.debug(f"Adjusted positions: {positions} -> {adjusted}")

        return adjusted

    def validate_positions_no_gaps(
        self,
        aligned_positions: list[int],
        reference: ReferenceSequence,
    ) -> bool:
        """
        Validate that adjusted positions don't fall on gaps.

        Parameters
        ----------
        aligned_positions : List[int]
            1-based positions in the aligned sequence.
        reference : ReferenceSequence
            The reference sequence.

        Returns
        -------
        bool
            True if no positions fall on gaps.
        """
        # Convert to 0-based for checking
        for pos in aligned_positions:
            idx = pos - 1
            if idx >= len(reference.sequence):
                return False
            if reference.sequence[idx] == "-":
                return False
        return True

    def extract_variant(
        self,
        record: SeqRecord,
        aligned_positions: list[int],
    ) -> str:
        """
        Extract the variant type (substring) at specified positions.

        Parameters
        ----------
        record : SeqRecord
            The sequence record.
        aligned_positions : List[int]
            1-based positions in the alignment.

        Returns
        -------
        str
            The extracted substring (variant type).
        """
        sequence = str(record.seq)
        # Convert to 0-based indexing
        variant = "".join(sequence[pos - 1] for pos in aligned_positions)
        return variant

    def extract_all_variants(
        self,
        aligned_positions: list[int],
    ) -> Iterator[tuple[str, str]]:
        """
        Extract variants for all sequences in the alignment.

        Parameters
        ----------
        aligned_positions : List[int]
            1-based positions in the alignment.

        Yields
        ------
        Tuple[str, str]
            Pairs of (sequence_id, variant_type).
        """
        for record in self.alignment:
            variant = self.extract_variant(record, aligned_positions)
            yield (record.id, variant)
