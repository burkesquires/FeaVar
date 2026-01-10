"""
Custom exceptions for FeaVar package.
"""


class FeaVarError(Exception):
    """Base exception for FeaVar package."""
    pass


class PositionParseError(FeaVarError):
    """Raised when position string cannot be parsed."""
    pass


class AlignmentError(FeaVarError):
    """Raised when there's an issue with the alignment file."""
    pass


class ReferenceNotFoundError(AlignmentError):
    """Raised when reference sequence is not found in alignment."""
    pass


class InvalidPositionError(FeaVarError):
    """Raised when positions are invalid for the given sequence."""
    pass


class MetadataError(FeaVarError):
    """Raised when there's an issue with metadata file."""
    pass
