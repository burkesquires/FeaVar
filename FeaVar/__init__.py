"""
FeaVar - Sequence Feature Variant Type Analysis

A Python package for computing clusters of unique substrings or "variant types"
(SFVTs) for user-selected sequence features in a set of aligned sequences.
"""

__version__ = "2.0.0"
__author__ = "R. Burke Squires"
__email__ = "burkesquires@gmail.com"

from feavar.alignment import AlignmentHandler, infer_format_from_extension
from feavar.core import FeaVarAnalysis
from feavar.naming import NamingScheme, generate_delta_name, generate_ranked_name
from feavar.parser import PositionParser
from feavar.visualization import VariantPlotter

__all__ = [
    "FeaVarAnalysis",
    "PositionParser",
    "AlignmentHandler",
    "VariantPlotter",
    "infer_format_from_extension",
    "NamingScheme",
    "generate_delta_name",
    "generate_ranked_name",
]
