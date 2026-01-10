"""
Position parsing utilities for FeaVar.

Handles parsing of position strings like "100-110" or "100-110,120,130"
into lists of integer positions.
"""

import re
import logging
from typing import List, Optional

from feavar.exceptions import PositionParseError

logger = logging.getLogger(__name__)


class PositionParser:
    """
    Parser for sequence position strings.
    
    Supports formats:
    - Single positions: "100"
    - Comma-separated: "100,110,120"
    - Ranges: "100-110"
    - Mixed: "100-110,120,130-140"
    
    Examples
    --------
    >>> parser = PositionParser()
    >>> parser.parse("100-105")
    [100, 101, 102, 103, 104, 105]
    >>> parser.parse("100,110,120")
    [100, 110, 120]
    >>> parser.parse("100-102,110")
    [100, 101, 102, 110]
    """
    
    # Pattern to validate input: only digits, commas, dashes, and whitespace
    VALID_PATTERN = re.compile(r'^[\d,\-\s]+$')
    
    def __init__(self):
        """Initialize the position parser."""
        pass
    
    def parse(self, raw_positions: str) -> List[int]:
        """
        Parse a position string into a sorted list of integers.
        
        Parameters
        ----------
        raw_positions : str
            Position string in format "100-110" or "100,110,120" or mixed.
            
        Returns
        -------
        List[int]
            Sorted list of position integers (1-based).
            
        Raises
        ------
        PositionParseError
            If the position string is empty, malformed, or contains invalid values.
            
        Examples
        --------
        >>> parser = PositionParser()
        >>> parser.parse("10-15")
        [10, 11, 12, 13, 14, 15]
        """
        if not raw_positions or not raw_positions.strip():
            raise PositionParseError("Position string cannot be empty")
        
        # Clean whitespace
        cleaned = raw_positions.strip()
        
        # Validate characters
        if not self.VALID_PATTERN.match(cleaned):
            raise PositionParseError(
                f"Invalid characters in position string: '{raw_positions}'. "
                "Only digits, commas, dashes, and spaces are allowed."
            )
        
        position_coordinates = []
        
        # Split by comma
        position_groupings = cleaned.split(",")
        
        for group in position_groupings:
            group = group.strip()
            if not group:
                continue
                
            positions = self._parse_group(group)
            position_coordinates.extend(positions)
        
        if not position_coordinates:
            raise PositionParseError(
                f"No valid positions found in: '{raw_positions}'"
            )
        
        # Sort and remove duplicates
        result = sorted(set(position_coordinates))
        
        logger.debug(f"Parsed positions '{raw_positions}' -> {result}")
        
        return result
    
    def _parse_group(self, group: str) -> List[int]:
        """
        Parse a single group (either a range or single position).
        
        Parameters
        ----------
        group : str
            A position group like "100" or "100-110".
            
        Returns
        -------
        List[int]
            List of positions from this group.
        """
        if "-" in group:
            return self._parse_range(group)
        else:
            return [self._parse_single(group)]
    
    def _parse_range(self, range_str: str) -> List[int]:
        """
        Parse a range string like "100-110".
        
        Parameters
        ----------
        range_str : str
            Range string in format "start-end".
            
        Returns
        -------
        List[int]
            List of integers from start to end (inclusive).
        """
        parts = range_str.split("-")
        
        # Filter out empty parts (handles cases like "-100" or "100-")
        parts = [p.strip() for p in parts if p.strip()]
        
        if len(parts) != 2:
            raise PositionParseError(
                f"Invalid range format: '{range_str}'. Expected format: 'start-end'"
            )
        
        try:
            start = int(parts[0])
            end = int(parts[1])
        except ValueError as e:
            raise PositionParseError(
                f"Invalid numbers in range: '{range_str}'"
            ) from e
        
        if start <= 0 or end <= 0:
            raise PositionParseError(
                f"Positions must be positive (1-based): '{range_str}'"
            )
        
        if start > end:
            raise PositionParseError(
                f"Range start ({start}) cannot be greater than end ({end})"
            )
        
        return list(range(start, end + 1))
    
    def _parse_single(self, pos_str: str) -> int:
        """
        Parse a single position string.
        
        Parameters
        ----------
        pos_str : str
            Single position string like "100".
            
        Returns
        -------
        int
            The position as an integer.
        """
        try:
            pos = int(pos_str.strip())
        except ValueError as e:
            raise PositionParseError(
                f"Invalid position value: '{pos_str}'"
            ) from e
        
        if pos <= 0:
            raise PositionParseError(
                f"Position must be positive (1-based): {pos}"
            )
        
        return pos
    
    def validate_positions(
        self, 
        positions: List[int], 
        sequence_length: int
    ) -> bool:
        """
        Validate that all positions are within sequence bounds.
        
        Parameters
        ----------
        positions : List[int]
            List of 1-based positions.
        sequence_length : int
            Length of the sequence (excluding gaps).
            
        Returns
        -------
        bool
            True if all positions are valid.
            
        Raises
        ------
        InvalidPositionError
            If any position is out of bounds.
        """
        from feavar.exceptions import InvalidPositionError
        
        invalid = [p for p in positions if p > sequence_length]
        
        if invalid:
            raise InvalidPositionError(
                f"Positions {invalid} exceed sequence length ({sequence_length})"
            )
        
        return True
