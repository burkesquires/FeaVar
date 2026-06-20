"""
Tests for the position parser module.
"""

import pytest

from feavar.exceptions import InvalidPositionError, PositionParseError
from feavar.parser import PositionParser


class TestPositionParser:
    """Tests for PositionParser class."""

    @pytest.fixture
    def parser(self):
        """Create a parser instance."""
        return PositionParser()

    # Basic parsing tests

    def test_parse_single_position(self, parser):
        """Test parsing a single position."""
        assert parser.parse("10") == [10]
        assert parser.parse("100") == [100]
        assert parser.parse("1") == [1]

    def test_parse_comma_separated(self, parser):
        """Test parsing comma-separated positions."""
        assert parser.parse("10,20,30") == [10, 20, 30]
        assert parser.parse("10, 20, 30") == [10, 20, 30]
        assert parser.parse("10 ,20 ,30") == [10, 20, 30]

    def test_parse_range(self, parser):
        """Test parsing a range of positions."""
        assert parser.parse("10-15") == [10, 11, 12, 13, 14, 15]
        assert parser.parse("1-5") == [1, 2, 3, 4, 5]

    def test_parse_range_with_spaces(self, parser):
        """Test parsing ranges with spaces."""
        assert parser.parse("10 - 15") == [10, 11, 12, 13, 14, 15]
        assert parser.parse(" 10-15 ") == [10, 11, 12, 13, 14, 15]

    def test_parse_mixed(self, parser):
        """Test parsing mixed ranges and positions."""
        result = parser.parse("10-12,20,30-32")
        assert result == [10, 11, 12, 20, 30, 31, 32]

    def test_parse_removes_duplicates(self, parser):
        """Test that duplicates are removed."""
        assert parser.parse("10,10,10") == [10]
        assert parser.parse("10-12,11") == [10, 11, 12]

    def test_parse_sorts_output(self, parser):
        """Test that output is sorted."""
        assert parser.parse("30,10,20") == [10, 20, 30]
        assert parser.parse("100,1,50") == [1, 50, 100]

    # Edge cases

    def test_parse_single_element_range(self, parser):
        """Test parsing a range where start equals end."""
        assert parser.parse("10-10") == [10]

    def test_parse_large_positions(self, parser):
        """Test parsing large position numbers."""
        assert parser.parse("1000,2000,3000") == [1000, 2000, 3000]
        assert parser.parse("10000-10002") == [10000, 10001, 10002]

    # Error cases

    def test_parse_empty_string_raises(self, parser):
        """Test that empty string raises PositionParseError."""
        with pytest.raises(PositionParseError, match="cannot be empty"):
            parser.parse("")

        with pytest.raises(PositionParseError, match="cannot be empty"):
            parser.parse("   ")

    def test_parse_invalid_characters_raises(self, parser):
        """Test that invalid characters raise PositionParseError."""
        with pytest.raises(PositionParseError, match="Invalid characters"):
            parser.parse("10;20;30")

        with pytest.raises(PositionParseError, match="Invalid characters"):
            parser.parse("10.20.30")

        with pytest.raises(PositionParseError, match="Invalid characters"):
            parser.parse("abc")

    def test_parse_invalid_range_raises(self, parser):
        """Test that invalid ranges raise PositionParseError."""
        with pytest.raises(PositionParseError, match="Invalid range"):
            parser.parse("10-")

        with pytest.raises(PositionParseError, match="Invalid range"):
            parser.parse("-10")

    def test_parse_reversed_range_raises(self, parser):
        """Test that reversed ranges raise PositionParseError."""
        with pytest.raises(PositionParseError, match="cannot be greater"):
            parser.parse("20-10")

    def test_parse_zero_position_raises(self, parser):
        """Test that zero position raises PositionParseError."""
        with pytest.raises(PositionParseError, match="must be positive"):
            parser.parse("0")

        with pytest.raises(PositionParseError, match="must be positive"):
            parser.parse("0-10")

    def test_parse_negative_position_raises(self, parser):
        """Test that a negative position raises PositionParseError."""
        with pytest.raises(PositionParseError, match="Invalid range"):
            parser.parse("-5")  # This looks like a range start

    # Validation tests

    def test_validate_positions_valid(self, parser):
        """Test validation of valid positions."""
        assert parser.validate_positions([1, 5, 10], 100) is True
        assert parser.validate_positions([1], 1) is True
        assert parser.validate_positions([100], 100) is True

    def test_validate_positions_out_of_range(self, parser):
        """Test validation of out-of-range positions."""
        with pytest.raises(InvalidPositionError, match="exceed sequence length"):
            parser.validate_positions([1, 5, 150], 100)

        with pytest.raises(InvalidPositionError, match="exceed sequence length"):
            parser.validate_positions([101], 100)


class TestPositionParserIntegration:
    """Integration tests for PositionParser."""

    def test_real_world_position_string(self):
        """Test with a real-world position string format."""
        parser = PositionParser()

        # Example from FeaVar README
        result = parser.parse("124-142")
        assert len(result) == 19
        assert result[0] == 124
        assert result[-1] == 142

        # Complex example
        result = parser.parse("124,142,143,268,331,332")
        assert result == [124, 142, 143, 268, 331, 332]

        # Mixed
        result = parser.parse("124-130,142,143,268,331")
        assert 124 in result
        assert 130 in result
        assert 142 in result
        assert len(result) == 11  # 7 from range + 4 singles
