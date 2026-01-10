"""
Tests for the alignment handler module.
"""

import pytest
from pathlib import Path

from feavar.alignment import AlignmentHandler, ReferenceSequence
from feavar.exceptions import AlignmentError, ReferenceNotFoundError, InvalidPositionError


class TestReferenceSequence:
    """Tests for ReferenceSequence dataclass."""
    
    def test_reference_sequence_creation(self):
        """Test creating a ReferenceSequence."""
        ref = ReferenceSequence(
            identifier="SEQ001",
            sequence="--ABC--DEF",
            raw_sequence="ABCDEF",
            gap_positions=[0, 1, 5, 6],
        )
        
        assert ref.identifier == "SEQ001"
        assert ref.length == 6
        assert ref.aligned_length == 10
        assert ref.gap_positions == [0, 1, 5, 6]
    
    def test_reference_sequence_no_gaps(self):
        """Test ReferenceSequence with no gaps."""
        ref = ReferenceSequence(
            identifier="SEQ001",
            sequence="ABCDEF",
            raw_sequence="ABCDEF",
            gap_positions=[],
        )
        
        assert ref.length == 6
        assert ref.aligned_length == 6


class TestAlignmentHandler:
    """Tests for AlignmentHandler class."""
    
    def test_load_alignment(self, sample_alignment_file):
        """Test loading an alignment file."""
        handler = AlignmentHandler(sample_alignment_file)
        
        assert handler.num_sequences == 3
        assert handler.alignment_length == 120
    
    def test_load_nonexistent_file_raises(self, tmp_path):
        """Test that loading nonexistent file raises AlignmentError."""
        with pytest.raises(AlignmentError, match="not found"):
            AlignmentHandler(tmp_path / "nonexistent.clw")
    
    def test_get_sequence_ids(self, sample_alignment_file):
        """Test getting sequence IDs."""
        handler = AlignmentHandler(sample_alignment_file)
        ids = handler.get_sequence_ids()
        
        assert "SEQ001" in ids
        assert "SEQ002" in ids
        assert "SEQ003" in ids
        assert len(ids) == 3
    
    def test_find_reference_exact(self, sample_alignment_file):
        """Test finding reference by exact match."""
        handler = AlignmentHandler(sample_alignment_file)
        record = handler.find_reference("SEQ001")
        
        assert record is not None
        assert record.id == "SEQ001"
    
    def test_find_reference_partial(self, sample_alignment_file):
        """Test finding reference by partial match."""
        handler = AlignmentHandler(sample_alignment_file)
        record = handler.find_reference("001")  # Partial match
        
        assert record is not None
        assert "001" in record.id
    
    def test_find_reference_not_found(self, sample_alignment_file):
        """Test that missing reference returns None."""
        handler = AlignmentHandler(sample_alignment_file)
        record = handler.find_reference("NONEXISTENT")
        
        assert record is None
    
    def test_get_reference(self, sample_alignment_file):
        """Test getting reference sequence object."""
        handler = AlignmentHandler(sample_alignment_file)
        ref = handler.get_reference("SEQ002")
        
        assert isinstance(ref, ReferenceSequence)
        assert ref.identifier == "SEQ002"
        assert len(ref.gap_positions) > 0  # SEQ002 has gaps
    
    def test_get_reference_not_found_raises(self, sample_alignment_file):
        """Test that missing reference raises error."""
        handler = AlignmentHandler(sample_alignment_file)
        
        with pytest.raises(ReferenceNotFoundError, match="not found"):
            handler.get_reference("NONEXISTENT")


class TestAlignmentHandlerPositionAdjustment:
    """Tests for position adjustment functionality."""
    
    def test_build_position_map_no_gaps(self, sample_alignment_file):
        """Test building position map with no gaps."""
        handler = AlignmentHandler(sample_alignment_file)
        ref = handler.get_reference("SEQ001")  # No gaps
        
        position_map = handler.build_position_map(ref)
        
        # With no gaps, positions should map 1:1
        assert position_map[1] == 1
        assert position_map[10] == 10
    
    def test_build_position_map_with_gaps(self, complex_alignment_file):
        """Test building position map with gaps."""
        handler = AlignmentHandler(complex_alignment_file)
        ref = handler.get_reference("REFERENCE")
        
        position_map = handler.build_position_map(ref)
        
        # Sequence: --ABC--DEF---GHI----JKL
        # Position 1 (A) should map to alignment position 3
        assert position_map[1] == 3
        # Position 4 (D) should map to alignment position 8
        assert position_map[4] == 8
    
    def test_adjust_positions(self, complex_alignment_file):
        """Test adjusting positions for gaps."""
        handler = AlignmentHandler(complex_alignment_file)
        ref = handler.get_reference("REFERENCE")
        
        adjusted = handler.adjust_positions([1, 2, 3], ref)
        
        # Positions should be adjusted for leading gaps
        assert adjusted[0] > 1  # First position should be shifted
    
    def test_adjust_positions_invalid_raises(self, sample_alignment_file):
        """Test that invalid positions raise error."""
        handler = AlignmentHandler(sample_alignment_file)
        ref = handler.get_reference("SEQ001")
        
        with pytest.raises(InvalidPositionError, match="out of range"):
            handler.adjust_positions([1, 5, 10000], ref)  # 10000 is too large


class TestAlignmentHandlerVariantExtraction:
    """Tests for variant extraction functionality."""
    
    def test_extract_variant(self, sample_alignment_file):
        """Test extracting variant from a sequence."""
        handler = AlignmentHandler(sample_alignment_file)
        record = handler.find_reference("SEQ001")
        
        # Extract first 5 characters
        variant = handler.extract_variant(record, [1, 2, 3, 4, 5])
        
        assert len(variant) == 5
        assert variant == "MFQAF"
    
    def test_extract_all_variants(self, sample_alignment_file):
        """Test extracting variants from all sequences."""
        handler = AlignmentHandler(sample_alignment_file)
        
        variants = list(handler.extract_all_variants([1, 2, 3, 4, 5]))
        
        assert len(variants) == 3
        # Each variant is (id, variant_type)
        assert all(len(v) == 2 for v in variants)
        # All variants should have length 5
        assert all(len(v[1]) == 5 for v in variants)
    
    def test_extract_variant_with_gaps(self, sample_alignment_file):
        """Test extracting variant that includes gap characters."""
        handler = AlignmentHandler(sample_alignment_file)
        ref = handler.get_reference("SEQ002")
        
        # SEQ002 starts with gaps
        variant = handler.extract_variant(
            handler.find_reference("SEQ002"),
            [1, 2, 3, 4, 5]
        )
        
        assert "-" in variant  # Should contain gap characters
    
    def test_validate_positions_no_gaps(self, sample_alignment_file):
        """Test validation of positions not on gaps."""
        handler = AlignmentHandler(sample_alignment_file)
        ref = handler.get_reference("SEQ001")  # No leading gaps
        
        # Positions in ungapped region should be valid
        assert handler.validate_positions_no_gaps([6, 7, 8, 9, 10], ref) is True
    
    def test_validate_positions_on_gaps(self, sample_alignment_file):
        """Test validation fails for positions on gaps."""
        handler = AlignmentHandler(sample_alignment_file)
        ref = handler.get_reference("SEQ002")  # Has leading gaps
        
        # Position 1 falls on a gap
        assert handler.validate_positions_no_gaps([1], ref) is False
