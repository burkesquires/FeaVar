"""Tests for the naming module."""

import pytest
import pandas as pd

from feavar.naming import (
    NamingScheme,
    generate_delta_name,
    generate_ranked_name,
    compute_variant_names,
    get_reference_variant_sequence,
)


class TestGenerateDeltaName:
    """Tests for generate_delta_name function."""
    
    def test_identical_sequences_returns_ref(self):
        """Identical sequences should return VT-REF."""
        result = generate_delta_name("GAAGACAGG", "GAAGACAGG")
        assert result == "VT-REF"
    
    def test_single_difference_at_end(self):
        """Single difference at last position."""
        result = generate_delta_name("GAAGACAGT", "GAAGACAGG")
        assert result == "VT-9T"
    
    def test_single_difference_at_start(self):
        """Single difference at first position."""
        result = generate_delta_name("TAAGACAGG", "GAAGACAGG")
        assert result == "VT-1T"
    
    def test_multiple_differences(self):
        """Multiple differences should be dot-separated."""
        result = generate_delta_name("TAAGGCAGT", "GAAGACAGG")
        # Differences at: 1(G->T), 5(A->G), 9(G->T)
        assert result == "VT-1T.5G.9T"
    
    def test_all_different(self):
        """All positions different."""
        result = generate_delta_name("TTTT", "AAAA")
        assert result == "VT-1T.2T.3T.4T"
    
    def test_with_custom_positions(self):
        """Custom positions should be used in naming."""
        # Positions 100-108 instead of 1-9
        positions = [100, 101, 102, 103, 104, 105, 106, 107, 108]
        result = generate_delta_name("GAAGACAGT", "GAAGACAGG", positions)
        assert result == "VT-108T"
    
    def test_with_custom_positions_multiple_diffs(self):
        """Custom positions with multiple differences."""
        positions = [10, 20, 30, 40]
        result = generate_delta_name("TCGT", "ACGA", positions)
        # Differences at: 10(A->T), 40(A->T)
        assert result == "VT-10T.40T"
    
    def test_length_mismatch_raises_error(self):
        """Mismatched sequence lengths should raise ValueError."""
        with pytest.raises(ValueError, match="lengths must match"):
            generate_delta_name("GAAG", "GAAGACAGG")
    
    def test_positions_length_mismatch_raises_error(self):
        """Mismatched positions length should raise ValueError."""
        with pytest.raises(ValueError, match="Positions length"):
            generate_delta_name("GAAG", "GAAG", [1, 2, 3])  # 4 chars, 3 positions
    
    def test_positions_sorted_in_output(self):
        """Output should have positions sorted numerically."""
        # Even if differences found in different order internally
        positions = [50, 10, 30, 20]  # Unsorted positions
        result = generate_delta_name("TTTT", "AAAA", positions)
        # Should be sorted: 10T.20T.30T.50T
        assert result == "VT-10T.20T.30T.50T"


class TestGenerateRankedName:
    """Tests for generate_ranked_name function."""
    
    def test_rank_1(self):
        """Rank 1 should be VT-001."""
        assert generate_ranked_name(1) == "VT-001"
    
    def test_rank_42(self):
        """Rank 42 should be VT-042."""
        assert generate_ranked_name(42) == "VT-042"
    
    def test_rank_100(self):
        """Rank 100 should be VT-100."""
        assert generate_ranked_name(100) == "VT-100"
    
    def test_rank_999(self):
        """Rank 999 should be VT-999."""
        assert generate_ranked_name(999) == "VT-999"
    
    def test_rank_1000(self):
        """Rank 1000 should be VT-1000 (4 digits)."""
        assert generate_ranked_name(1000) == "VT-1000"


class TestComputeVariantNames:
    """Tests for compute_variant_names function."""
    
    @pytest.fixture
    def sample_variants_df(self):
        """Sample variants DataFrame."""
        return pd.DataFrame({
            'accession': ['seq1', 'seq2', 'seq3', 'seq4', 'seq5'],
            'variant_type': ['GAAG', 'GAAG', 'GAAT', 'GAAG', 'TAAG'],
        })
    
    def test_delta_scheme(self, sample_variants_df):
        """Delta scheme should generate proper names."""
        ref_variant = "GAAG"
        positions = [1, 2, 3, 4]
        
        names = compute_variant_names(
            sample_variants_df,
            ref_variant,
            positions,
            NamingScheme.DELTA,
        )
        
        assert names["GAAG"] == "VT-REF"
        assert names["GAAT"] == "VT-4T"
        assert names["TAAG"] == "VT-1T"
    
    def test_ranked_scheme(self, sample_variants_df):
        """Ranked scheme should generate proper names."""
        ref_variant = "GAAG"
        positions = [1, 2, 3, 4]
        
        names = compute_variant_names(
            sample_variants_df,
            ref_variant,
            positions,
            NamingScheme.RANKED,
        )
        
        # GAAG appears 3 times (most common) -> VT-001
        # GAAT appears 1 time -> VT-002 or VT-003
        # TAAG appears 1 time -> VT-002 or VT-003
        assert names["GAAG"] == "VT-001"
        assert names["GAAT"] in ["VT-002", "VT-003"]
        assert names["TAAG"] in ["VT-002", "VT-003"]


class TestGetReferenceVariantSequence:
    """Tests for get_reference_variant_sequence function."""
    
    @pytest.fixture
    def sample_df(self):
        """Sample DataFrame with accessions and variants."""
        return pd.DataFrame({
            'accession': ['CY021716', 'CY020292', 'CY083917'],
            'variant_type': ['GAAGACAGG', 'GAAGACAGT', 'GAAGGCAGG'],
        })
    
    def test_exact_match(self, sample_df):
        """Exact accession match should return variant."""
        result = get_reference_variant_sequence(sample_df, 'CY021716')
        assert result == 'GAAGACAGG'
    
    def test_partial_match(self, sample_df):
        """Partial accession match should work."""
        result = get_reference_variant_sequence(sample_df, 'CY0207')
        # CY020292 contains 'CY0207'? No, but CY021716 does not either
        # Actually CY020292 does NOT contain CY0207
        # Let's test with a substring that works
        result = get_reference_variant_sequence(sample_df, 'CY02')
        # Should match first one containing 'CY02'
        assert result in ['GAAGACAGG', 'GAAGACAGT']
    
    def test_no_match_returns_none(self, sample_df):
        """No match should return None."""
        result = get_reference_variant_sequence(sample_df, 'NONEXISTENT')
        assert result is None


class TestNamingSchemeEnum:
    """Tests for NamingScheme enum."""
    
    def test_ranked_value(self):
        """RANKED should have value 'ranked'."""
        assert NamingScheme.RANKED.value == "ranked"
    
    def test_delta_value(self):
        """DELTA should have value 'delta'."""
        assert NamingScheme.DELTA.value == "delta"
    
    def test_from_string(self):
        """Should be able to create from string."""
        assert NamingScheme("ranked") == NamingScheme.RANKED
        assert NamingScheme("delta") == NamingScheme.DELTA
    
    def test_invalid_string_raises_error(self):
        """Invalid string should raise ValueError."""
        with pytest.raises(ValueError):
            NamingScheme("invalid")
