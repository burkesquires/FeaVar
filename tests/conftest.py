"""
Test configuration and fixtures for FeaVar.
"""

import pytest


@pytest.fixture
def sample_alignment_content():
    """Sample Clustal alignment content for testing."""
    # Note: Clustal format should NOT have sequence numbers at end of lines
    # The spacing must be consistent (use spaces, not tabs)
    return """CLUSTAL W (1.82) multiple sequence alignment


SEQ001          MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA
SEQ002          -----PGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAA---ECAGLGEMPGSFVPTVTA
SEQ003          MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA
                     *******************************************************

SEQ001          ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS
SEQ002          ITTSQDLQWLVQPTLISSMAQSQGQPLASQPP-----DMPGTSYSTPGMSGYSSGGASGS
SEQ003          ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS
                ********************************     ***********  **  ******
"""


@pytest.fixture
def sample_alignment_file(sample_alignment_content, tmp_path):
    """Create a temporary alignment file."""
    alignment_path = tmp_path / "test_alignment.clw"
    alignment_path.write_text(sample_alignment_content)
    return alignment_path


@pytest.fixture
def sample_metadata_content():
    """Sample metadata content for testing."""
    return """accession\thost\tyear\tlocation
SEQ001\thuman\t2020\tUSA
SEQ002\thuman\t2021\tUK
SEQ003\tavian\t2020\tChina
"""


@pytest.fixture
def sample_metadata_file(sample_metadata_content, tmp_path):
    """Create a temporary metadata file."""
    metadata_path = tmp_path / "test_metadata.tsv"
    metadata_path.write_text(sample_metadata_content)
    return metadata_path


@pytest.fixture
def complex_alignment_content():
    """More complex alignment with gaps for position adjustment testing."""
    # Proper Clustal format without sequence numbers
    return """CLUSTAL W (1.82) multiple sequence alignment


REFERENCE       --ABC--DEF---GHI----JKL
SEQ001          XXABCXXDEFXXXGHIXXXXJKL
SEQ002          YYABCYYDEFYYYGHIYYYYJKL
SEQ003          ZZABCZZDEFZZZGHIZZZZKL-
"""


@pytest.fixture
def complex_alignment_file(complex_alignment_content, tmp_path):
    """Create a temporary complex alignment file."""
    alignment_path = tmp_path / "complex_alignment.clw"
    alignment_path.write_text(complex_alignment_content)
    return alignment_path
