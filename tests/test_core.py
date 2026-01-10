"""
Tests for the core analysis module.
"""

import pytest
import pandas as pd
from pathlib import Path

from feavar.core import FeaVarAnalysis, AnalysisResult
from feavar.exceptions import FeaVarError, MetadataError


class TestAnalysisResult:
    """Tests for AnalysisResult dataclass."""
    
    def test_analysis_result_creation(self):
        """Test creating an AnalysisResult."""
        variants_df = pd.DataFrame({
            'accession': ['SEQ001', 'SEQ002', 'SEQ003'],
            'variant_type': ['ABC', 'ABC', 'DEF'],
        })
        summary_df = pd.DataFrame({
            'variant_type': ['ABC', 'DEF'],
            'count': [2, 1],
            'VT': ['VT-001', 'VT-002'],
        })
        
        result = AnalysisResult(
            variants_df=variants_df,
            summary_df=summary_df,
            reference_id='SEQ001',
            positions=[1, 2, 3],
            adjusted_positions=[1, 2, 3],
        )
        
        assert result.num_sequences == 3
        assert result.num_variant_types == 2
    
    def test_get_top_variants(self):
        """Test getting top variants."""
        summary_df = pd.DataFrame({
            'variant_type': ['A', 'B', 'C', 'D'],
            'count': [100, 50, 25, 10],
            'VT': ['VT-001', 'VT-002', 'VT-003', 'VT-004'],
        })
        
        result = AnalysisResult(
            variants_df=pd.DataFrame({'accession': [], 'variant_type': []}),
            summary_df=summary_df,
            reference_id='REF',
            positions=[1],
            adjusted_positions=[1],
        )
        
        top2 = result.get_top_variants(2)
        assert len(top2) == 2
        assert top2.iloc[0]['VT'] == 'VT-001'
    
    def test_get_variant_for_accession(self):
        """Test getting variant for specific accession."""
        variants_df = pd.DataFrame({
            'accession': ['SEQ001', 'SEQ002', 'SEQ003'],
            'variant_type': ['ABC', 'DEF', 'GHI'],
        })
        
        result = AnalysisResult(
            variants_df=variants_df,
            summary_df=pd.DataFrame(),
            reference_id='REF',
            positions=[1],
            adjusted_positions=[1],
        )
        
        assert result.get_variant_for_accession('SEQ002') == 'DEF'
        assert result.get_variant_for_accession('NONEXISTENT') is None


class TestFeaVarAnalysis:
    """Tests for FeaVarAnalysis class."""
    
    def test_init(self, sample_alignment_file):
        """Test initialization."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        assert analysis.alignment_path == str(sample_alignment_file)
        assert analysis.reference_id == 'SEQ001'
        assert analysis.positions_str == '10-20'
    
    def test_validate(self, sample_alignment_file):
        """Test validation step."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        assert analysis.validate() is True
    
    def test_validate_invalid_reference(self, sample_alignment_file):
        """Test validation with invalid reference."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='NONEXISTENT',
            positions='10-20',
        )
        
        with pytest.raises(FeaVarError):
            analysis.validate()
    
    def test_run_analysis(self, sample_alignment_file):
        """Test running full analysis."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        result = analysis.run()
        
        assert isinstance(result, AnalysisResult)
        assert result.num_sequences == 3
        assert result.num_variant_types > 0
        assert len(result.variants_df) == 3
        assert 'accession' in result.variants_df.columns
        assert 'variant_type' in result.variants_df.columns
    
    def test_run_analysis_summary_format(self, sample_alignment_file):
        """Test that summary has correct format."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        result = analysis.run()
        
        assert 'variant_type' in result.summary_df.columns
        assert 'count' in result.summary_df.columns
        assert 'VT' in result.summary_df.columns
        
        # VT labels should be in format VT-001, VT-002, etc.
        vt_labels = result.summary_df['VT'].tolist()
        assert vt_labels[0] == 'VT-001'
    
    def test_run_analysis_sorted_by_count(self, sample_alignment_file):
        """Test that summary is sorted by count descending."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        result = analysis.run()
        
        counts = result.summary_df['count'].tolist()
        assert counts == sorted(counts, reverse=True)
    
    def test_get_result_before_run(self, sample_alignment_file):
        """Test that get_result returns None before run."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        assert analysis.get_result() is None
    
    def test_get_result_after_run(self, sample_alignment_file):
        """Test that get_result returns result after run."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        analysis.run()
        result = analysis.get_result()
        
        assert result is not None
        assert isinstance(result, AnalysisResult)


class TestFeaVarAnalysisMetadata:
    """Tests for metadata merging functionality."""
    
    def test_merge_metadata(self, sample_alignment_file, sample_metadata_file):
        """Test merging with metadata."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        analysis.run()
        merged = analysis.merge_metadata(str(sample_metadata_file))
        
        assert isinstance(merged, pd.DataFrame)
        assert 'host' in merged.columns
        assert 'year' in merged.columns
        assert 'VT' in merged.columns
    
    def test_merge_metadata_before_run_raises(
        self, sample_alignment_file, sample_metadata_file
    ):
        """Test that merging before run raises error."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        with pytest.raises(FeaVarError, match="Must run analysis"):
            analysis.merge_metadata(str(sample_metadata_file))
    
    def test_merge_metadata_nonexistent_file_raises(self, sample_alignment_file):
        """Test that nonexistent metadata file raises error."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
        )
        
        analysis.run()
        
        with pytest.raises(MetadataError, match="Failed to load"):
            analysis.merge_metadata("nonexistent.tsv")


class TestFeaVarAnalysisSaveResults:
    """Tests for saving results functionality."""
    
    def test_save_results(self, sample_alignment_file, tmp_path):
        """Test saving results to files."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
            output_dir=str(tmp_path),
        )
        
        analysis.run()
        saved = analysis.save_results(prefix="test")
        
        assert 'variants' in saved
        assert 'summary' in saved
        assert saved['variants'].exists()
        assert saved['summary'].exists()
    
    def test_save_results_content(self, sample_alignment_file, tmp_path):
        """Test content of saved files."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
            output_dir=str(tmp_path),
        )
        
        analysis.run()
        saved = analysis.save_results()
        
        # Load and verify variants file
        variants_df = pd.read_csv(saved['variants'])
        assert 'accession' in variants_df.columns
        assert 'variant_type' in variants_df.columns
        
        # Load and verify summary file
        summary_df = pd.read_csv(saved['summary'])
        assert 'VT' in summary_df.columns
        assert 'count' in summary_df.columns
    
    def test_save_results_before_run_raises(self, sample_alignment_file, tmp_path):
        """Test that saving before run raises error."""
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
            output_dir=str(tmp_path),
        )
        
        with pytest.raises(FeaVarError, match="Must run analysis"):
            analysis.save_results()
    
    def test_save_results_creates_directory(self, sample_alignment_file, tmp_path):
        """Test that save creates output directory if needed."""
        output_dir = tmp_path / "new_dir" / "nested"
        
        analysis = FeaVarAnalysis(
            alignment_path=str(sample_alignment_file),
            reference_id='SEQ001',
            positions='10-20',
            output_dir=str(output_dir),
        )
        
        analysis.run()
        analysis.save_results()
        
        assert output_dir.exists()
