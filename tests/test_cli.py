"""
Tests for the command-line interface.
"""

import pytest
from click.testing import CliRunner

from feavar import __version__
from feavar.cli import main


class TestCLI:
    """Tests for the CLI."""

    @pytest.fixture
    def runner(self):
        """Create a CLI runner."""
        return CliRunner()

    def test_version_flag(self, runner):
        """Test --version flag."""
        result = runner.invoke(main, ["--version"])

        assert result.exit_code == 0
        assert __version__ in result.output

    def test_help_flag(self, runner):
        """Test --help flag."""
        result = runner.invoke(main, ["--help"])

        assert result.exit_code == 0
        assert "FeaVar" in result.output
        assert "--alignment" in result.output
        assert "--reference" in result.output
        assert "--positions" in result.output

    def test_missing_required_args(self, runner):
        """Test that missing required args shows error."""
        result = runner.invoke(main, [])

        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    def test_basic_analysis(self, runner, sample_alignment_file, tmp_path):
        """Test basic analysis run."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
                "-o",
                str(tmp_path),
                "--no-plot",  # Skip plotting to avoid matplotlib dependency
            ],
        )

        assert result.exit_code == 0
        assert "Analysis complete" in result.output
        assert "Sequences analyzed" in result.output
        assert "variant types" in result.output.lower()

    def test_output_files_created(self, runner, sample_alignment_file, tmp_path):
        """Test that output files are created."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
                "-o",
                str(tmp_path),
                "--no-plot",
            ],
        )

        assert result.exit_code == 0

        # Check output files exist
        variants_file = tmp_path / "feavar_variants.csv"
        summary_file = tmp_path / "feavar_summary.csv"

        assert variants_file.exists()
        assert summary_file.exists()

    def test_with_metadata(
        self, runner, sample_alignment_file, sample_metadata_file, tmp_path
    ):
        """Test analysis with metadata."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
                "-m",
                str(sample_metadata_file),
                "-o",
                str(tmp_path),
                "--no-plot",
            ],
        )

        assert result.exit_code == 0
        assert "Merging with metadata" in result.output

        merged_file = tmp_path / "feavar_merged.csv"
        assert merged_file.exists()

    def test_verbose_flag(self, runner, sample_alignment_file, tmp_path):
        """Test verbose output."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
                "-o",
                str(tmp_path),
                "-v",
                "--no-plot",
            ],
        )

        assert result.exit_code == 0

    def test_invalid_alignment_file(self, runner, tmp_path):
        """Test with invalid alignment file."""
        fake_file = tmp_path / "fake.clw"

        result = runner.invoke(
            main,
            [
                "-a",
                str(fake_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
            ],
        )

        assert result.exit_code != 0

    def test_invalid_reference(self, runner, sample_alignment_file, tmp_path):
        """Test with invalid reference sequence."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "NONEXISTENT",
                "-p",
                "10-20",
                "-o",
                str(tmp_path),
                "--no-plot",
            ],
        )

        assert result.exit_code != 0
        assert "not found" in result.output.lower() or "error" in result.output.lower()

    def test_top_n_option(self, runner, sample_alignment_file, tmp_path):
        """Test --top option."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
                "-o",
                str(tmp_path),
                "-t",
                "5",
                "--no-plot",
            ],
        )

        assert result.exit_code == 0
        assert "Top 5" in result.output or "Top" in result.output

    def test_custom_format(self, runner, sample_alignment_file, tmp_path):
        """Test custom alignment format option."""
        result = runner.invoke(
            main,
            [
                "-a",
                str(sample_alignment_file),
                "-r",
                "SEQ001",
                "-p",
                "10-20",
                "-f",
                "clustal",
                "-o",
                str(tmp_path),
                "--no-plot",
            ],
        )

        assert result.exit_code == 0
