"""
Tests for visualization utilities.
"""

import pandas as pd
import pytest

from feavar.visualization import HAS_MATPLOTLIB, VariantPlotter


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason="matplotlib not installed")
def test_plot_stacked_by_field_with_delta_vts(tmp_path):
    """Ensure stacked plots work with delta-style VT labels."""
    merged_df = pd.DataFrame(
        {
            "accession": ["A", "B", "C", "D"],
            "variant_type": ["AAA", "AAA", "BBB", "CCC"],
            "VT": ["VT-REF", "VT-1A", "VT-1A", "VT-2G"],
            "host": ["human", "human", "avian", "human"],
        }
    )

    plotter = VariantPlotter(output_dir=tmp_path)
    plotter.plot_stacked_by_field(merged_df, field="host", top_n=2)

    assert plotter._current_ax is not None
    assert len(plotter._current_ax.patches) > 0

    plotter.close()
