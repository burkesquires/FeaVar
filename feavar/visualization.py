"""
Visualization utilities for FeaVar.

Provides plotting functions for variant type analysis results.
"""

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

# Try to import matplotlib - it's optional
try:
    import matplotlib

    # FeaVar only writes image files. A non-interactive backend makes that work
    # in CI, servers, notebooks, and macOS without an app context.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available - plotting disabled")


class VariantPlotter:
    """
    Plotter for variant type analysis results.

    Provides methods for creating bar charts and stacked bar plots
    of variant type distributions.

    Parameters
    ----------
    output_dir : str or Path, optional
        Directory for saving plots (default: current directory).
    figsize : Tuple[int, int], optional
        Default figure size (default: (12, 6)).
    dpi : int, optional
        Resolution for saved figures (default: 100).

    Examples
    --------
    >>> plotter = VariantPlotter(output_dir="plots")
    >>> plotter.plot_variant_counts(summary_df, top_n=10)
    >>> plotter.save("variant_counts.png")
    """

    def __init__(
        self,
        output_dir: str | Path | None = None,
        figsize: tuple[int, int] = (12, 6),
        dpi: int = 100,
    ):
        """Initialize the plotter."""
        if not HAS_MATPLOTLIB:
            logger.warning("Plotting functionality requires matplotlib")

        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        self.figsize = figsize
        self.dpi = dpi

        # Matplotlib is optional, so concrete plot types are unavailable when
        # plotting dependencies are omitted. Plot methods assign Figure/Axes.
        self._current_fig: Any | None = None
        self._current_ax: Any | None = None

    def _check_matplotlib(self) -> None:
        """Check if matplotlib is available."""
        if not HAS_MATPLOTLIB:
            raise ImportError(
                "matplotlib is required for plotting. "
                "Install with: pip install matplotlib"
            )

    def plot_variant_counts(
        self,
        summary_df: pd.DataFrame,
        top_n: int = 10,
        title: str | None = None,
        xlabel: str = "Variant Type",
        ylabel: str = "Count",
        color: str = "#3498db",
        show_values: bool = True,
    ) -> "VariantPlotter":
        """
        Create a bar chart of variant type counts.

        Parameters
        ----------
        summary_df : pd.DataFrame
            Summary DataFrame with 'VT' and 'count' columns.
        top_n : int, optional
            Number of top variants to show (default: 10).
        title : str, optional
            Plot title (default: auto-generated).
        xlabel : str, optional
            X-axis label.
        ylabel : str, optional
            Y-axis label.
        color : str, optional
            Bar color (default: blue).
        show_values : bool, optional
            Whether to show count values on bars (default: True).

        Returns
        -------
        VariantPlotter
            Self, for method chaining.
        """
        self._check_matplotlib()

        # Select top N variants
        plot_data = summary_df.head(top_n).copy()

        # Create figure
        self._current_fig, self._current_ax = plt.subplots(figsize=self.figsize)
        ax = self._current_ax

        # Create bar chart
        bars = ax.bar(
            plot_data["VT"],
            plot_data["count"],
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )

        # Add value labels on bars
        if show_values:
            for bar, count in zip(bars, plot_data["count"], strict=True):
                height = bar.get_height()
                ax.annotate(
                    f"{int(count)}",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

        # Set labels
        if title is None:
            title = f"Top {top_n} Variant Types by Count"
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)

        # Rotate x-axis labels
        plt.xticks(rotation=45, ha="right")

        # Add grid
        ax.yaxis.grid(True, linestyle="--", alpha=0.7)
        ax.set_axisbelow(True)

        plt.tight_layout()

        return self

    def plot_stacked_by_field(
        self,
        merged_df: pd.DataFrame,
        field: str,
        top_n: int = 10,
        title: str | None = None,
        colormap: str = "Set3",
    ) -> "VariantPlotter":
        """
        Create a stacked bar chart grouped by a metadata field.

        Parameters
        ----------
        merged_df : pd.DataFrame
            Merged DataFrame with 'VT' column and metadata fields.
        field : str
            Metadata field to group by.
        top_n : int, optional
            Number of top variants to show (default: 10).
        title : str, optional
            Plot title (default: auto-generated).
        colormap : str, optional
            Matplotlib colormap name (default: 'Set3').

        Returns
        -------
        VariantPlotter
            Self, for method chaining.
        """
        self._check_matplotlib()

        if field not in merged_df.columns:
            raise ValueError(
                f"Field '{field}' not found in DataFrame. "
                f"Available columns: {list(merged_df.columns)}"
            )

        if "VT" not in merged_df.columns:
            raise ValueError(
                "Column 'VT' not found in DataFrame. "
                f"Available columns: {list(merged_df.columns)}"
            )

        # Filter to top N variants
        top_vts = merged_df["VT"].value_counts().index[:top_n].tolist()

        if not top_vts:
            raise ValueError("No variant types available to plot.")

        plot_data = merged_df[merged_df["VT"].isin(top_vts)].copy()

        # Pivot for stacked bar
        pivot = pd.crosstab(plot_data["VT"], plot_data[field])

        # Sort by VT order
        pivot = pivot.reindex(top_vts).dropna(how="all")

        # Create figure
        self._current_fig, self._current_ax = plt.subplots(
            figsize=(self.figsize[0] * 1.5, self.figsize[1])
        )

        # Create stacked bar chart
        pivot.plot(
            kind="bar",
            stacked=True,
            ax=self._current_ax,
            colormap=colormap,
            edgecolor="white",
            linewidth=0.5,
        )

        # Set labels
        if title is None:
            title = f"Variant Types by {field.replace('_', ' ').title()}"
        self._current_ax.set_title(title, fontsize=14, fontweight="bold")
        self._current_ax.set_xlabel("Variant Type", fontsize=12)
        self._current_ax.set_ylabel("Count", fontsize=12)

        # Rotate x-axis labels
        plt.xticks(rotation=45, ha="right")

        # Adjust legend
        self._current_ax.legend(
            title=field.replace("_", " ").title(),
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )

        plt.tight_layout()

        return self

    def plot_variant_distribution(
        self,
        summary_df: pd.DataFrame,
        top_n: int = 20,
        title: str | None = None,
    ) -> "VariantPlotter":
        """
        Create a horizontal bar chart showing variant distribution.

        Parameters
        ----------
        summary_df : pd.DataFrame
            Summary DataFrame with 'VT', 'variant_type', and 'count' columns.
        top_n : int, optional
            Number of variants to show (default: 20).
        title : str, optional
            Plot title.

        Returns
        -------
        VariantPlotter
            Self, for method chaining.
        """
        self._check_matplotlib()

        plot_data = summary_df.head(top_n).copy()

        # Create figure
        fig_height = max(6, top_n * 0.4)
        self._current_fig, self._current_ax = plt.subplots(figsize=(10, fig_height))
        ax = self._current_ax

        # Create horizontal bar chart
        y_pos = range(len(plot_data))
        ax.barh(
            y_pos,
            plot_data["count"],
            color="#2ecc71",
            edgecolor="white",
        )

        # Set y-axis labels (VT and sequence)
        labels = [
            f"{row['VT']}: {row['variant_type'][:20]}..."
            if len(row["variant_type"]) > 20
            else f"{row['VT']}: {row['variant_type']}"
            for _, row in plot_data.iterrows()
        ]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontfamily="monospace", fontsize=9)

        # Invert y-axis so top variant is at top
        ax.invert_yaxis()

        # Set labels
        if title is None:
            title = f"Top {top_n} Variant Type Distribution"
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.set_xlabel("Count", fontsize=12)

        # Add grid
        ax.xaxis.grid(True, linestyle="--", alpha=0.7)
        ax.set_axisbelow(True)

        plt.tight_layout()

        return self

    def save(
        self,
        filename: str,
        format: str | None = None,
    ) -> Path:
        """
        Save the current plot to a file.

        Parameters
        ----------
        filename : str
            Output filename (can include path).
        format : str, optional
            File format (inferred from filename if not specified).

        Returns
        -------
        Path
            Path to the saved file.
        """
        self._check_matplotlib()

        if self._current_fig is None:
            raise RuntimeError("No plot to save. Create a plot first.")

        self.output_dir.mkdir(parents=True, exist_ok=True)

        filepath = self.output_dir / filename

        self._current_fig.savefig(
            filepath,
            format=format,
            dpi=self.dpi,
            bbox_inches="tight",
        )

        logger.info(f"Saved plot to {filepath}")

        return filepath

    def show(self) -> None:
        """Display the current plot."""
        self._check_matplotlib()

        if self._current_fig is None:
            raise RuntimeError("No plot to show. Create a plot first.")

        plt.show()

    def close(self) -> None:
        """Close the current figure."""
        if HAS_MATPLOTLIB and self._current_fig is not None:
            plt.close(self._current_fig)
            self._current_fig = None
            self._current_ax = None
