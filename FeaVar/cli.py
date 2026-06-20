"""
Command-line interface for FeaVar.

Provides the `feavar` command for running variant type analysis.
"""

import logging
import sys
from pathlib import Path

import click

from feavar import __version__
from feavar.core import FeaVarAnalysis
from feavar.exceptions import FeaVarError
from feavar.visualization import VariantPlotter


def setup_logging(level: str) -> None:
    """Configure logging based on verbosity level."""
    log_level = getattr(logging, level.upper(), logging.INFO)

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )


@click.command()
@click.version_option(__version__, prog_name="FeaVar")
@click.option(
    "-a",
    "--alignment",
    required=True,
    type=click.Path(exists=True),
    help="Path to the multiple sequence alignment file.",
)
@click.option(
    "-r",
    "--reference",
    required=True,
    help="Reference sequence identifier (e.g., CY021716).",
)
@click.option(
    "-p",
    "--positions",
    required=True,
    help='Positions to analyze. Examples: "100-110" or "100-110,120,130".',
)
@click.option(
    "-f",
    "--format",
    "alignment_format",
    default="auto",
    help='Alignment file format. Use "auto" to infer from extension (default: auto). '
    "Supported: clustal, fasta, phylip, nexus, stockholm.",
)
@click.option(
    "-m",
    "--metadata",
    type=click.Path(exists=True),
    help="Path to metadata file (tab-delimited) for merging.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    default=".",
    help="Output directory for results (default: current directory).",
)
@click.option(
    "-t",
    "--top",
    default=10,
    type=click.IntRange(min=1),
    help="Number of top variant types to plot (default: 10).",
)
@click.option(
    "-n",
    "--naming",
    "naming_scheme",
    default="delta",
    type=click.Choice(["delta", "ranked"], case_sensitive=False),
    help="Variant naming scheme (default: delta). "
    '"delta" uses stable names based on differences from reference (e.g., VT-5G.12C). '
    '"ranked" uses traditional frequency-based names (e.g., VT-001).',
)
@click.option(
    "--plot/--no-plot",
    default=True,
    help="Generate plots (default: True).",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (-v for INFO, -vv for DEBUG).",
)
def main(
    alignment: str,
    reference: str,
    positions: str,
    alignment_format: str,
    metadata: str | None,
    output: str,
    top: int,
    naming_scheme: str,
    plot: bool,
    verbose: int,
) -> None:
    """
    FeaVar - Sequence Feature Variant Type Analysis

    Compute variant types for specified positions in a multiple sequence alignment.

    \b
    Examples:
        feavar -a sequences.clw -r CY021716 -p "124-142"
        feavar -a sequences.clw -r CY021716 -p "100-110,120,130" -m metadata.tsv
        feavar -a sequences.fasta -r REF001 -p "50-100" -n ranked
    """
    # Set up logging
    log_level = ["WARNING", "INFO", "DEBUG"][min(verbose, 2)]
    setup_logging(log_level)

    logger = logging.getLogger(__name__)

    try:
        # Create output directory
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Run analysis
        click.echo("Running FeaVar analysis...")
        click.echo(f"  Alignment: {alignment}")
        click.echo(f"  Reference: {reference}")
        click.echo(f"  Positions: {positions}")
        click.echo(f"  Naming scheme: {naming_scheme}")

        analysis = FeaVarAnalysis(
            alignment_path=alignment,
            reference_id=reference,
            positions=positions,
            alignment_format=alignment_format,
            output_dir=output,
            naming_scheme=naming_scheme,
        )

        result = analysis.run()

        # Display summary
        click.echo("\nAnalysis complete!")
        click.echo(f"  Sequences analyzed: {result.num_sequences}")
        click.echo(f"  Unique variant types: {result.num_variant_types}")
        click.echo(f"\nTop {min(top, result.num_variant_types)} variant types:")

        top_variants = result.get_top_variants(top)
        for _, row in top_variants.iterrows():
            click.echo(f"  {row['VT']}: {row['variant_type']} (n={row['count']})")

        # Save results
        saved = analysis.save_results(prefix="feavar")
        click.echo("\nResults saved to:")
        for name, path in saved.items():
            click.echo(f"  {name}: {path}")

        # Merge metadata if provided
        if metadata:
            click.echo(f"\nMerging with metadata from: {metadata}")
            merged_df = analysis.merge_metadata(metadata)
            merged_path = output_dir / "feavar_merged.csv"
            merged_df.to_csv(merged_path, index=False)
            click.echo(f"  Saved merged data to: {merged_path}")

            # Plot by metadata fields if available
            if plot:
                plotter = VariantPlotter(output_dir=output_dir)

                # Find categorical columns for plotting
                categorical_cols = merged_df.select_dtypes(
                    include=["object", "category"]
                ).columns
                categorical_cols = [
                    c
                    for c in categorical_cols
                    if c not in ["accession", "variant_type", "VT"]
                ]

                for col in categorical_cols[:3]:  # Plot up to 3 fields
                    try:
                        plotter.plot_stacked_by_field(merged_df, col, top_n=top)
                        plot_path = plotter.save(f"feavar_by_{col}.png")
                        click.echo(f"  Saved plot: {plot_path}")
                        plotter.close()
                    except Exception as e:
                        logger.warning(f"Could not plot by {col}: {e}")

        # Generate standard plots
        if plot:
            try:
                plotter = VariantPlotter(output_dir=output_dir)

                # Bar chart of counts
                plotter.plot_variant_counts(result.summary_df, top_n=top)
                plot_path = plotter.save("feavar_counts.png")
                click.echo(f"\nPlot saved: {plot_path}")
                plotter.close()

                # Distribution chart
                plotter.plot_variant_distribution(result.summary_df, top_n=top)
                plot_path = plotter.save("feavar_distribution.png")
                click.echo(f"Plot saved: {plot_path}")
                plotter.close()

            except ImportError:
                click.echo("\nNote: matplotlib not installed, skipping plots")
            except Exception as e:
                logger.warning(f"Could not generate plots: {e}")

        click.echo("\nDone!")

    except FeaVarError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        logger.exception("Unexpected error")
        click.echo(f"Unexpected error: {e}", err=True)
        sys.exit(1)


def run_cli() -> None:
    """Entry point for the CLI."""
    main()


if __name__ == "__main__":
    run_cli()
