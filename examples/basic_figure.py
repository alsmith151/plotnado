"""
Basic example of Plotnado usage.
"""

from plotnado import Figure
from plotnado.tracks import GenomicRegion
import matplotlib.pyplot as plt


def main():
    # Define a genomic region
    gr = GenomicRegion(chromosome="chr1", start=1000000, end=1100000)

    # Create a figure
    fig = Figure(width=12, track_height=1.5)

    # Add tracks using aliases
    fig.add_track("scalebar")
    fig.add_track("axis")

    # Add a genes track (using mock data or path if available)
    # For demonstration, we'll just show the API
    fig.add_track("genes", title="RefSeq Genes", genome="hg38")

    # Add a spacer
    fig.add_track("spacer", height=0.5)

    # Add a BigWig track
    # fig.add_track("bigwig", data="path/to/signal.bw", title="ChIP-seq")

    # Plot the region
    # Note: This will try to fetch data, so it might fail if files don't exist
    # Here we just show the setup
    print(f"Figure setup complete with {len(fig.tracks)} tracks.")
    print(fig)
    fig.save("basic_figure.png", region=gr)


if __name__ == "__main__":
    main()
