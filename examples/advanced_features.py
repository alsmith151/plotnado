"""
Advanced features example: Scaling, Overlays, and Highlights.
"""

from plotnado import Figure
from plotnado.tracks import BigWigTrack, BigwigOverlay, GenomicRegion


def main():
    # Define a genomic region
    gr = GenomicRegion(chromosome="chr1", start=1000000, end=1100000)

    # 1. Automatic Scaling
    # Enable autoscale(True) share y-axis limits across all BigWig tracks
    fig_scale = Figure()
    fig_scale.autoscale(True)
    fig_scale.add_track("bigwig", data="signal1.bw", title="Track 1")
    fig_scale.add_track("bigwig", data="signal2.bw", title="Track 2")

    # 2. Overlays
    # Plot multiple bigwig signals on the same set of axes
    overlay = BigwigOverlay(
        tracks=[
            BigWigTrack(data="control.bw", title="Control"),
            BigWigTrack(data="treatment.bw", title="Treatment"),
        ],
        title="Signal Overlay",
    )

    fig_overlay = Figure()
    fig_overlay.add_track(overlay)

    # 3. Highlights and Region Extension
    # Extend the region by 10kb upstream and downstream
    extended_gr = gr.extend(upstream=10000, downstream=10000)

    fig_features = Figure()
    fig_features.highlight("chr1:1045000-1055000")  # Highlight a ROI
    fig_features.autocolor(palette="Set1")  # Automatically color tracks

    fig_features.add_track("axis")
    fig_features.add_track("bigwig", data="signal.bw")
    fig_features.add_track("genes", genome="hg38")

    print("Advanced setup examples created.")
    print(f"Original region: {gr}")
    print(f"Extended region: {extended_gr}")
    fig_features.save("advanced_features.png", region=extended_gr)


if __name__ == "__main__":
    main()
