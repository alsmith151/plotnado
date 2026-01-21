"""
Figure class for composing and plotting genomic tracks.
"""

from pathlib import Path
from typing import List, Optional, Union

import matplotlib.pyplot as plt
import matplotlib.figure
from loguru import logger

from .tracks import (
    BedTrack,
    BigWigTrack,
    Genes,
    GenomicAxis,
    GenomicRegion,
    HighlightsFromFile,
    ScaleBar,
    Spacer,
    BigwigOverlay,
    Track,
)


class Figure:
    """
    Compose and plot multiple genomic tracks.

    Example:
        >>> fig = Figure()
        >>> fig.add_track('scalebar')
        >>> fig.add_track('genes', genome='hg38')
        >>> fig.add_track(BigWigTrack(data='signal.bw', title='ChIP-seq'))
        >>> fig.plot('chr1:1000000-2000000')
    """

    def __init__(
        self,
        tracks: Optional[List[Track]] = None,
        width: float = 12,
        track_height: float = 2.0,
    ):
        """
        Initialize a Figure.

        Args:
            tracks: List of Track instances
            width: Figure width in inches
            track_height: Default height per track in inches
        """
        self.tracks: List[Track] = tracks or []
        self.width = width
        self.track_height = track_height
        self._highlight_regions: List[GenomicRegion] = []
        self._autocolor_palette: Optional[str] = None
        self._autoscale: bool = False

    def add_track(
        self,
        track: Union[str, Track],
        **kwargs,
    ) -> "Figure":
        """
        Add a track to the figure.

        Args:
            track: Either a Track instance or a string alias
            **kwargs: Additional arguments passed to track constructor

        Returns:
            Self for method chaining

        Example:
            >>> fig.add_track('scalebar')
            >>> fig.add_track('genes', genome='hg38')
            >>> fig.add_track(BigWigTrack(data='file.bw', title='Signal'))
        """
        if isinstance(track, str):
            track = self._create_track_from_alias(track, **kwargs)

        self.tracks.append(track)
        return self

    def autoscale(self, enable: bool = True) -> "Figure":
        """Enable or disable automatic scaling of y-axes."""
        self._autoscale = enable
        return self

    def autocolor(self, palette: str = "tab10") -> "Figure":
        """Enable automatic coloring of tracks using a palette."""
        self._autocolor_palette = palette
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        cmap = plt.get_cmap(palette)
        for i, track in enumerate(self.tracks):
            if hasattr(track, "aesthetics") and hasattr(track.aesthetics, "color"):
                # Use modulo to cycle through colors if more tracks than colors
                color = mcolors.to_hex(cmap(i % cmap.N))
                track.aesthetics.color = color
        return self

    def highlight(self, region: Union[str, GenomicRegion]) -> "Figure":
        """Highight a genomic region across all tracks."""
        if isinstance(region, str):
            region = GenomicRegion.from_str(region)
        self._highlight_regions.append(region)
        return self

    def _create_track_from_alias(self, alias: str, **kwargs) -> Track:
        """Create a track from a string alias."""
        alias_map = {
            "scalebar": ScaleBar,
            "scale": ScaleBar,
            "genes": Genes,
            "spacer": Spacer,
            "bigwig": BigWigTrack,
            "bed": BedTrack,
            "axis": GenomicAxis,
            "highlight": HighlightsFromFile,
            "bigwig_overlay": BigwigOverlay,
        }

        if alias.lower() not in alias_map:
            raise ValueError(
                f"Unknown track alias: {alias}. Available: {list(alias_map.keys())}"
            )

        track_class = alias_map[alias.lower()]
        return track_class(**kwargs)

    def plot(
        self,
        region: Union[str, GenomicRegion],
        show: bool = True,
        **kwargs,
    ) -> matplotlib.figure.Figure:
        """
        Plot all tracks for the given genomic region.

        Args:
            region: Genomic region as string ('chr1:1000-2000') or GenomicRegion
            show: Whether to display the plot
            **kwargs: Additional matplotlib figure kwargs

        Returns:
            matplotlib Figure object
        """
        if isinstance(region, str):
            gr = GenomicRegion.from_str(region)
        else:
            gr = region

        if not self.tracks:
            logger.warning("No tracks to plot")
            return None

        # Calculate heights
        heights = [track.height for track in self.tracks]
        total_height = sum(heights) * self.track_height

        # Create figure and axes
        fig, axes = plt.subplots(
            len(self.tracks),
            1,
            figsize=(self.width, total_height),
            gridspec_kw={"height_ratios": heights},
            **kwargs,
        )

        # Handle single track case
        if len(self.tracks) == 1:
            axes = [axes]

        # Apply autoscaling if enabled
        if self._autoscale:
            from .tracks.scaling import Autoscaler

            autoscaler = Autoscaler(tracks=self.tracks, gr=gr)
            autoscaler.apply()

        # Plot each track
        for ax, track in zip(axes, self.tracks):
            try:
                track.plot(ax, gr)

                # Draw highlights
                for highlight_gr in self._highlight_regions:
                    if highlight_gr.chromosome == gr.chromosome:
                        ax.axvspan(
                            highlight_gr.start,
                            highlight_gr.end,
                            color="yellow",
                            alpha=0.3,
                            zorder=-1,
                        )
            except Exception as e:
                logger.error(f"Error plotting track {track.__class__.__name__}: {e}")
                ax.text(
                    0.5,
                    0.5,
                    f"Error: {str(e)}",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                )

        plt.tight_layout()

        if show:
            plt.show()

        return fig

    def save(
        self,
        path: Union[str, Path],
        region: Union[str, GenomicRegion],
        dpi: int = 300,
        **kwargs,
    ) -> None:
        """
        Plot and save the figure to a file.

        Args:
            path: Output file path
            region: Genomic region to plot
            dpi: Resolution in dots per inch
            **kwargs: Additional arguments passed to savefig
        """
        fig = self.plot(region, show=False)
        if fig:
            fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
            logger.info(f"Saved figure to {path}")
            plt.close(fig)

    def __repr__(self) -> str:
        track_names = [t.__class__.__name__ for t in self.tracks]
        return f"Figure(tracks={track_names})"
