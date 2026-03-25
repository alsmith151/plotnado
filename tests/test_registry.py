"""Tests for the TrackRegistry and unified TrackType enum."""

import pytest
from plotnado.tracks.enums import TrackType


def test_tracktype_covers_all_former_templatetracktype_values():
    """All values that were in TemplateTrackType must exist in TrackType."""
    former_template_types = [
        "bigwig", "bed", "narrowpeak", "bedgraph",
        "gene", "links", "annotation", "overlay", "unknown",
    ]
    for v in former_template_types:
        assert TrackType(v), f"Missing TrackType value: {v}"


def test_tracktype_covers_all_former_alias_map_keys():
    """All keys from GenomicFigure._alias_map() must exist in TrackType."""
    former_alias_keys = [
        "scalebar", "scale", "genes", "spacer", "bigwig", "bed",
        "axis", "highlight", "bigwig_overlay", "overlay", "narrowpeak",
        "links", "hline", "vline", "cooler", "capcruncher",
        "cooler_average", "bigwig_collection", "bigwig_diff",
        "quantnado_coverage", "quantnado_stranded_coverage",
        "quantnado_methylation", "quantnado_variant",
    ]
    all_values = {t.value for t in TrackType}
    missing = [k for k in former_alias_keys if k not in all_values]
    assert not missing, f"TrackType missing values: {missing}"


def test_base_aesthetics_exists():
    """BaseAesthetics and BaseMultiColorAesthetics must be defined."""
    from pydantic import BaseModel
    from plotnado.tracks.aesthetics import BaseAesthetics, BaseMultiColorAesthetics

    assert issubclass(BaseAesthetics, BaseModel)
    assert issubclass(BaseMultiColorAesthetics, BaseModel)
    # Must have the shared fields
    assert "color" in BaseAesthetics.model_fields
    assert "alpha" in BaseAesthetics.model_fields
    assert "linewidth" in BaseAesthetics.model_fields
    assert "colors" in BaseMultiColorAesthetics.model_fields
    assert "color" not in BaseMultiColorAesthetics.model_fields  # no single color


def test_registry_module_imports():
    """Registry module must be importable and expose the singleton."""
    from plotnado.tracks.registry import TrackRegistry, TrackEntry, registry
    assert isinstance(registry, TrackRegistry)


def test_registry_register_decorator():
    """@registry.register must store the class under the given TrackType."""
    from plotnado.tracks.registry import TrackRegistry
    from plotnado.tracks.enums import TrackType
    from plotnado.tracks.base import Track

    local_registry = TrackRegistry()

    @local_registry.register(TrackType.BIGWIG, aliases=["bw"])
    class FakeTrack(Track):
        def fetch_data(self, gr): ...
        def plot(self, ax, gr): ...

    entry = local_registry.get("bigwig")
    assert entry.cls is FakeTrack
    assert entry.track_type is TrackType.BIGWIG

    alias_entry = local_registry.get("bw")
    assert alias_entry.cls is FakeTrack


def test_registry_unknown_key_raises():
    from plotnado.tracks.registry import TrackRegistry
    r = TrackRegistry()
    with pytest.raises(KeyError, match="Unknown track type"):
        r.get("nonexistent_type")


def test_signal_tracks_registered():
    """Signal track classes are registered after tracks/__init__.py is imported."""
    import plotnado.tracks  # triggers all @registry.register decorators
    from plotnado.tracks.registry import registry
    from plotnado.tracks.bigwig import BigWigTrack
    from plotnado.tracks.bigwig_diff import BigWigDiff
    from plotnado.tracks.bigwig_collection import BigWigCollection
    from plotnado.tracks.overlay import OverlayTrack

    assert registry.get("bigwig").cls is BigWigTrack
    assert registry.get("bedgraph").cls is BigWigTrack   # alias
    assert registry.get("bigwig_diff").cls is BigWigDiff
    assert registry.get("bigwig_collection").cls is BigWigCollection
    assert registry.get("overlay").cls is OverlayTrack
    assert registry.get("bigwig_overlay").cls is OverlayTrack  # alias


def test_bigwig_aesthetics_inherits_base():
    from plotnado.tracks.bigwig import BigwigAesthetics
    from plotnado.tracks.aesthetics import BaseAesthetics
    assert issubclass(BigwigAesthetics, BaseAesthetics)
    # Base fields are present
    assert "color" in BigwigAesthetics.model_fields
    assert "alpha" in BigwigAesthetics.model_fields
    # Track-specific field still present
    assert "style" in BigwigAesthetics.model_fields


def test_overlay_aesthetics_inherits_multi_color():
    from plotnado.tracks.overlay import OverlayTrackAesthetics
    from plotnado.tracks.aesthetics import BaseMultiColorAesthetics
    assert issubclass(OverlayTrackAesthetics, BaseMultiColorAesthetics)
    assert "colors" in OverlayTrackAesthetics.model_fields
    assert "color" not in OverlayTrackAesthetics.model_fields


def test_interval_and_guide_tracks_registered():
    """Interval and guide track classes are registered after tracks/__init__.py is imported."""
    import plotnado.tracks  # triggers all @registry.register decorators
    from plotnado.tracks.registry import registry
    from plotnado.tracks.bed import BedTrack
    from plotnado.tracks.peaks import NarrowPeakTrack
    from plotnado.tracks.genes import Genes
    from plotnado.tracks.links import LinksTrack
    from plotnado.tracks.axis import GenomicAxis
    from plotnado.tracks.scalebar import ScaleBar
    from plotnado.tracks.spacer import Spacer
    from plotnado.tracks.highlight import HighlightsFromFile
    from plotnado.tracks.annotations import HLineTrack, VLineTrack

    assert registry.get("bed").cls is BedTrack
    assert registry.get("annotation").cls is BedTrack      # alias
    assert registry.get("unknown").cls is BedTrack          # fallback alias
    assert registry.get("narrowpeak").cls is NarrowPeakTrack
    assert registry.get("gene").cls is Genes
    assert registry.get("genes").cls is Genes               # alias
    assert registry.get("links").cls is LinksTrack
    assert registry.get("axis").cls is GenomicAxis
    assert registry.get("scalebar").cls is ScaleBar
    assert registry.get("scale").cls is ScaleBar            # alias
    assert registry.get("spacer").cls is Spacer
    assert registry.get("highlight").cls is HighlightsFromFile
    assert registry.get("hline").cls is HLineTrack
    assert registry.get("vline").cls is VLineTrack


def test_bed_aesthetics_inherits_base():
    from plotnado.tracks.bed import BedAesthetics
    from plotnado.tracks.aesthetics import BaseAesthetics
    assert issubclass(BedAesthetics, BaseAesthetics)


def test_cooler_and_quantnado_tracks_registered():
    import plotnado.tracks
    from plotnado.tracks.registry import registry

    # Hi-C tracks
    from plotnado.tracks.cooler_track import CoolerTrack, CapcruncherTrack, CoolerAverage
    assert registry.get("cooler").cls is CoolerTrack
    assert registry.get("capcruncher").cls is CapcruncherTrack
    assert registry.get("cooler_average").cls is CoolerAverage

    # QuantNado tracks (optional dependency — skip if not installed)
    try:
        from plotnado.tracks.quantnado import (
            QuantNadoCoverageTrack, QuantNadoStrandedCoverageTrack,
            QuantNadoMethylationTrack, QuantNadoVariantTrack,
        )
        assert registry.get("quantnado_coverage").cls is QuantNadoCoverageTrack
        assert registry.get("quantnado_stranded_coverage").cls is QuantNadoStrandedCoverageTrack
        assert registry.get("quantnado_methylation").cls is QuantNadoMethylationTrack
        assert registry.get("quantnado_variant").cls is QuantNadoVariantTrack
    except ImportError:
        pytest.skip("quantnado not installed")


def test_cooler_aesthetics_stays_base_model():
    """CoolerAesthetics must NOT inherit BaseAesthetics (no single color field)."""
    from plotnado.tracks.cooler_track import CoolerAesthetics
    from plotnado.tracks.aesthetics import BaseAesthetics
    from pydantic import BaseModel
    assert issubclass(CoolerAesthetics, BaseModel)
    assert not issubclass(CoolerAesthetics, BaseAesthetics)


def test_registry_fully_populated_after_tracks_import():
    """Importing plotnado.tracks must populate the registry with all canonical types."""
    import plotnado.tracks  # noqa: F401 — side-effect import
    from plotnado.tracks.registry import registry

    # All canonical (non-alias) types must be registered
    canonical_types = [
        "bigwig", "bigwig_diff", "bigwig_collection", "overlay",
        "bed", "narrowpeak", "gene", "links",
        "scalebar", "axis", "spacer", "highlight", "hline", "vline",
        "cooler",
    ]
    for t in canonical_types:
        entry = registry.get(t)
        assert entry is not None, f"Not registered: {t}"
