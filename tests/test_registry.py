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
