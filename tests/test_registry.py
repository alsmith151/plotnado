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
