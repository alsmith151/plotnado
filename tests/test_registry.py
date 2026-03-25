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
