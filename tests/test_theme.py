"""Tests for Theme.apply() method."""

import pytest


def test_theme_apply_sets_unset_fields():
    """Theme.apply() only modifies aesthetics fields not explicitly set."""
    from plotnado.theme import Theme
    from plotnado.tracks.bigwig import BigWigTrack

    theme = Theme.publication()
    track = BigWigTrack(data="fake.bw")  # no explicit aesthetics
    theme.apply(track, palette_color="#123456")

    assert track.aesthetics.color == "#123456"  # palette_color applied


def test_theme_apply_respects_explicit_color():
    """Theme.apply() must not overwrite an explicitly set color."""
    from plotnado.theme import Theme
    from plotnado.tracks.bigwig import BigWigTrack

    theme = Theme.publication()
    track = BigWigTrack(data="fake.bw", aesthetics={"color": "#aabbcc"})
    theme.apply(track, palette_color="#999999")

    assert track.aesthetics.color == "#aabbcc"  # explicit value preserved


def test_theme_apply_sets_alpha_when_unset():
    from plotnado.theme import Theme
    from plotnado.tracks.bigwig import BigWigTrack

    theme = Theme(alpha=0.5)
    track = BigWigTrack(data="fake.bw")
    theme.apply(track)

    assert track.aesthetics.alpha == 0.5
