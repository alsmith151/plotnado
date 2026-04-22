"""Tests for IGV session XML parsing."""

import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from plotnado.igv import IgvSession, parse_igv_session, _igv_color_to_hex, _infer_track_type
from plotnado.template import Template, TrackSpec
from plotnado.tracks.enums import TrackType


IGV_XML = Path(__file__).parent / "data" / "igv_session_CAT_and_MPAL_v2.xml"


class TestColorConversion:
    def test_red(self):
        assert _igv_color_to_hex("255,0,0") == "#ff0000"

    def test_black(self):
        assert _igv_color_to_hex("0,0,0") == "#000000"

    def test_empty(self):
        assert _igv_color_to_hex("") is None

    def test_invalid(self):
        assert _igv_color_to_hex("not,a,color,string") is None


class TestTrackTypeInference:
    def test_bigwig(self):
        assert _infer_track_type("/path/to/signal.bw") == TrackType.BIGWIG

    def test_bed(self):
        assert _infer_track_type("peaks.bed") == TrackType.BED

    def test_narrowpeak(self):
        assert _infer_track_type("peaks.narrowPeak") == TrackType.NARROWPEAK

    def test_unknown(self):
        assert _infer_track_type("file.xyz") == TrackType.UNKNOWN


class TestParseIgvSession:
    def test_returns_igv_session(self):
        result = parse_igv_session(IGV_XML)
        assert isinstance(result, IgvSession)
        assert isinstance(result.template, Template)

    def test_genome_extracted(self):
        result = parse_igv_session(IGV_XML)
        assert result.genome == "hg38"
        assert result.template.genome == "hg38"

    def test_locus_extracted(self):
        result = parse_igv_session(IGV_XML)
        assert result.locus == "chr11:6211915-6260069"

    def test_gene_track_becomes_guide(self):
        result = parse_igv_session(IGV_XML)
        assert result.template.guides.genes is True

    def test_bigwig_tracks_parsed(self):
        result = parse_igv_session(IGV_XML)
        bw_tracks = [t for t in result.template.tracks if t.type == TrackType.BIGWIG]
        # 1 CAT track + 23 MPAL tracks = 24 total
        assert len(bw_tracks) == 24

    def test_cat_track_color(self):
        result = parse_igv_session(IGV_XML)
        cat_track = next(t for t in result.template.tracks if t.title == "CAT-MV411_H3K27ac")
        assert cat_track.color == "#ff0000"

    def test_mpal_tracks_black(self):
        result = parse_igv_session(IGV_XML)
        mpal_track = next(t for t in result.template.tracks if t.title == "25_NK")
        assert mpal_track.color == "#000000"

    def test_fixed_scale_options(self):
        result = parse_igv_session(IGV_XML)
        cat_track = next(t for t in result.template.tracks if t.title == "CAT-MV411_H3K27ac")
        assert cat_track.options["min_value"] == 0.0
        assert cat_track.options["max_value"] == pytest.approx(1769.38)

    def test_autoscale_group_assigned(self):
        result = parse_igv_session(IGV_XML)
        # MPAL tracks share autoscaleGroup="3"
        mpal_track = next(t for t in result.template.tracks if t.title == "25_NK")
        assert mpal_track.group == "igv_group_3"

    def test_group_specs_created(self):
        result = parse_igv_session(IGV_XML)
        group_names = {g.name for g in result.template.groups}
        assert "igv_group_3" in group_names

    def test_group_autoscale_true_autocolor_false(self):
        result = parse_igv_session(IGV_XML)
        mpal_group = next(g for g in result.template.groups if g.name == "igv_group_3")
        assert mpal_group.autoscale is True
        assert mpal_group.autocolor is False

    def test_axis_and_scalebar_guides(self):
        result = parse_igv_session(IGV_XML)
        assert result.template.guides.axis is True
        assert result.template.guides.scalebar is True


class TestFromIgvSession:
    def test_returns_figure_and_locus(self):
        from plotnado import GenomicFigure

        fig, locus = GenomicFigure.from_igv_session(IGV_XML)
        assert locus == "chr11:6211915-6260069"
        assert len(fig.tracks) > 0

    def test_figure_has_expected_track_count(self):
        from plotnado import GenomicFigure

        fig, _ = GenomicFigure.from_igv_session(IGV_XML)
        # 24 bigwig + scalebar + axis + genes = 27
        assert len(fig.tracks) == 27


class TestTrackEditing:
    @pytest.fixture
    def fig(self):
        from plotnado import GenomicFigure
        fig, _ = GenomicFigure.from_igv_session(IGV_XML)
        return fig

    def test_getitem_by_title(self, fig):
        track = fig["25_NK"]
        assert track.title == "25_NK"

    def test_getitem_case_insensitive(self, fig):
        assert fig["25_nk"] is fig["25_NK"]

    def test_getitem_by_index(self, fig):
        assert fig[0] is fig.tracks[0]

    def test_getitem_missing_title_raises(self, fig):
        with pytest.raises(KeyError):
            fig["nonexistent_track"]

    def test_getitem_out_of_range_raises(self, fig):
        with pytest.raises(IndexError):
            fig[9999]

    def test_direct_mutation_via_getitem(self, fig):
        fig["25_NK"].color = "steelblue"
        assert fig["25_NK"].color == "steelblue"

    def test_update_track_by_title(self, fig):
        fig.update_track("25_NK", color="coral", max_value=100.0)
        track = fig["25_NK"]
        assert track.color == "coral"
        assert track.max_value == 100.0

    def test_update_track_returns_self_for_chaining(self, fig):
        result = fig.update_track("25_NK", color="blue")
        assert result is fig

    def test_update_track_by_index(self, fig):
        original_title = fig[0].title
        fig.update_track(0, height=2.0)
        assert fig[0].height == 2.0
        assert fig[0].title == original_title

    def test_remove_track_by_title(self, fig):
        count = len(fig.tracks)
        fig.remove_track("25_NK")
        assert len(fig.tracks) == count - 1
        with pytest.raises(KeyError):
            fig["25_NK"]

    def test_remove_track_by_index(self, fig):
        first_title = fig[1].title
        fig.remove_track(0)
        assert fig[0].title == first_title

    def test_remove_track_returns_self(self, fig):
        assert fig.remove_track("25_NK") is fig

    def test_remove_missing_track_raises(self, fig):
        with pytest.raises(KeyError):
            fig.remove_track("nonexistent")

    def test_add_track_position_bottom(self, fig):
        from plotnado.tracks import BigWigTrack
        count = len(fig.tracks)
        fig.bigwig("/tmp/new.bw", title="new_track")
        assert len(fig.tracks) == count + 1
        assert fig.tracks[-1].title == "new_track"

    def test_add_track_position_top(self, fig):
        from plotnado.tracks import BigWigTrack
        fig.bigwig("/tmp/new.bw", title="new_track", position="top")
        assert fig.tracks[0].title == "new_track"


class TestUpdateTracks:
    @pytest.fixture
    def fig(self):
        from plotnado import GenomicFigure
        fig, _ = GenomicFigure.from_igv_session(IGV_XML)
        return fig

    def test_update_all_tracks(self, fig):
        fig.update_track(height=0.2)
        assert all(t.height == 0.2 for t in fig.tracks)

    def test_update_returns_self(self, fig):
        assert fig.update_track(height=0.5) is fig

    def test_filter_by_track_type_string(self, fig):
        from plotnado.tracks import BigWigTrack, ScaleBar
        non_bw_heights_before = {id(t): t.height for t in fig.tracks if not isinstance(t, BigWigTrack)}
        fig.update_track(track_type="bigwig", height=0.3)
        assert all(t.height == 0.3 for t in fig.tracks if isinstance(t, BigWigTrack))
        assert all(t.height == non_bw_heights_before[id(t)] for t in fig.tracks if not isinstance(t, BigWigTrack))

    def test_filter_by_track_type_class(self, fig):
        from plotnado.tracks import BigWigTrack
        fig.update_track(track_type=BigWigTrack, color="#aaaaaa")
        assert all(t.color == "#aaaaaa" for t in fig.tracks if isinstance(t, BigWigTrack))

    def test_filter_by_group(self, fig):
        from plotnado.tracks import BigWigTrack
        heights_before = {id(t): t.height for t in fig.tracks}
        fig.update_track(group="igv_group_3", height=0.5)
        for t in fig.tracks:
            if t.autoscale_group == "igv_group_3":
                assert t.height == 0.5
            else:
                assert t.height == heights_before[id(t)]

    def test_filter_by_where(self, fig):
        heights_before = {id(t): t.height for t in fig.tracks}
        fig.update_track(where=lambda t: t.title and t.title.startswith("0"), height=2.0)
        for t in fig.tracks:
            if t.title and t.title.startswith("0"):
                assert t.height == 2.0
            else:
                assert t.height == heights_before[id(t)]

    def test_combined_filters(self, fig):
        from plotnado.tracks import BigWigTrack
        heights_before = {id(t): t.height for t in fig.tracks}
        fig.update_track(track_type="bigwig", group="igv_group_3", height=0.4)
        for t in fig.tracks:
            if isinstance(t, BigWigTrack) and t.autoscale_group == "igv_group_3":
                assert t.height == 0.4
            else:
                assert t.height == heights_before[id(t)]
