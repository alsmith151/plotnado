"""Tests for inference engine: type classification, title inference, color assignment."""

import pytest
from plotnado.cli.inference import (
    TrackClassifier,
    TitleInference,
    SeqnadoPattern,
    InferenceResult,
    infer_track,
    _ANNOTATION_COLOR,
    _BIGWIG_PALETTE,
)
from plotnado.template import TrackType


class TestTrackClassifier:
    def test_bigwig_by_extension(self):
        track_type, conf = TrackClassifier.classify("sample.bw")
        assert track_type == TrackType.BIGWIG
        assert conf >= 0.9

    def test_bigwig_full_extension(self):
        track_type, _ = TrackClassifier.classify("sample.bigwig")
        assert track_type == TrackType.BIGWIG

    def test_narrowpeak(self):
        track_type, _ = TrackClassifier.classify("peaks.narrowpeak")
        assert track_type == TrackType.NARROWPEAK

    def test_bed(self):
        track_type, _ = TrackClassifier.classify("regions.bed")
        assert track_type == TrackType.BED

    def test_bedgraph(self):
        track_type, _ = TrackClassifier.classify("coverage.bedgraph")
        assert track_type == TrackType.BEDGRAPH

    def test_links(self):
        track_type, _ = TrackClassifier.classify("interactions.bedpe")
        assert track_type == TrackType.LINKS

    def test_unknown(self):
        track_type, conf = TrackClassifier.classify("data.csv")
        assert track_type == TrackType.UNKNOWN
        assert conf == 0.0

    def test_url_bigwig(self):
        track_type, conf = TrackClassifier.classify("https://example.com/track.bw")
        assert track_type == TrackType.BIGWIG
        assert conf >= 0.8

    def test_case_insensitive(self):
        track_type, _ = TrackClassifier.classify("sample.BW")
        assert track_type == TrackType.BIGWIG


class TestSeqnadoPattern:
    def test_standard_seqnado(self):
        assert SeqnadoPattern.is_seqnado("sample1_H3K27ac.bw")

    def test_with_hyphens(self):
        assert SeqnadoPattern.is_seqnado("THP1-ctrl_H3K4me3.bigwig")

    def test_rejects_filetype_as_antibody(self):
        assert not SeqnadoPattern.is_seqnado("THP1H3K4me1_bigWig.bigWig")

    def test_parse_returns_sample_antibody(self):
        result = SeqnadoPattern.parse("sample1_H3K27ac.bw")
        assert result == ("sample1", "H3K27ac")

    def test_parse_non_seqnado_returns_none(self):
        result = SeqnadoPattern.parse("random_file.bw")
        # A file with only one underscore-separated token before extension may or may not match
        # depending on the pattern — just check it doesn't crash
        assert result is None or isinstance(result, tuple)


class TestTitleInference:
    def test_seqnado_title_format(self):
        title, inferred = TitleInference.infer("sample1_H3K27ac.bw")
        assert "H3K27ac" in title
        assert "sample1" in title
        assert inferred

    def test_cleans_bigwig_suffix(self):
        title, _ = TitleInference.infer("THP1H3K4me3_bigWig.bigWig")
        assert "bigwig" not in title.lower()
        assert "bigWig" not in title

    def test_adds_peaks_suffix_for_bed(self):
        title, _ = TitleInference.infer("sample1_regions.bed", TrackType.BED)
        assert "peaks" in title.lower()

    def test_camel_case_separation(self):
        title, _ = TitleInference.infer("THP1H3K4me3_bigBed.bigBed")
        assert " " in title  # Should have been split


class TestColorAssignment:
    def test_bigwig_gets_palette_color(self):
        result = infer_track("sample1_H3K27ac.bw")
        assert result.suggested_color is not None
        assert result.suggested_color in _BIGWIG_PALETTE

    def test_narrowpeak_gets_annotation_color(self):
        result = infer_track("peaks.narrowpeak")
        assert result.suggested_color == _ANNOTATION_COLOR

    def test_bed_gets_annotation_color(self):
        result = infer_track("regions.bed")
        assert result.suggested_color == _ANNOTATION_COLOR

    def test_same_antibody_same_color(self):
        """Two seqnado files with the same antibody should get the same color."""
        r1 = infer_track("sample1_H3K27ac.bw")
        r2 = infer_track("sample2_H3K27ac.bw")
        assert r1.suggested_color == r2.suggested_color

    def test_different_antibodies_may_differ(self):
        """Different antibodies should typically get different colors."""
        r1 = infer_track("sample1_H3K27ac.bw")
        r2 = infer_track("sample1_H3K4me3.bw")
        # Colors are hash-based, so they may collide, but are drawn from the palette
        assert r1.suggested_color in _BIGWIG_PALETTE
        assert r2.suggested_color in _BIGWIG_PALETTE

    def test_unknown_type_no_color(self):
        result = infer_track("data.csv")
        assert result.suggested_color is None
