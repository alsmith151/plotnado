"""Tests for grouping strategies."""

import pytest
from plotnado.cli.grouping import (
    SeqnadoSampleStrategy,
    SeqnadoAntibodyStrategy,
    HistoneMarkStrategy,
    RegexGroupingStrategy,
    detect_and_apply_grouping,
    PredefinedGroupingStrategies,
)


SEQNADO_FILES = [
    "sample1_H3K27ac.bw",
    "sample1_H3K4me3.bw",
    "sample2_H3K27ac.bw",
    "sample2_H3K4me3.bw",
]

HISTONE_FILES = [
    "ctrl_H3K27ac.bw",
    "treat_H3K27ac.bw",
    "ctrl_H3K4me3.bw",
]


class TestSeqnadoSampleStrategy:
    def test_groups_by_sample(self):
        result = SeqnadoSampleStrategy().apply(SEQNADO_FILES)
        assert result is not None
        assert "sample1_autoscale" in result.groups
        assert "sample2_autoscale" in result.groups
        assert set(result.groups["sample1_autoscale"]) == {0, 1}
        assert set(result.groups["sample2_autoscale"]) == {2, 3}

    def test_returns_none_for_non_seqnado(self):
        result = SeqnadoSampleStrategy().apply(["random.bw", "other.bw"])
        assert result is None

    def test_no_group_for_single_member(self):
        result = SeqnadoSampleStrategy().apply(["sample1_H3K27ac.bw"])
        assert result is None  # Single member, no group worth creating


class TestSeqnadoAntibodyStrategy:
    def test_groups_by_antibody(self):
        result = SeqnadoAntibodyStrategy().apply(SEQNADO_FILES)
        assert result is not None
        assert "H3K27ac_autoscale" in result.groups
        assert "H3K4me3_autoscale" in result.groups
        assert set(result.groups["H3K27ac_autoscale"]) == {0, 2}
        assert set(result.groups["H3K4me3_autoscale"]) == {1, 3}

    def test_returns_none_for_non_seqnado(self):
        result = SeqnadoAntibodyStrategy().apply(["random.bw"])
        assert result is None


class TestHistoneMarkStrategy:
    def test_groups_by_mark(self):
        result = HistoneMarkStrategy().apply(HISTONE_FILES)
        assert result is not None
        # Both H3K27ac files should be grouped
        h27ac_group = next(
            (k for k in result.groups if "h3k27ac" in k.lower()), None
        )
        assert h27ac_group is not None

    def test_returns_none_if_no_marks(self):
        result = HistoneMarkStrategy().apply(["sample_rep1.bw", "sample_rep2.bw"])
        assert result is None


class TestRegexGroupingStrategy:
    def test_groups_by_prefix(self):
        paths = ["control_rep1.bw", "control_rep2.bw", "treat_rep1.bw", "treat_rep2.bw"]
        strategy = RegexGroupingStrategy(r"([^_]+)_rep[0-9]+")
        result = strategy.apply(paths)
        assert result is not None
        assert any("control" in k for k in result.groups)
        assert any("treat" in k for k in result.groups)

    def test_no_match_returns_none(self):
        strategy = RegexGroupingStrategy(r"NOMATCH")
        result = strategy.apply(["a.bw", "b.bw"])
        assert result is None


class TestDetectAndApplyGrouping:
    def test_detects_seqnado(self):
        result = detect_and_apply_grouping(SEQNADO_FILES)
        assert result is not None
        assert result.strategy_name == "seqnado-sample"

    def test_detects_histone(self):
        result = detect_and_apply_grouping(HISTONE_FILES)
        assert result is not None

    def test_returns_none_for_unrelated_files(self):
        result = detect_and_apply_grouping(["alpha.bw", "beta.bw"])
        assert result is None


class TestPredefinedGroupingStrategies:
    def test_parse_sample(self):
        strategy = PredefinedGroupingStrategies.parse_group_by("sample")
        assert isinstance(strategy, SeqnadoSampleStrategy)

    def test_parse_antibody(self):
        strategy = PredefinedGroupingStrategies.parse_group_by("antibody")
        assert isinstance(strategy, SeqnadoAntibodyStrategy)

    def test_parse_regex(self):
        strategy = PredefinedGroupingStrategies.parse_group_by(r"([^_]+)_rep")
        assert isinstance(strategy, RegexGroupingStrategy)

    def test_invalid_regex_raises(self):
        with pytest.raises(ValueError):
            PredefinedGroupingStrategies.parse_group_by("[invalid")
