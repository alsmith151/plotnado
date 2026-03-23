"""Tests for Template model: round-trip, optional groups, name field, annotated YAML."""

import pytest
import yaml
from pathlib import Path

from plotnado.template import Template, TrackSpec, GroupSpec, TrackType


def make_template() -> Template:
    t = Template()
    t.genome = "hg38"
    t.tracks = [
        TrackSpec(path="sample1.bw", type=TrackType.BIGWIG, title="H3K27ac (Sample 1)", group="s1"),
        TrackSpec(path="sample2.bw", type=TrackType.BIGWIG, title="H3K27ac (Sample 2)", group="s2"),
        TrackSpec(path="peaks.narrowpeak", type=TrackType.NARROWPEAK, title="Peaks"),
    ]
    return t


class TestRoundTrip:
    def test_to_yaml_and_load(self, tmp_path):
        t = make_template()
        yaml_path = tmp_path / "template.yaml"
        t.save(yaml_path)

        loaded = Template.load(yaml_path)
        assert loaded.genome == "hg38"
        assert len(loaded.tracks) == 3
        assert loaded.tracks[0].title == "H3K27ac (Sample 1)"
        assert loaded.tracks[0].group == "s1"

    def test_name_field_preserved(self, tmp_path):
        t = Template()
        t.tracks = [
            TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="Track A", name="track-a"),
        ]
        yaml_path = tmp_path / "t.yaml"
        t.save(yaml_path)

        loaded = Template.load(yaml_path)
        assert loaded.tracks[0].name == "track-a"

    def test_guides_round_trip(self, tmp_path):
        t = Template()
        t.guides.genes = True
        t.guides.axis = False
        yaml_path = tmp_path / "t.yaml"
        t.save(yaml_path)

        loaded = Template.load(yaml_path)
        assert loaded.guides.genes is True
        assert loaded.guides.axis is False


class TestOptionalGroupsSection:
    def test_groups_omitted_when_all_defaults(self, tmp_path):
        """Groups with default autoscale/autocolor should not appear in YAML."""
        t = make_template()
        # All-default groups — should be omitted
        t.groups = [GroupSpec(name="s1", tracks=["H3K27ac (Sample 1)"], autoscale=True, autocolor=True)]
        yaml_str = t.to_yaml()
        data = yaml.safe_load(yaml_str)
        assert "groups" not in data

    def test_groups_included_when_non_default(self, tmp_path):
        """Groups with non-default flags should appear in YAML."""
        t = make_template()
        t.groups = [GroupSpec(name="s1", tracks=["H3K27ac (Sample 1)"], autoscale=True, autocolor=False)]
        yaml_str = t.to_yaml()
        data = yaml.safe_load(yaml_str)
        assert "groups" in data
        assert data["groups"][0]["autocolor"] is False

    def test_template_with_no_groups_section_loads(self, tmp_path):
        """A template YAML with no groups key should load without error."""
        yaml_content = """
genome: hg38
tracks:
  - path: a.bw
    type: bigwig
    title: Track A
    group: my_group
"""
        yaml_path = tmp_path / "minimal.yaml"
        yaml_path.write_text(yaml_content)
        t = Template.load(yaml_path)
        assert t.tracks[0].group == "my_group"
        assert t.groups == []


class TestAnnotatedYAML:
    def test_header_added_when_header_args_given(self):
        t = make_template()
        yaml_str = t.to_yaml(header_args="sample1.bw sample2.bw")
        assert yaml_str.startswith("# PlotNado template")
        assert "plotnado init sample1.bw sample2.bw" in yaml_str

    def test_no_header_by_default(self):
        t = make_template()
        yaml_str = t.to_yaml()
        assert not yaml_str.startswith("#")
