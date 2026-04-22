"""Tests for TemplateCompiler: no-mutation regression, autocolor, group defaults."""

import pytest
from plotnado.template import Template, TrackSpec, GroupSpec
from plotnado.tracks.enums import TrackType
from plotnado.render import TemplateCompiler


def make_template_with_groups(autocolor: bool = True) -> Template:
    t = Template()
    t.tracks = [
        TrackSpec(path="s1.bw", type=TrackType.BIGWIG, title="H3K27ac (Sample 1)", group="s1"),
        TrackSpec(path="s2.bw", type=TrackType.BIGWIG, title="H3K27ac (Sample 2)", group="s2"),
        TrackSpec(path="peaks.narrowpeak", type=TrackType.NARROWPEAK, title="Peaks (Sample 1)", group="s1"),
    ]
    t.groups = [
        GroupSpec(name="s1", tracks=["H3K27ac (Sample 1)", "Peaks (Sample 1)"], autoscale=True, autocolor=autocolor),
        GroupSpec(name="s2", tracks=["H3K27ac (Sample 2)"], autoscale=True, autocolor=autocolor),
    ]
    return t


class TestNoMutationRegression:
    def test_compile_twice_preserves_group_tracks(self):
        """Compiling a template must not overwrite group.tracks with indices."""
        t = make_template_with_groups()
        original_tracks = [list(g.tracks) for g in t.groups]

        TemplateCompiler.compile(t)
        after_first = [list(g.tracks) for g in t.groups]

        TemplateCompiler.compile(t)
        after_second = [list(g.tracks) for g in t.groups]

        assert after_first == original_tracks, "First compile mutated group.tracks"
        assert after_second == original_tracks, "Second compile mutated group.tracks"

    def test_resolved_indices_stored_in_plan(self):
        """Resolved numeric indices should be in RenderPlan, not the Template."""
        t = make_template_with_groups()
        plan = TemplateCompiler.compile(t)

        assert "s1" in plan.resolved_group_indices
        assert plan.resolved_group_indices["s1"] == [0, 2]


class TestAutocolorFlag:
    def test_autocolor_true_sets_color_group(self):
        t = make_template_with_groups(autocolor=True)
        plan = TemplateCompiler.compile(t)
        resolved = plan.tracks[0]  # s1 track
        assert resolved.color_group == "s1"

    def test_autocolor_false_clears_color_group(self):
        t = make_template_with_groups(autocolor=False)
        plan = TemplateCompiler.compile(t)
        resolved = plan.tracks[0]  # s1 track
        assert resolved.color_group is None

    def test_figure_kwargs_include_autoscale_but_not_color_when_autocolor_false(self):
        t = make_template_with_groups(autocolor=False)
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert "autoscale_group" in kwargs
        assert "color_group" not in kwargs


class TestGroupDefaults:
    def test_track_with_group_but_no_groupspec_still_gets_autoscale(self):
        """Track.group with no matching GroupSpec should still set autoscale_group."""
        t = Template()
        t.tracks = [
            TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="Track A", group="my_group"),
        ]
        # No groups section
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert kwargs.get("autoscale_group") == "my_group"
        # autocolor defaults to True since no GroupSpec overrides it
        assert kwargs.get("color_group") == "my_group"


class TestCompilerEdgeCases:
    def test_empty_template_compiles(self):
        t = Template()
        plan = TemplateCompiler.compile(t)
        assert plan.tracks == []
        assert plan.resolved_group_indices == {}

    def test_unresolvable_group_track_reference_raises(self):
        t = Template()
        t.tracks = [
            TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="Track A"),
        ]
        t.groups = [
            GroupSpec(name="g1", tracks=["Does Not Exist"], autoscale=True),
        ]
        with pytest.raises(ValueError, match="not found"):
            TemplateCompiler.compile(t)

    def test_numeric_track_reference_resolves(self):
        """Group tracks list supports bare integer index strings."""
        t = Template()
        t.tracks = [
            TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="Track A"),
            TrackSpec(path="b.bw", type=TrackType.BIGWIG, title="Track B"),
        ]
        t.groups = [
            GroupSpec(name="g1", tracks=["0", "1"], autoscale=True),
        ]
        plan = TemplateCompiler.compile(t)
        assert plan.resolved_group_indices["g1"] == [0, 1]

    def test_compile_carries_through_genome_and_guides(self):
        t = Template()
        t.genome = "mm10"
        t.guides.genes = True
        t.guides.axis = False
        t.tracks = []
        plan = TemplateCompiler.compile(t)
        assert plan.genome == "mm10"
        assert plan.add_genes is True
        assert plan.add_axis is False


class TestResolvedTrackToFigureKwargs:
    def test_autoscale_group_in_kwargs(self):
        t = Template()
        t.tracks = [
            TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="T", group="g1"),
        ]
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert kwargs["autoscale_group"] == "g1"

    def test_height_omitted_at_default(self):
        """height=1.0 should not appear in kwargs (default is implicit)."""
        t = Template()
        t.tracks = [TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="T")]
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert "height" not in kwargs

    def test_nondefault_height_in_kwargs(self):
        t = Template()
        t.tracks = [TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="T", height=2.5)]
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert kwargs["height"] == 2.5

    def test_color_in_kwargs_when_set(self):
        t = Template()
        t.tracks = [TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="T", color="red")]
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert kwargs["color"] == "red"

    def test_options_merged_into_kwargs(self):
        t = Template()
        t.tracks = [
            TrackSpec(path="a.bw", type=TrackType.BIGWIG, title="T", options={"alpha": 0.3}),
        ]
        plan = TemplateCompiler.compile(t)
        kwargs = plan.tracks[0].to_figure_kwargs()
        assert kwargs["alpha"] == 0.3


class TestCaseInsensitiveLookup:
    def test_case_mismatch_in_group_tracks_resolves(self):
        """Group track references with wrong case should still resolve."""
        t = Template()
        t.tracks = [
            TrackSpec(path="s1.bw", type=TrackType.BIGWIG, title="H3K27ac (Sample 1)", group="s1"),
        ]
        t.groups = [
            GroupSpec(name="s1", tracks=["H3K27ac (sample 1)"], autoscale=True, autocolor=False),
        ]
        # Should not raise ValueError
        plan = TemplateCompiler.compile(t)
        assert plan.resolved_group_indices["s1"] == [0]
