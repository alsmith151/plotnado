import tempfile
import importlib.resources
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from plotnado import GenomicFigure
from plotnado.theme import Theme
from plotnado.tracks import (
    BedAesthetics,
    BedTrack,
    BigWigTrack,
    BigwigAesthetics,
    LabelConfig,

    QuantNadoCoverageTrack,
    QuantNadoMethylationTrack,
    QuantNadoStrandedCoverageTrack,
    QuantNadoVariantTrack,
    ScaleBar,
)


class TestFigureRefactor:
    def test_default_figure_uses_publication_theme(self):
        fig = GenomicFigure()

        assert fig.theme == Theme.publication()
        assert fig.highlight_color == Theme.publication().highlight_color

    def test_theme_none_opt_out(self):
        fig = GenomicFigure(theme=None)

        assert fig.theme is None
        assert fig.highlight_color == "#ffd700"
        assert fig.highlight_alpha == 0.15

    def test_default_publication_theme_disables_label_boxes(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )
        fig = GenomicFigure()
        fig.add_track(BigWigTrack(data=df))

        assert fig.tracks[0].label.label_box_enabled is False

    def test_theme_font_family_syncs_label_fonts(self):
        theme = Theme(font_family="Arial")

        assert theme.label.title_font == "Arial"
        assert theme.label.scale_font == "Arial"

    def test_theme_applies_default_track_color(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )
        fig = GenomicFigure(theme=Theme(color="#111111"))
        fig.add_track(BigWigTrack(data=df))

        assert fig.tracks[0].color == "#111111"

    def test_theme_does_not_override_explicit_track_color(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )
        fig = GenomicFigure(theme=Theme(color="#111111"))
        fig.add_track(BigWigTrack(data=df, aesthetics=BigwigAesthetics(color="#ff0000")))

        assert fig.tracks[0].color == "#ff0000"

    def test_theme_controls_highlight_defaults(self):
        fig = GenomicFigure(theme=Theme(highlight_color="#00ff00", highlight_alpha=0.33))
        assert fig.highlight_color == "#00ff00"
        assert fig.highlight_alpha == 0.33

    def test_theme_string_resolves_builtin(self):
        fig = GenomicFigure(theme="minimal")
        assert fig.theme is not None
        assert fig.theme == Theme.minimal()
        assert fig.highlight_color == Theme.minimal().highlight_color

    def test_theme_string_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown builtin theme"):
            GenomicFigure(theme="unknown")

    def test_highlight_style_configuration(self):
        fig = GenomicFigure().highlight_style(color="#123456", alpha=0.2)
        assert fig.highlight_color == "#123456"
        assert fig.highlight_alpha == 0.2

        fig.highlight_style(color="#abcdef", alpha=0.4)
        assert fig.highlight_color == "#abcdef"
        assert fig.highlight_alpha == 0.4

    def test_available_track_aliases_contains_bigwig(self):
        aliases = GenomicFigure.available_track_aliases()
        assert "bigwig" in aliases
        assert aliases["bigwig"] == "BigWigTrack"


    def test_available_track_aliases_contains_quantnado(self):
        aliases = GenomicFigure.available_track_aliases()
        assert aliases["quantnado_coverage"] == "QuantNadoCoverageTrack"
        assert aliases["quantnado_stranded_coverage"] == "QuantNadoStrandedCoverageTrack"
        assert aliases["quantnado_methylation"] == "QuantNadoMethylationTrack"
        assert aliases["quantnado_variant"] == "QuantNadoVariantTrack"
    def test_track_options_by_alias(self):
        options = GenomicFigure.track_options("bigwig")
        assert "aesthetics" in options
        assert "track" in options
        assert "color" in options["aesthetics"]

    def test_track_options_unknown_alias_raises(self):
        with pytest.raises(ValueError, match="Unknown track alias"):
            GenomicFigure.track_options("missing")

    def test_track_options_markdown_by_alias(self):
        markdown = GenomicFigure.track_options_markdown("bigwig")
        assert "## BigWigTrack options" in markdown
        assert "### Aesthetics fields" in markdown
        assert "| color |" in markdown

    def test_alias_constructor_routes_aesthetic_shorthand_kwargs(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure().add_track("bigwig", data=df, color="#ff0000", alpha=0.25)
        track = fig.tracks[0]

        assert track.aesthetics.color == "#ff0000"
        assert track.aesthetics.alpha == 0.25

    def test_alias_constructor_routes_label_shorthand_kwargs(self):
        fig = GenomicFigure().add_track("scalebar", plot_title=False, scale_size=11)
        track = fig.tracks[0]

        assert track.label.plot_title is False
        assert track.label.scale_size == 11


    def test_quantnado_alias_constructor_routes_aesthetic_shorthand_kwargs(self):
        xr = pytest.importorskip("xarray")
        coverage_data = xr.DataArray(
            [[1.0, 2.0, 3.0]],
            dims=("sample", "position"),
            coords={"sample": ["s1"], "position": [100, 101, 102]},
        )

        fig = GenomicFigure(theme=None).add_track(
            "quantnado_coverage",
            sample="s1",
            coverage_data=coverage_data,
            color="#ff00ff",
            alpha=0.4,
        )
        track = fig.tracks[0]

        assert isinstance(track, QuantNadoCoverageTrack)
        assert track.color == "#ff00ff"
        assert track.alpha == 0.4
    def test_explicit_nested_model_and_shorthand_merge(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure().add_track(
            "bigwig",
            data=df,
            aesthetics=BigwigAesthetics(color="#00aa00", alpha=0.5),
            alpha=0.2,
        )
        track = fig.tracks[0]

        assert track.aesthetics.color == "#00aa00"
        assert track.aesthetics.alpha == 0.2


    def test_quantnado_helper_methods_append_expected_tracks(self):
        xr = pytest.importorskip("xarray")
        base = xr.DataArray(
            [[1.0, 2.0, 3.0]],
            dims=("sample", "position"),
            coords={"sample": ["s1"], "position": [100, 101, 102]},
        )
        fig = GenomicFigure(theme=None)
        fig.quantnado_coverage("s1", coverage_data=base)
        fig.quantnado_stranded_coverage("s1", coverage_fwd_data=base, coverage_rev_data=base)
        fig.quantnado_methylation("s1", methylation_data=base)
        fig.quantnado_variant("s1", allele_depth_ref_data=base, allele_depth_alt_data=base)

        assert isinstance(fig.tracks[0], QuantNadoCoverageTrack)
        assert isinstance(fig.tracks[1], QuantNadoStrandedCoverageTrack)
        assert isinstance(fig.tracks[2], QuantNadoMethylationTrack)
        assert isinstance(fig.tracks[3], QuantNadoVariantTrack)
    def test_extend_parameter(self):
        fig = GenomicFigure().add_track(ScaleBar())
        out = fig.plot("chr1:100-200", show=False, extend=0.5)

        xlim = out.axes[0].get_xlim()
        assert xlim[0] == 50
        assert xlim[1] == 250

    def test_plot_regions_list(self):
        fig = GenomicFigure().add_track(ScaleBar())
        plots = fig.plot_regions(["chr1:100-200", "chr1:300-400"], show=False)

        assert len(plots) == 2

    def test_plot_regions_grid_mode(self):
        fig = GenomicFigure().add_track(ScaleBar())
        plots = fig.plot_regions(
            ["chr1:100-200", "chr1:300-400", "chr1:500-600"],
            ncols=2,
            show=False,
        )

        assert len(plots) == 1
        assert len(plots[0].axes) == 4

    def test_plot_regions_invalid_ncols(self):
        fig = GenomicFigure().add_track(ScaleBar())

        with pytest.raises(ValueError, match="ncols must be >= 1"):
            fig.plot_regions(["chr1:100-200"], ncols=0, show=False)

    def test_plot_regions_bed_path(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as handle:
            handle.write("chr1\t100\t200\n")
            handle.write("chr1\t300\t400\n")
            bed_path = handle.name

        fig = GenomicFigure().add_track(ScaleBar())
        plots = fig.plot_regions(bed_path, show=False)

        assert len(plots) == 2
        Path(bed_path).unlink(missing_ok=True)

    def test_autoscale_group(self):
        df1 = pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [100, 150], "end": [150, 200], "value": [1.0, 2.0]}
        )
        df2 = pd.DataFrame(
            {"chrom": ["chr1", "chr1"], "start": [100, 150], "end": [150, 200], "value": [10.0, 20.0]}
        )

        fig = GenomicFigure()
        fig.add_track(BigWigTrack(data=df1, autoscale_group="g1"))
        fig.add_track(BigWigTrack(data=df2, autoscale_group="g1"))

        out = fig.plot("chr1:90-210", show=False)
        ylim1 = out.axes[0].get_ylim()
        ylim2 = out.axes[1].get_ylim()
        assert ylim1 == ylim2

    def test_publication_theme_auto_palette_applies_to_multiple_signal_tracks(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure()
        fig.add_track(BigWigTrack(data=df))
        fig.add_track(BigWigTrack(data=df))

        palette = Theme.publication().palette
        assert palette is not None
        assert fig.tracks[0].color == palette[0]
        assert fig.tracks[1].color == palette[1]

    def test_publication_theme_auto_palette_respects_explicit_color(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure()
        fig.add_track(BigWigTrack(data=df, aesthetics=BigwigAesthetics(color="#101010")))
        fig.add_track(BigWigTrack(data=df))

        palette = Theme.publication().palette
        assert palette is not None
        assert fig.tracks[0].color == "#101010"
        assert fig.tracks[1].color == palette[0]

    def test_publication_theme_auto_palette_respects_color_group(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure()
        fig.add_track(BigWigTrack(data=df, color_group="histone"))
        fig.add_track(BigWigTrack(data=df, color_group="histone"))
        fig.add_track(BigWigTrack(data=df, color_group="atac"))

        assert fig.tracks[0].color == fig.tracks[1].color
        assert fig.tracks[2].color != fig.tracks[0].color

    def test_autocolor_applies_when_enabled_before_add_track(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure().autocolor("tab10")
        fig.add_track(BigWigTrack(data=df))

        expected = mcolors.to_hex(plt.get_cmap("tab10")(0))
        assert fig.tracks[0].color == expected

    def test_autocolor_without_palette_uses_theme_palette(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure(theme="publication")
        fig.add_track(BigWigTrack(data=df))
        fig.add_track(BigWigTrack(data=df))
        fig.autocolor()

        palette = Theme.publication().palette
        assert palette is not None
        assert fig.tracks[0].color == palette[0]
        assert fig.tracks[1].color == palette[1]

    def test_autocolor_respects_color_group(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )
        fig = GenomicFigure(theme=None)
        fig.add_track(BigWigTrack(data=df, color_group="group_a"))
        fig.add_track(BigWigTrack(data=df, color_group="group_a"))
        fig.add_track(BigWigTrack(data=df, color_group="group_b"))

        fig.autocolor("tab10")

        assert fig.tracks[0].color == fig.tracks[1].color
        assert fig.tracks[2].color != fig.tracks[0].color

    def test_autocolor_overrides_theme_default_for_new_tracks(self):
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )

        fig = GenomicFigure(theme="publication").autocolor("tab10")
        fig.add_track(BigWigTrack(data=df))

        expected = mcolors.to_hex(plt.get_cmap("tab10")(0))
        assert fig.tracks[0].color == expected

    def test_autocolor_skips_meta_tracks_and_keeps_data_palette_contiguous(self):
        signal_df = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1"],
                "start": [100, 200],
                "end": [200, 300],
                "value": [1.0, 2.0],
            }
        )
        highlight_df = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [120],
                "end": [180],
            }
        )

        fig = GenomicFigure().autocolor("tab10")
        fig.add_track("axis")
        fig.add_track("scale")
        fig.add_track("highlight", data=highlight_df)
        fig.add_track(BigWigTrack(data=signal_df))
        fig.add_track(BigWigTrack(data=signal_df))

        cmap = plt.get_cmap("tab10")
        expected_first = mcolors.to_hex(cmap(0))
        expected_second = mcolors.to_hex(cmap(1))

        assert fig.tracks[0].color == "#666666"
        assert fig.tracks[1].color == "#333333"
        assert fig.tracks[2].color == "yellow"
        assert fig.tracks[3].color == expected_first
        assert fig.tracks[4].color == expected_second

    def test_toml_roundtrip(self):
        pytest.importorskip("tomli_w")
        fig = GenomicFigure().add_track(ScaleBar())

        with tempfile.NamedTemporaryFile(suffix=".toml", delete=False) as handle:
            path = handle.name

        fig.to_toml(path)
        loaded = GenomicFigure.from_toml(path)

        assert isinstance(loaded, GenomicFigure)
        assert len(loaded.tracks) == 1
        assert loaded.tracks[0].__class__.__name__ == "ScaleBar"
        Path(path).unlink(missing_ok=True)

    def test_toml_uses_track_tables(self):
        pytest.importorskip("tomli_w")
        fig = GenomicFigure().add_track(ScaleBar())

        with tempfile.NamedTemporaryFile(suffix=".toml", delete=False) as handle:
            path = handle.name

        fig.to_toml(path)
        content = Path(path).read_text()

        assert "[[tracks.ScaleBar]]" in content
        Path(path).unlink(missing_ok=True)

    def test_plot_gene_extend_math_with_bundled_hg38(self):
        bed_prefix = importlib.resources.files("plotnado.data.gene_bed_files")
        with open(bed_prefix / "genes.json") as handle:
            mapping = json.load(handle)

        genes_df = pd.read_csv(bed_prefix / mapping["hg38"], sep="\t", header=None)
        genes_df.columns = [
            "chrom",
            "start",
            "end",
            "name",
            *[f"field_{i}" for i in range(max(0, genes_df.shape[1] - 4))],
        ]

        gene = str(genes_df.iloc[0]["name"])
        match = genes_df.loc[genes_df["name"].astype(str).str.upper() == gene.upper()].iloc[0]
        expected_length = int(match["end"]) - int(match["start"])
        expected_extend = int(expected_length * 0.5)

        fig = GenomicFigure().add_track(ScaleBar())
        out = fig.plot_gene(gene, extend=0.5, show=False)

        xlim = out.axes[0].get_xlim()
        assert xlim[0] == int(match["start"]) - expected_extend
        assert xlim[1] == int(match["end"]) + expected_extend

    def test_plot_gene_adds_genes_track_when_missing(self):
        fig = GenomicFigure().add_track(ScaleBar())
        out = fig.plot_gene("DDX11L1", show=False)

        assert out is not None
        assert len(out.axes) >= 2

    def test_plot_normalizes_font_family_across_text_artists(self):
        data = pd.DataFrame(
            {
                "chrom": ["chr1"],
                "start": [100],
                "end": [200],
                "name": ["interval_1"],
            }
        )
        theme = Theme(font_family="Arial", label=LabelConfig(scale_font="Times New Roman"))
        fig = GenomicFigure(theme=theme)
        fig.add_track(BedTrack(data=data, aesthetics=BedAesthetics(show_labels=True)))

        out = fig.plot("chr1:50-250", show=False)

        fonts = []
        for axis in out.axes:
            for text_artist in axis.texts:
                families = text_artist.get_fontfamily()
                if isinstance(families, (list, tuple)):
                    fonts.extend([str(font) for font in families])
                else:
                    fonts.append(str(families))

        assert fonts
        assert set(fonts) == {"Arial"}

    def test_save_default_dpi_is_600(self):
        fig = GenomicFigure().add_track(ScaleBar())

        with patch.object(fig, "plot") as mock_plot:
            mock_matplotlib_figure = plt.figure()
            mock_plot.return_value = mock_matplotlib_figure
            with patch.object(mock_matplotlib_figure, "savefig") as mock_savefig:
                fig.save("out.png", "chr1:100-200")
                mock_savefig.assert_called_once()
                assert mock_savefig.call_args.kwargs["dpi"] == 600
            plt.close(mock_matplotlib_figure)
