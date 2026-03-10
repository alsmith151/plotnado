"""Generate static type stub for GenomicFigure hover/autocomplete in IDEs.

This script writes `plotnado/figure.pyi` with explicit keyword-only parameters
for fluent track methods so hover shows available options directly.

Run from repository root:

    python scripts/generate_figure_stub.py
"""

from __future__ import annotations

import types
from enum import Enum
from pathlib import Path
from typing import Any, Union, get_args, get_origin

from pydantic import BaseModel
from pydantic_core import PydanticUndefined

from plotnado.tracks import (
    BigWigCollection,
    BigWigDiff,
    BigWigTrack,
    BedTrack,
    CoolerTrack,
    Genes,
    GenomicAxis,
    HLineTrack,
    HighlightsFromFile,
    LabelConfig,
    LinksTrack,
    NarrowPeakTrack,
    OverlayTrack,
    QuantNadoCoverageTrack,
    QuantNadoMethylationTrack,
    QuantNadoStrandedCoverageTrack,
    QuantNadoVariantTrack,
    ScaleBar,
    Spacer,
    VLineTrack,
)


MODULE_PATH = Path(__file__).resolve().parents[1] / "plotnado" / "figure.pyi"


def _compose_fields(
    track_cls: type[BaseModel],
    *,
    positional_args: set[str],
    include_aesthetics: bool = True,
) -> list[tuple[str, str, str]]:
    track_fields: list[tuple[str, str]] = []
    aesthetics_fields: list[tuple[str, str]] = []
    label_fields: list[tuple[str, str]] = []
    seen: set[str] = set()

    def _annotation_to_str(annotation: Any) -> str:
        if annotation is None:
            return "None"

        origin = get_origin(annotation)
        if origin is None:
            if isinstance(annotation, type):
                if annotation is type(None):
                    return "None"
                module = annotation.__module__
                name = annotation.__name__
                if module == "pandas.core.frame" and name == "DataFrame":
                    return "pd.DataFrame"
                if module in {"builtins", "typing"}:
                    return name
                return name

            rendered = str(annotation)
            rendered = rendered.replace("typing.", "")
            rendered = rendered.replace("types.", "")
            rendered = rendered.replace("pathlib.Path", "Path")
            rendered = rendered.replace("pandas.core.frame.DataFrame", "pd.DataFrame")
            return rendered

        if origin is list:
            args = get_args(annotation)
            arg_str = _annotation_to_str(args[0]) if args else "Any"
            return f"list[{arg_str}]"

        if origin is dict:
            key_type, value_type = get_args(annotation)
            return f"dict[{_annotation_to_str(key_type)}, {_annotation_to_str(value_type)}]"

        if origin is tuple:
            args = get_args(annotation)
            if len(args) == 2 and args[1] is Ellipsis:
                return f"tuple[{_annotation_to_str(args[0])}, ...]"
            return f"tuple[{', '.join(_annotation_to_str(arg) for arg in args)}]"

        if origin is set:
            args = get_args(annotation)
            arg_str = _annotation_to_str(args[0]) if args else "Any"
            return f"set[{arg_str}]"

        if origin in (Union, types.UnionType):
            args = [arg for arg in get_args(annotation)]
            if len(args) == 2 and type(None) in args:
                non_none = args[0] if args[1] is type(None) else args[1]
                return f"{_annotation_to_str(non_none)} | None"
            return " | ".join(_annotation_to_str(arg) for arg in args)

        rendered = str(annotation)
        rendered = rendered.replace("typing.", "")
        rendered = rendered.replace("types.", "")
        rendered = rendered.replace("pathlib.Path", "Path")
        rendered = rendered.replace("pandas.core.frame.DataFrame", "pd.DataFrame")
        return rendered

    def _default_to_str(field_info: Any) -> str:
        if field_info.is_required():
            return "..."

        default = field_info.default
        if default is PydanticUndefined:
            return "..."
        if default is None:
            return "None"
        if isinstance(default, Enum):
            return f"{default.__class__.__name__}.{default.name}"
        if isinstance(default, str):
            return repr(default)
        if isinstance(default, (bool, int, float)):
            return repr(default)
        if isinstance(default, (list, dict, tuple, set)):
            return repr(default)
        return "..."

    for name, field_info in track_cls.model_fields.items():
        if name in {"aesthetics", "label"}:
            continue
        if name in positional_args:
            continue
        if name in seen:
            continue
        track_fields.append(
            (name, _annotation_to_str(field_info.annotation), _default_to_str(field_info))
        )
        seen.add(name)

    if include_aesthetics:
        aesthetics_model = getattr(track_cls, "aesthetics_model")()
        if aesthetics_model is not None:
            for name, field_info in aesthetics_model.model_fields.items():
                if name in seen:
                    continue
                aesthetics_fields.append(
                    (name, _annotation_to_str(field_info.annotation), _default_to_str(field_info))
                )
                seen.add(name)

    for name, field_info in LabelConfig.model_fields.items():
        if name in seen:
            continue
        label_fields.append(
            (name, _annotation_to_str(field_info.annotation), _default_to_str(field_info))
        )
        seen.add(name)

    return [*track_fields, *aesthetics_fields, *label_fields]


def _render_track_method(
    name: str,
    positional_params: list[str],
    kwargs: list[tuple[str, str, str]],
) -> str:
    lines: list[str] = []

    if kwargs:
        kw_lines = [
            f"{kw_name}: {kw_type} = {kw_default},"
            for kw_name, kw_type, kw_default in kwargs
        ]
        if positional_params:
            signature = (
                f"def {name}(self, {', '.join(positional_params)}, /, *,\n"
                + "\n".join(f"        {line}" for line in kw_lines)
                + "\n    ) -> Self: ..."
            )
        else:
            signature = (
                f"def {name}(self, *,\n"
                + "\n".join(f"        {line}" for line in kw_lines)
                + "\n    ) -> Self: ..."
            )
        lines.append(signature)
    else:
        if positional_params:
            lines.append(f"def {name}(self, {', '.join(positional_params)}, /) -> Self: ...")
        else:
            lines.append(f"def {name}(self) -> Self: ...")

    return "\n".join(lines)


def generate() -> str:
    method_specs: list[tuple[str, type[BaseModel], list[str], set[str], bool]] = [
        ("bigwig", BigWigTrack, ["data: Any"], {"data"}, True),
        ("genes", Genes, ["genome: str = ..."], {"genome"}, True),
        ("axis", GenomicAxis, [], set(), True),
        ("scalebar", ScaleBar, [], set(), True),
        ("spacer", Spacer, ["height: float = ..."], {"height"}, False),
        ("bed", BedTrack, ["data: Any"], {"data"}, True),
        ("cooler", CoolerTrack, ["file: str"], {"file"}, True),
        ("bigwig_collection", BigWigCollection, ["files: list[str]"], {"files"}, True),
        ("bigwig_diff", BigWigDiff, ["file_a: str", "file_b: str"], {"file_a", "file_b"}, True),
        ("bigwig_overlay", OverlayTrack, ["tracks: list[Any]"], {"tracks"}, True),
        ("overlay", OverlayTrack, ["tracks: list[Any]"], {"tracks"}, True),
        ("narrowpeak", NarrowPeakTrack, ["data: Any"], {"data"}, True),
        ("links", LinksTrack, ["data: Any"], {"data"}, True),
        ("highlights", HighlightsFromFile, ["data: Any"], {"data"}, True),
        ("hline", HLineTrack, ["y_value: float"], {"y_value"}, True),
        ("vline", VLineTrack, ["x_position: int | str"], {"x_position"}, True),
        ("quantnado_coverage", QuantNadoCoverageTrack, ["sample: str"], {"sample"}, True),
        (
            "quantnado_stranded_coverage",
            QuantNadoStrandedCoverageTrack,
            ["sample: str"],
            {"sample"},
            True,
        ),
        ("quantnado_methylation", QuantNadoMethylationTrack, ["sample: str"], {"sample"}, True),
        ("quantnado_variant", QuantNadoVariantTrack, ["sample: str"], {"sample"}, True),
    ]

    lines: list[str] = [
        "from __future__ import annotations",
        "",
        "from pathlib import Path",
        "from typing import Any, Self",
        "",
        "import matplotlib.figure",
        "import pandas as pd",
        "",
        "from .tracks.enums import (",
        "    BigWigDiffMethod,",
        "    CollectionStyle,",
        "    CoolerTransform,",
        "    DataRangeStyle,",
        "    DisplayMode,",
        "    FontWeight,",
        "    GeneLabelOverlapStrategy,",
        "    GeneLabelStyle,",
        "    NarrowPeakColorBy,",
        "    PlotStyle,",
        "    Position,",
        ")",
        "from .theme import BuiltinTheme, Theme",
        "from .tracks import GenomicRegion, Track",
        "",
        "",
        "class GenomicFigure:",
        "    def __init__(",
        "        self,",
        "        tracks: list[Track] | None = ...,",
        "        width: float = ...,",
        "        track_height: float = ...,",
        "        theme: Theme | BuiltinTheme | str | None = ...,",
        "    ) -> None: ...",
        "",
        "    def add_track(self, track: str | Track, **kwargs: Any) -> Self: ...",
        "",
    ]

    for name, track_cls, positional_params, positional_arg_names, include_aesthetics in method_specs:
        kwargs = _compose_fields(
            track_cls,
            positional_args=positional_arg_names,
            include_aesthetics=include_aesthetics,
        )
        method_block = _render_track_method(name, positional_params, kwargs)
        lines.append("    " + method_block.replace("\n", "\n    "))
        lines.append("")

    lines.extend(
        [
            "    @classmethod",
            "    def available_track_aliases(cls) -> dict[str, str]: ...",
            "",
            "    @classmethod",
            "    def track_options(cls, track: str | type[Track]) -> dict[str, dict]: ...",
            "",
            "    @classmethod",
            "    def track_options_markdown(cls, track: str | type[Track]) -> str: ...",
            "",
            "    def autoscale(self, enable: bool = ...) -> Self: ...",
            "    def autocolor(self, palette: str = ...) -> Self: ...",
            "    def highlight(self, region: str | GenomicRegion) -> Self: ...",
            "    def highlight_style(self, color: str | None = ..., alpha: float | None = ...) -> Self: ...",
            "",
            "    def plot(",
            "        self,",
            "        region: str | GenomicRegion,",
            "        show: bool = ...,",
            "        extend: float | int | None = ...,",
            "        **kwargs: Any,",
            "    ) -> matplotlib.figure.Figure | None: ...",
            "",
            "    def plot_gene(self, gene: str, extend: float = ..., **kwargs: Any) -> matplotlib.figure.Figure | None: ...",
            "",
            "    def plot_regions(",
            "        self,",
            "        regions: list[str] | str,",
            "        ncols: int = ...,",
            "        **kwargs: Any,",
            "    ) -> list[matplotlib.figure.Figure]: ...",
            "",
            "    def to_toml(self, path: str) -> None: ...",
            "",
            "    @classmethod",
            "    def from_toml(cls, path: str) -> GenomicFigure: ...",
            "",
            "    def save(",
            "        self,",
            "        path: str | Path,",
            "        region: str | GenomicRegion,",
            "        dpi: int = ...,",
            "        **kwargs: Any,",
            "    ) -> None: ...",
            "",
            "    def __repr__(self) -> str: ...",
            "    def _repr_html_(self) -> str: ...",
        ]
    )

    return "\n".join(lines) + "\n"


def main() -> None:
    content = generate()
    MODULE_PATH.write_text(content)
    print(f"Wrote {MODULE_PATH}")


if __name__ == "__main__":
    main()
