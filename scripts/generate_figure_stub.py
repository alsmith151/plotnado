"""Generate static type stub for GenomicFigure hover/autocomplete in IDEs.

This script writes `plotnado/figure.pyi` with explicit keyword-only parameters
for fluent track methods so hover shows available options directly.

Run from repository root:

    python scripts/generate_figure_stub.py

Strategy
--------
* Track-builder methods (bigwig, genes, …) — rich parameter lists are generated
  by introspecting Pydantic model fields via ``_compose_fields()``.
* All other public methods — discovered automatically from the live
  ``GenomicFigure`` class via ``inspect`` + ``typing.get_type_hints()``, so
  new methods added to ``figure.py`` appear in the stub without requiring any
  changes here.
* A small set of methods with complex / non-inferrable signatures (``__init__``,
  ``from_igv_session``, ``update_track``, etc.) are still written as literal
  strings in the header section of ``generate()``.
"""

from __future__ import annotations

import inspect
import types
import typing
from enum import Enum
from pathlib import Path
from typing import Any, Literal, Union, get_args, get_origin

from pydantic import BaseModel
from pydantic_core import PydanticUndefined

from plotnado.tracks import (
    BigWigCollection,
    BigWigDiff,
    BigWigTrack,
    BedTrack,
    CapcruncherTrack,
    CoolerTrack,
    CoolerAverage,
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

# Imported here so _render_sig_stub can introspect it at generation time.
from plotnado.figure import GenomicFigure
import plotnado.figure as _figure_module


MODULE_PATH = Path(__file__).resolve().parents[1] / "plotnado" / "figure.pyi"


# ---------------------------------------------------------------------------
# Shared annotation-to-string helper (used by both Pydantic and inspect paths)
# ---------------------------------------------------------------------------

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
            if module in {"pandas", "pandas.core.frame"} and name == "DataFrame":
                return "pd.DataFrame"
            if module in {"pandas", "pandas.core.series"} and name == "Series":
                return "pd.Series"
            if module.startswith("matplotlib"):
                return f"{module}.{name}"
            if module in {"builtins", "typing"}:
                return name
            return name

        rendered = str(annotation)
        rendered = rendered.replace("typing.", "")
        rendered = rendered.replace("types.", "")
        rendered = rendered.replace("pathlib.Path", "Path")
        rendered = rendered.replace("DataFrame", "pd.DataFrame")
        rendered = rendered.replace("pandas.DataFrame", "pd.DataFrame")
        rendered = rendered.replace("pandas.core.frame.DataFrame", "pd.DataFrame")
        return rendered

    if origin is type:
        args = get_args(annotation)
        return f"type[{_annotation_to_str(args[0])}]" if args else "type"

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
    rendered = rendered.replace("DataFrame", "pd.DataFrame")
    rendered = rendered.replace("pandas.DataFrame", "pd.DataFrame")
    rendered = rendered.replace("pandas.core.frame.DataFrame", "pd.DataFrame")
    return rendered


# ---------------------------------------------------------------------------
# Pydantic field helpers (track-builder path)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# inspect-based helpers (non-track method path)
# ---------------------------------------------------------------------------

def _default_repr(default: Any) -> str:
    if default is None:
        return "None"
    if isinstance(default, Enum):
        return f"{default.__class__.__name__}.{default.name}"
    if isinstance(default, (bool, int, float, str)):
        return repr(default)
    return "..."


def _render_sig_stub(name: str) -> list[str]:
    """Emit stub lines for a single method on GenomicFigure using inspect."""
    raw = inspect.getattr_static(GenomicFigure, name)
    is_classmethod = isinstance(raw, classmethod)

    func = raw.__func__ if is_classmethod else raw

    try:
        hints = typing.get_type_hints(func, globalns=vars(_figure_module))
    except Exception:
        hints = {}

    method = getattr(GenomicFigure, name)
    try:
        sig = inspect.signature(method)
    except (ValueError, TypeError):
        return [f"    def {name}(self, *args: Any, **kwargs: Any) -> Any: ..."]

    params: list[str] = ["cls" if is_classmethod else "self"]
    has_star = False

    for pname, param in sig.parameters.items():
        # classmethod: cls already bound, excluded from sig
        # instance method: self present in sig — skip it, we added it above
        if pname in ("self", "cls"):
            continue
        ann = hints.get(pname)
        ann_str = f": {_annotation_to_str(ann)}" if ann is not None else ""

        if param.kind == inspect.Parameter.VAR_POSITIONAL:
            params.append(f"*{pname}{ann_str}")
            has_star = True
        elif param.kind == inspect.Parameter.VAR_KEYWORD:
            params.append(f"**{pname}{ann_str}")
        elif param.kind == inspect.Parameter.KEYWORD_ONLY:
            if not has_star:
                params.append("*")
                has_star = True
            default = param.default
            if default is inspect.Parameter.empty:
                params.append(f"{pname}{ann_str}")
            else:
                params.append(f"{pname}{ann_str} = {_default_repr(default)}")
        else:
            default = param.default
            if default is inspect.Parameter.empty:
                params.append(f"{pname}{ann_str}")
            else:
                params.append(f"{pname}{ann_str} = {_default_repr(default)}")

    return_ann = hints.get("return")
    return_str = f" -> {_annotation_to_str(return_ann)}" if return_ann is not None else ""

    output: list[str] = []
    if is_classmethod:
        output.append("    @classmethod")

    joined = ", ".join(params)
    one_liner = f"    def {name}({joined}){return_str}: ..."
    if len(one_liner) <= 100:
        output.append(one_liner)
    else:
        output.append(f"    def {name}(")
        for part in params:
            output.append(f"        {part},")
        output.append(f"    ){return_str}: ...")

    return output


def _collect_auto_method_stubs(*, skip: set[str]) -> list[str]:
    """Discover and render stubs for all public GenomicFigure methods not in skip."""
    # Include these even though they start with '_'
    dunder_include = {"__repr__", "_repr_html_"}

    lines: list[str] = []
    seen_names: list[str] = []

    for name in dir(GenomicFigure):
        if name in skip:
            continue
        if name.startswith("_") and name not in dunder_include:
            continue
        raw = inspect.getattr_static(GenomicFigure, name)
        if not callable(raw if not isinstance(raw, (classmethod, staticmethod)) else raw.__func__):
            continue
        seen_names.append(name)

    for name in seen_names:
        lines.extend(_render_sig_stub(name))
        lines.append("")

    return lines


# ---------------------------------------------------------------------------
# Main generator
# ---------------------------------------------------------------------------

def generate() -> str:
    method_specs: list[tuple[str, type[BaseModel], list[str], set[str], bool]] = [
        ("bigwig", BigWigTrack, ["data: Any"], {"data"}, True),
        ("genes", Genes, ["genome: str = ..."], {"genome"}, True),
        ("axis", GenomicAxis, [], set(), True),
        ("scalebar", ScaleBar, [], set(), True),
        ("spacer", Spacer, ["height: float = ..."], {"height"}, False),
        ("bed", BedTrack, ["data: Any"], {"data"}, True),
        ("cooler", CoolerTrack, ["file: str"], {"file"}, True),
        ("capcruncher", CapcruncherTrack, ["file: str"], {"file"}, True),
        ("cooler_average", CoolerAverage, ["files: list[str]"], {"files"}, True),
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

    track_method_names = {name for name, *_ in method_specs}

    # Methods with complex/non-inferrable signatures written as literal strings below.
    hand_written_names = {
        "__init__",
        "from_template",
        "from_igv_session",
        "__getitem__",
        "update_track",
        "update_tracks",
        "remove_track",
        "add_track",
    }

    lines: list[str] = [
        "from __future__ import annotations",
        "",
        "from pathlib import Path",
        "from typing import Any, Literal, Self",
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
        "    @classmethod",
        "    def from_template(",
        "        cls,",
        "        path: str | Path,",
        "        *,",
        "        width: float | None = ...,",
        "        theme: Theme | BuiltinTheme | str | None = ...,",
        "    ) -> GenomicFigure: ...",
        "",
        "    @classmethod",
        "    def from_igv_session(",
        "        cls,",
        "        path: str | Path,",
        "        *,",
        "        width: float | None = ...,",
        "        theme: Theme | BuiltinTheme | str | None = ...,",
        "    ) -> tuple[GenomicFigure, str | None]: ...",
        "",
        "    def __getitem__(self, key: int | str) -> Track: ...",
        "",
        "    def update_track(",
        "        self,",
        "        key: int | str | None = ...,",
        "        *,",
        "        track_type: str | type[Track] | None = ...,",
        "        group: str | None = ...,",
        "        where: Any = ...,",
        "        **kwargs: Any,",
        "    ) -> Self: ...",
        "",
        "    def update_tracks(",
        "        self,",
        "        *,",
        "        track_type: str | type[Track] | None = ...,",
        "        group: str | None = ...,",
        "        where: Any = ...,",
        "        **kwargs: Any,",
        "    ) -> Self: ...",
        "",
        "    def remove_track(self, key: int | str) -> Self: ...",
        "",
        "    def add_track(self, track: str | Track, *, position: str = ..., **kwargs: Any) -> Self: ...",
        "",
    ]

    # --- Track-builder methods (Pydantic introspection) ---
    for name, track_cls, positional_params, positional_arg_names, include_aesthetics in method_specs:
        kwargs = _compose_fields(
            track_cls,
            positional_args=positional_arg_names,
            include_aesthetics=include_aesthetics,
        )
        method_block = _render_track_method(name, positional_params, kwargs)
        lines.append("    " + method_block.replace("\n", "\n    "))
        lines.append("")

    # --- All remaining public methods (auto-detected via inspect) ---
    skip = track_method_names | hand_written_names
    auto_lines = _collect_auto_method_stubs(skip=skip)
    lines.extend("    " + ln if ln and not ln.startswith("    ") else ln for ln in auto_lines)

    return "\n".join(lines) + "\n"


def main() -> None:
    content = generate()
    MODULE_PATH.write_text(content)
    print(f"Wrote {MODULE_PATH}")


if __name__ == "__main__":
    main()
