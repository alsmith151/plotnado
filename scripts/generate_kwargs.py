"""Generate TypedDict kwargs types for GenomicFigure builder methods.

This script introspects the track registry and the same flattened kwarg routing
used by ``GenomicFigure._create_track_from_alias()`` to write
``plotnado/_kwargs.py``.

Run from repository root:

    python scripts/generate_kwargs.py
"""

from __future__ import annotations

import types
from pathlib import Path
from typing import Any, Literal, Union, get_args, get_origin

from pydantic import BaseModel

import plotnado.tracks  # noqa: F401
from plotnado.tracks import LabelConfig
from plotnado.tracks.enums import TrackType
from plotnado.tracks.registry import registry


MODULE_PATH = Path(__file__).resolve().parents[1] / "plotnado" / "_kwargs.py"


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


def _collect_fields(model: type[BaseModel]) -> dict[str, str]:
    return {
        name: _annotation_to_str(field_info.annotation)
        for name, field_info in model.model_fields.items()
    }


def _compose_fields(
    track_cls: type[BaseModel],
    *,
    positional_args: set[str],
    include_aesthetics: bool = True,
) -> dict[str, str]:
    fields: dict[str, str] = {}

    for name, field_info in track_cls.model_fields.items():
        if name in {"aesthetics", "label"} or name in positional_args:
            continue
        fields[name] = _annotation_to_str(field_info.annotation)

    if include_aesthetics:
        aesthetics_model = getattr(track_cls, "aesthetics_model")()
        if aesthetics_model is not None:
            for name, field_type in _collect_fields(aesthetics_model).items():
                fields.setdefault(name, field_type)

    for name, field_type in _collect_fields(LabelConfig).items():
        fields.setdefault(name, field_type)

    return fields


def _typed_dict_block(name: str, fields: dict[str, str]) -> str:
    lines = [f"class {name}(TypedDict, total=False):"]
    if not fields:
        lines.append("    pass")
        return "\n".join(lines)

    for field_name, field_type in fields.items():
        lines.append(f"    {field_name}: {field_type}")
    return "\n".join(lines)


def generate() -> str:
    method_specs: list[tuple[str, TrackType | str, set[str], bool]] = [
        ("BigwigKwargs", TrackType.BIGWIG, {"data"}, True),
        ("GenesKwargs", TrackType.GENE, {"genome"}, True),
        ("AxisKwargs", TrackType.AXIS, set(), True),
        ("ScalebarKwargs", TrackType.SCALEBAR, set(), True),
        ("SpacerKwargs", TrackType.SPACER, {"height"}, False),
        ("BedKwargs", TrackType.BED, {"data"}, True),
        ("CoolerKwargs", TrackType.COOLER, {"file"}, True),
        ("CapcruncherKwargs", TrackType.CAPCRUNCHER, {"file"}, True),
        ("CoolerAverageKwargs", TrackType.COOLER_AVERAGE, {"files"}, True),
        ("BigwigCollectionKwargs", TrackType.BIGWIG_COLLECTION, {"files"}, True),
        ("BigwigDiffKwargs", TrackType.BIGWIG_DIFF, {"file_a", "file_b"}, True),
        ("BigwigOverlayKwargs", "bigwig_overlay", {"tracks"}, True),
        ("OverlayKwargs", TrackType.OVERLAY, {"tracks"}, True),
        ("NarrowpeakKwargs", TrackType.NARROWPEAK, {"data"}, True),
        ("LinksKwargs", TrackType.LINKS, {"data"}, True),
        ("HighlightsKwargs", TrackType.HIGHLIGHT, {"data"}, True),
        ("HlineKwargs", TrackType.HLINE, {"y_value"}, True),
        ("VlineKwargs", TrackType.VLINE, {"x_position"}, True),
        ("QuantnadoCoverageKwargs", TrackType.QUANTNADO_COVERAGE, {"sample"}, True),
        (
            "QuantnadoStrandedCoverageKwargs",
            TrackType.QUANTNADO_STRANDED_COVERAGE,
            {"sample"},
            True,
        ),
        ("QuantnadoMethylationKwargs", TrackType.QUANTNADO_METHYLATION, {"sample"}, True),
        ("QuantnadoVariantKwargs", TrackType.QUANTNADO_VARIANT, {"sample"}, True),
    ]

    imports = [
        "# AUTO-GENERATED - do not edit manually.",
        "# Re-generate with: python scripts/generate_kwargs.py",
        "#",
        "# This file provides TypedDict definitions for GenomicFigure builder",
        "# methods, enabling IDE autocompletion and type checking.",
        "from __future__ import annotations",
        "",
        "from pathlib import Path",
        "from typing import Any, Literal, TypedDict",
        "",
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
        "",
    ]

    blocks: list[str] = []
    for td_name, key, positional_args, include_aesthetics in method_specs:
        entry = registry.get(key.value if isinstance(key, TrackType) else key)
        fields = _compose_fields(
            entry.cls,
            positional_args=positional_args,
            include_aesthetics=include_aesthetics,
        )
        blocks.append(_typed_dict_block(td_name, fields))

    return "\n".join(imports) + "\n\n" + "\n\n".join(blocks) + "\n"


def main() -> None:
    MODULE_PATH.write_text(generate())
    print(f"Wrote {MODULE_PATH}")


if __name__ == "__main__":
    main()
