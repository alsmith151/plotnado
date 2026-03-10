"""Generate TypedDict kwargs types for Figure fluent methods.

This script introspects Plotnado track models and writes `plotnado/_kwargs.py`.
Run from repository root:

    python scripts/generate_kwargs.py
"""

from __future__ import annotations

import types
from pathlib import Path
from typing import Any, Union, get_args, get_origin

from pydantic import BaseModel

from plotnado.tracks import (
    BigWigCollection,
    BigWigDiff,
    BigWigTrack,
    BedTrack,
    CapcruncherTrack,
    CoolerAverage,
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
 2273572867bad7bb6edf1bf8f5ecff6cd4752d5b
    ScaleBar,
    Spacer,
    VLineTrack,
)


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
        if name in {"aesthetics", "label"}:
            continue
        if name in positional_args:
            continue
        fields[name] = _annotation_to_str(field_info.annotation)

    if include_aesthetics:
        aesthetics_model = getattr(track_cls, "aesthetics_model")()
        if aesthetics_model is not None:
            for name, field_type in _collect_fields(aesthetics_model).items():
                if name in fields:
                    continue
                fields[name] = field_type

    for name, field_type in _collect_fields(LabelConfig).items():
        if name in fields:
            continue
        fields[name] = field_type

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
    specs: list[tuple[str, type[BaseModel], set[str], bool]] = [
        ("BigwigKwargs", BigWigTrack, {"data"}, True),
        ("GenesKwargs", Genes, {"genome"}, True),
        ("AxisKwargs", GenomicAxis, set(), True),
        ("ScalebarKwargs", ScaleBar, set(), True),
        ("SpacerKwargs", Spacer, {"height"}, False),
        ("BedKwargs", BedTrack, {"data"}, True),
        ("CoolerKwargs", CoolerTrack, {"file"}, True),
        ("BigwigCollectionKwargs", BigWigCollection, {"files"}, True),
        ("BigwigDiffKwargs", BigWigDiff, {"file_a", "file_b"}, True),
        ("BigwigOverlayKwargs", OverlayTrack, {"tracks"}, True),
        ("OverlayKwargs", OverlayTrack, {"tracks"}, True),
        ("NarrowpeakKwargs", NarrowPeakTrack, {"data"}, True),
        ("LinksKwargs", LinksTrack, {"data"}, True),
        ("HighlightsKwargs", HighlightsFromFile, {"data"}, True),
        ("HlineKwargs", HLineTrack, {"y_value"}, True),
        ("VlineKwargs", VLineTrack, {"x_position"}, True),
        ("CapcruncherKwargs", CapcruncherTrack, {"file"}, True),
        ("CoolerAverageKwargs", CoolerAverage, {"files"}, True),

        ("QuantnadoCoverageKwargs", QuantNadoCoverageTrack, {"sample"}, True),
        ("QuantnadoStrandedCoverageKwargs", QuantNadoStrandedCoverageTrack, {"sample"}, True),
        ("QuantnadoMethylationKwargs", QuantNadoMethylationTrack, {"sample"}, True),
        ("QuantnadoVariantKwargs", QuantNadoVariantTrack, {"sample"}, True),
 2273572867bad7bb6edf1bf8f5ecff6cd4752d5b
    ]

    imports = [
        "from __future__ import annotations",
        "",
        "from pathlib import Path",
        "from typing import Any, TypedDict",
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
        "from .tracks.region import GenomicRegion",
        "from .tracks.base import Track",
        "",
        "",
    ]

    blocks: list[str] = []
    for td_name, track_cls, positional_args, include_aesthetics in specs:
        fields = _compose_fields(
            track_cls,
            positional_args=positional_args,
            include_aesthetics=include_aesthetics,
        )
        blocks.append(_typed_dict_block(td_name, fields))

    content = "\n".join(imports) + "\n\n" + "\n\n".join(blocks) + "\n"
    return content


def main() -> None:
    content = generate()
    MODULE_PATH.write_text(content)
    print(f"Wrote {MODULE_PATH}")


if __name__ == "__main__":
    main()
