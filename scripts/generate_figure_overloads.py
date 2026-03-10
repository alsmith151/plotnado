"""Generate overload blocks in `plotnado/figure.py` for IDE hover signatures.

This script injects auto-generated `@overload` blocks above fluent Figure methods.
The overloads are generated from track + aesthetics + label model fields so keyword
names stay synchronized automatically.

Run from repository root:

    python scripts/generate_figure_overloads.py
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
    CoolerTrack,
    Genes,
    GenomicAxis,
    HLineTrack,
    HighlightsFromFile,
    LabelConfig,
    LinksTrack,
    NarrowPeakTrack,
    OverlayTrack,
<<<<<<< HEAD
    QuantNadoCoverageTrack,
    QuantNadoMethylationTrack,
    QuantNadoStrandedCoverageTrack,
    QuantNadoVariantTrack,
=======
>>>>>>> 2273572867bad7bb6edf1bf8f5ecff6cd4752d5b
    ScaleBar,
    Spacer,
    VLineTrack,
)


FIGURE_PATH = Path(__file__).resolve().parents[1] / "plotnado" / "figure.py"


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


def _compose_fields(
    track_cls: type[BaseModel],
    *,
    positional_args: set[str],
    include_aesthetics: bool = True,
) -> list[tuple[str, str]]:
    track_fields: list[tuple[str, str]] = []
    aesthetics_fields: list[tuple[str, str]] = []
    label_fields: list[tuple[str, str]] = []
    seen: set[str] = set()

    for name, field_info in sorted(track_cls.model_fields.items()):
        if name in {"aesthetics", "label"}:
            continue
        if name in positional_args:
            continue
        if name in seen:
            continue
        track_fields.append((name, _annotation_to_str(field_info.annotation)))
        seen.add(name)

    if include_aesthetics:
        aesthetics_model = getattr(track_cls, "aesthetics_model")()
        if aesthetics_model is not None:
            for name, field_info in sorted(aesthetics_model.model_fields.items()):
                if name in seen:
                    continue
                aesthetics_fields.append((name, _annotation_to_str(field_info.annotation)))
                seen.add(name)

    for name, field_info in sorted(LabelConfig.model_fields.items()):
        if name in seen:
            continue
        label_fields.append((name, _annotation_to_str(field_info.annotation)))
        seen.add(name)

    return [*track_fields, *aesthetics_fields, *label_fields]


def _render_method_block(
    name: str,
    positional_params: list[str],
    kwargs: list[tuple[str, str]],
) -> str:
    lines: list[str] = []
    marker_start = f"    # BEGIN AUTO-GENERATED OVERLOAD: {name}"
    marker_end = f"    # END AUTO-GENERATED OVERLOAD: {name}"

    lines.append(marker_start)
    lines.append("    @overload")
    lines.append(f"    def {name}(")
    lines.append("        self,")

    for param in positional_params:
        lines.append(f"        {param},")

    if positional_params:
        lines.append("        /,")

    lines.append("        *,")
    for kw_name, kw_type in kwargs:
        lines.append(f"        {kw_name}: {kw_type} = ...,")
    lines.append("    ) -> Self: ...")
    lines.append(marker_end)
    return "\n".join(lines)


def generate_blocks() -> dict[str, str]:
    method_specs: list[tuple[str, type[BaseModel], list[str], set[str], bool]] = [
        ("bigwig", BigWigTrack, ["data: Any"], {"data"}, True),
        ("genes", Genes, ["genome: str = \"hg38\""], {"genome"}, True),
        ("axis", GenomicAxis, [], set(), True),
        ("scalebar", ScaleBar, [], set(), True),
        ("spacer", Spacer, ["height: float = 0.5"], {"height"}, False),
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
<<<<<<< HEAD
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
=======
>>>>>>> 2273572867bad7bb6edf1bf8f5ecff6cd4752d5b
    ]

    blocks: dict[str, str] = {}
    for name, track_cls, positional_params, positional_arg_names, include_aesthetics in method_specs:
        kwargs = _compose_fields(
            track_cls,
            positional_args=positional_arg_names,
            include_aesthetics=include_aesthetics,
        )
        blocks[name] = _render_method_block(name, positional_params, kwargs)

    return blocks


def _remove_manual_overload_block(source: str, method_name: str) -> str:
    signature = f"    @overload\n    def {method_name}("
    if signature not in source:
        return source

    start = source.index(signature)
    marker = "    ) -> Self: ...\n\n"
    end = source.find(marker, start)
    if end == -1:
        return source
    end += len(marker)
    return source[:start] + source[end:]


def apply_to_source(source: str) -> str:
    updated = source
    blocks = generate_blocks()

    for method_name, block in blocks.items():
        updated = _remove_manual_overload_block(updated, method_name)
        marker_start = f"    # BEGIN AUTO-GENERATED OVERLOAD: {method_name}"
        marker_end = f"    # END AUTO-GENERATED OVERLOAD: {method_name}"
        method_def = f"    def {method_name}("

        if marker_start in updated and marker_end in updated:
            start = updated.index(marker_start)
            end = updated.index(marker_end, start) + len(marker_end)
            updated = updated[:start] + block + updated[end:]
            continue

        insert_at = updated.find(method_def)
        if insert_at == -1:
            raise RuntimeError(f"Could not locate method definition for {method_name}")

        updated = updated[:insert_at] + block + "\n" + updated[insert_at:]

    return updated


def generate() -> str:
    source = FIGURE_PATH.read_text()
    return apply_to_source(source)


def main() -> None:
    content = generate()
    FIGURE_PATH.write_text(content)
    print(f"Wrote {FIGURE_PATH}")


if __name__ == "__main__":
    main()
