"""Generate code and track-option reference pages."""

from pathlib import Path

import mkdocs_gen_files

from plotnado.figure import Figure

nav = mkdocs_gen_files.Nav()

for path in sorted(Path("plotnado").rglob("*.py")):

    module_path = path.relative_to(".").with_suffix("")
    doc_path = module_path.with_suffix(".md")
    full_doc_path = Path("reference", doc_path)

    parts = list(module_path.parts)

    if module_path.stem == "__init__":
        continue
    elif parts[-1] == "__init__":
        parts = parts[:-1]
        doc_path = doc_path.with_name("index.md")
        full_doc_path = full_doc_path.with_name("index.md")
    elif parts[-1] == "__main__":
        continue

    nav[parts] = doc_path.as_posix()

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        identifier = ".".join(parts)
        print("::: " + identifier, file=fd)

    mkdocs_gen_files.set_edit_path(full_doc_path, path)

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())


def _format_section_rows(options: dict[str, dict], section: str) -> list[str]:
    rows = [
        f"### {section.title()} fields",
        "",
        "| Name | Type | Default | Choices | Required | Description |",
        "|---|---|---|---|---|---|",
    ]
    section_data = options.get(section, {})
    for field_name, meta in sorted(section_data.items()):
        description = meta.get("description") or "—"
        choices = meta.get("choices") or []
        choices_text = ", ".join(str(choice) for choice in choices) if choices else "—"
        rows.append(
            f"| {field_name} | {meta['type']} | {meta['default']} | {choices_text} | {meta['required']} | {description} |"
        )
    rows.append("")
    return rows


def _generate_aesthetics_reference() -> None:
    alias_map = Figure.available_track_aliases()
    class_by_aliases: dict[str, list[str]] = {}
    class_by_name: dict[str, type] = {}
    for alias, class_name in alias_map.items():
        class_by_aliases.setdefault(class_name, []).append(alias)
        track_cls = Figure._alias_map()[alias]
        class_by_name[class_name] = track_cls

    lines = [
        "# Aesthetics Reference",
        "",
        "This page is auto-generated from runtime model metadata (`Track.options()`).",
        "",
        "Track styles are configured through nested `aesthetics=...` models.",
        "",
    ]

    for class_name in sorted(class_by_name):
        track_cls = class_by_name[class_name]
        aliases = ", ".join(f"`{alias}`" for alias in sorted(class_by_aliases[class_name]))
        lines.append(f"## {class_name}")
        lines.append("")
        lines.append(f"Aliases: {aliases}")
        lines.append("")
        options = track_cls.options()
        lines.extend(_format_section_rows(options, "track"))
        lines.extend(_format_section_rows(options, "aesthetics"))
        lines.extend(_format_section_rows(options, "label"))

    with mkdocs_gen_files.open("aesthetics_reference.md", "w") as handle:
        handle.write("\n".join(lines))


_generate_aesthetics_reference()
