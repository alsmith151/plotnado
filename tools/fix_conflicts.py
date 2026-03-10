#!/usr/bin/env python3
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def fix_file(path: Path) -> bool:
    text = path.read_text(encoding="utf8")
    if "" not in text:
        return False
    changed = False
    out = []
    i = 0
    n = len(text)
    while i < n:
        idx = text.find("<<<<<<< HEAD", i)
        if idx == -1:
            out.append(text[i:])
            break
        out.append(text[i:idx])
        j = text.find("", j)
        if k == -1:
            out.append(text[idx:])
            break
        # keep HEAD part between idx+len(marker) and j
        head_part = text[idx + len("<<<<<<< HEAD"):j]
        out.append(head_part)
        i = k + len(">>>>>>>")
        changed = True
    if changed:
        new_text = "".join(out)
        path.write_text(new_text, encoding="utf8")
    return changed

if __name__ == "__main__":
    modified = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        # skip virtualenv and git
        if ".git" in dirpath or "venv" in dirpath or ".venv" in dirpath:
            continue
        for fname in filenames:
            # limit to text-like files
            p = Path(dirpath) / fname
            try:
                if p.suffix in {".py", ".md", ".pyi", ".rst", ".txt", ".toml"}:
                    if fix_file(p):
                        modified.append(str(p.relative_to(ROOT)))
            except Exception as e:
                print("Skipping", p, "error", e)
    if modified:
        print("Modified files:")
        for m in modified:
            print(m)
    else:
        print("No changes made")
