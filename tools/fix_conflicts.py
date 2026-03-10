#!/usr/bin/env python3
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def fix_file(path: Path) -> bool:
    text = path.read_text(encoding="utf8")
    if "<<<<<<< HEAD" not in text:
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
        j = text.find("=======", idx)
        if j == -1:
            out.append(text[idx:])
            break
        k = text.find(">>>>>>>", j)
        if k == -1:
            out.append(text[idx:])
            break
        # find end of the line that contains the >>>>>>> marker to skip commit id
        line_end = text.find("\n", k)
        if line_end == -1:
            line_end = k + len(">>>>>>>")
        # keep HEAD part between idx+len(marker) and j
        head_part = text[idx + len("<<<<<<< HEAD"):j]
        out.append(head_part)
        i = line_end + 1
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
