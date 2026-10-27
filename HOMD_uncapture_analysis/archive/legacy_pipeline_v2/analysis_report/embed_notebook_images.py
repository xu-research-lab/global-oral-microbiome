#!/usr/bin/env python3
"""
Embed notebook markdown images (figs_v2/*.png) as base64 data URIs so the
notebook is SELF-CONTAINED and renders on GitHub / nbviewer without uploading
the figs_v2 directory.

Usage: python embed_notebook_images.py [src.ipynb] [dst.ipynb]
"""
import base64
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "novel_oral_taxa_report_v2.ipynb"
DST = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "novel_oral_taxa_report_v2_embedded.ipynb"

IMG_RE = re.compile(r"(!\[[^\]]*\]\()([^)]+)(\))")
n_embedded = 0


def to_data_uri(path: str) -> str:
    p = (HERE / path).resolve()
    ext = p.suffix.lower()
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg",
            "pdf": "application/pdf", "svg": "image/svg+xml"}.get(ext.lstrip("."), "image/png")
    b64 = base64.b64encode(p.read_bytes()).decode()
    return f"data:{mime};base64,{b64}"


def main():
    nb = json.loads(SRC.read_text())
    n_embedded = 0
    for cell in nb["cells"]:
        if cell["cell_type"] != "markdown":
            continue
        src = "".join(cell["source"])

        def repl(m):
            global n_embedded
            path = m.group(2)
            if path.startswith("data:"):
                return m.group(0)
            n_embedded += 1
            return f"{m.group(1)}{to_data_uri(path)}{m.group(3)}"

        cell["source"] = [IMG_RE.sub(repl, src)]
    DST.write_text(json.dumps(nb))
    size_mb = DST.stat().st_size / 1e6
    print(f"embedded {n_embedded} images -> {DST} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
