#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compose_alpha_pair_vs_random_3metrics.py

Horizontal strip of the three per-metric paired-vs-random alpha figures,
in the order Observed OTUs -> Shannon -> Inverse Simpson
(user request 2026-10-07). Pure Pillow composition — the source PNGs are
unchanged; all three source figures are the same height already.

Inputs : results/new/figures/06_alpha_pair_vs_random/
           alpha_paired_vs_random_filter_observed_otus.png
           alpha_paired_vs_random_filter_shannon.png
           alpha_paired_vs_random_filter_invsimpson.png
Outputs: results/new/figures/06_alpha_pair_vs_random/
           alpha_paired_vs_random_filter_3metrics.{png,pdf}
"""

from pathlib import Path

from PIL import Image

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
DIR = ROOT / "results/new/figures/06_alpha_pair_vs_random"

PANELS = [
    "alpha_paired_vs_random_filter_observed_otus.png",
    "alpha_paired_vs_random_filter_shannon.png",
    "alpha_paired_vs_random_filter_invsimpson.png",
]

MARGIN = 20
GAP = 40


def main() -> None:
    imgs = [Image.open(DIR / p).convert("RGB") for p in PANELS]
    h = min(im.height for im in imgs)
    imgs = [im.resize((round(im.width * h / im.height), h), Image.LANCZOS)
            for im in imgs]

    w = sum(im.width for im in imgs) + GAP * (len(imgs) - 1)
    canvas = Image.new("RGB", (MARGIN * 2 + w, MARGIN * 2 + h), "white")
    x = MARGIN
    for im in imgs:
        canvas.paste(im, (x, MARGIN))
        x += im.width + GAP

    canvas.save(DIR / "alpha_paired_vs_random_filter_3metrics.png", dpi=(300, 300))
    canvas.save(DIR / "alpha_paired_vs_random_filter_3metrics.pdf", resolution=300)
    print(f"canvas {canvas.width}x{canvas.height} -> {DIR}")


if __name__ == "__main__":
    main()