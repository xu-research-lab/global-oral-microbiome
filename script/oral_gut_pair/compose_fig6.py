#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compose_fig6.py  —  rebuild of the 04 composite (replaces the old
make_fig6_composite.py, whose panel D input had been deleted and whose
panels are superseded; that script and its old output were removed).

Layout (all panels pasted as PNG renders, aspect preserved):
    row 1:  a (05 beta paired-vs-random, unweighted UniFrac)
            b (07 alpha Oral–Gut Spearman, raw Observed)
            c (11 pooled PCoA density)
    row 2:  d (region-AUC heatmap)   e (03 oral fraction LME forest), centered

Each row uses the same panel height. Row 2 is centered beneath row 1.
PIL composition (Pillow), with panel letters in each image's top-left corner.

Inputs (all under results/new/figures/):
    05_beta_pair_vs_random/beta_pooled_pair_vs_random_filter_unweighted_unifrac.png
    07_alpha_spearman_by_disease/alpha_oral_gut_spearman_by_disease_raw_observed.png
    11_pooled_pcoa_density/pcoa_pooled_paired_adonis2_unweighted_unifrac_density.png
    Disease_classification/微信图片_20260930155216_405_52.png
    03_oral_fraction/oral_fraction_by_disease_lme.png
Outputs: results/new/figures/04_fig6_composite/fig6_composite.{png,pdf}
"""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path("/home/hexincheng/projects/hxc/oral_project/oral_gut_pair")
FIG = ROOT / "results/new/figures"
OUT_DIR = ROOT / "results/new/figures/04_fig6_composite"

PANELS = {
    "a": FIG / "05_beta_pair_vs_random/beta_pooled_pair_vs_random_filter_unweighted_unifrac.png",
    "b": FIG / "07_alpha_spearman_by_disease/alpha_oral_gut_spearman_by_disease_raw_observed.png",
    "c": FIG / "11_pooled_pcoa_density/pcoa_pooled_paired_adonis2_unweighted_unifrac_density.png",
    "d": FIG / "Disease_classification/微信图片_20260930155216_405_52.png",
    "e": FIG / "03_oral_fraction/oral_fraction_by_disease_lme.png",
}

# user request (2026-10-06): switch to a 2-above / 3-below layout — a, b on
# the top row, c, d, e on the bottom row
ROW1 = ["a", "b"]
ROW2 = ["c", "d", "e"]

PANEL_H = 900
GAP = 40
MARGIN = 20
FONT_SIZE = 70
# Shift letter chips into the inter-panel gap where they would cover titles.
# 'b' shifted further left so its chip stops clipping the "HC (n=2103)"
# legend marker in panel b (user report: b往左移一点点).
LETTER_DX = {"b": -70, "e": -46}

FONT_PATH = Path(
    "/home/hexincheng/miniconda3/envs/py39/lib/python3.9/site-packages"
    "/matplotlib/mpl-data/fonts/ttf/DejaVuSans-Bold.ttf"
)


def fit_height(img: Image.Image, h: int) -> Image.Image:
    w = round(img.width * h / img.height)
    return img.resize((w, h), Image.LANCZOS)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    font = ImageFont.truetype(str(FONT_PATH), FONT_SIZE)

    imgs = {k: Image.open(v).convert("RGB") for k, v in PANELS.items()}

    # --- two panels above, three below, both rows centred ----------------
    # canvas spans the WIDER row (with 2-up/3-down the bottom row is wider,
    # so the old "canvas = row width" rule would clip the bottom panels)
    row1 = [fit_height(imgs[k], PANEL_H) for k in ROW1]
    row2 = [fit_height(imgs[k], PANEL_H) for k in ROW2]
    row1_width = sum(im.width for im in row1) + GAP * (len(row1) - 1)
    row2_width = sum(im.width for im in row2) + GAP * (len(row2) - 1)

    canvas_w = MARGIN * 2 + max(row1_width, row2_width)
    canvas = Image.new(
        "RGB",
        (canvas_w, MARGIN * 2 + PANEL_H * 2 + GAP),
        "white",
    )
    x0, y0 = MARGIN, MARGIN
    y2 = y0 + PANEL_H + GAP
    x = x0 + (canvas_w - row1_width) // 2
    for im in row1:
        canvas.paste(im, (x, y0))
        x += im.width + GAP
    x = x0 + (canvas_w - row2_width) // 2
    for im in row2:
        canvas.paste(im, (x, y2))
        x += im.width + GAP

    dr = ImageDraw.Draw(canvas)
    offs = {}
    x = x0 + (canvas_w - row1_width) // 2
    for k, im in zip(ROW1, row1):
        offs[k] = (x, y0)
        x += im.width + GAP
    x = x0 + (canvas_w - row2_width) // 2
    for k, im in zip(ROW2, row2):
        offs[k] = (x, y2)
        x += im.width + GAP
    for letter, (ox, oy) in offs.items():
        # white background chip so the letter stays readable on any panel
        # ('e' shifted LEFT into the border gap so its chip stops covering
        # panel e's title — user report + request: e往左移一点点)
        # user request (2026-10-08): panel letters in UPPER CASE
        dx = LETTER_DX.get(letter, 0)
        letter = letter.upper()
        tb = dr.textbbox((0, 0), letter, font=font)
        dr.rectangle([ox + dx + 18, oy + 14, ox + dx + 18 + (tb[2] - tb[0]) + 22,
                      oy + 14 + (tb[3] - tb[1]) + 22], fill="white")
        dr.text((ox + dx + 28, oy + 24), letter, font=font, fill="black")

    canvas.save(OUT_DIR / "fig6_composite.png", dpi=(300, 300))
    canvas.save(OUT_DIR / "fig6_composite.pdf", resolution=300)
    print(f"canvas {canvas.width}x{canvas.height} -> {OUT_DIR}")
    print("row1 " + "+".join(str(im.width) for im in row1)
          + f" h={PANEL_H}, row2 " + "+".join(str(im.width) for im in row2)
          + f" h={PANEL_H} (centered)")


if __name__ == "__main__":
    main()