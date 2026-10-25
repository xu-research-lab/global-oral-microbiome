#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_paired_2x2_fontscale.py — uniform-fontscale variant of the paired 2x2
module networks (trial; does NOT touch the standard deliverables).

Same mechanism as make_paired_2x2.py, but sets NET_FONT_SCALE so that
module_network_figures multiplies EVERY fontsize (genus labels, M# badges,
title band, shared legend) by the factor while keeping the canvas
unchanged — text therefore reads proportionally larger when the figure is
scaled to page width. All outputs carry a _fs{S} suffix.

Usage (from network/):
  python script_hc_pd_dc_oscc/make_paired_2x2_fontscale.py [SCALE]
  (default SCALE = 1.25)

Outputs published into result_paired/backbone_v2/paired_figures/:
  paired_module_networks_2x2_{mode}_fs{S}.{png,pdf}
The standard paired_module_networks_2x2_{mode}.*  are untouched.
"""

import os
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

SCALE = sys.argv[1] if len(sys.argv) > 1 else "1.25"
SUFFIX = f"_fs{float(SCALE):g}"
MODES = ("plain", "font")


def main():
    res = os.path.join(ROOT, "result_paired")
    pub = os.path.join(ROOT, "result_paired", "backbone_v2",
                       "paired_figures")

    os.environ["NET_FONT_SCALE"] = SCALE   # must be set BEFORE the import
    import module_network_figures as fig
    assert abs(fig.FS - float(SCALE)) < 1e-9, "NET_FONT_SCALE not picked up"
    fig.ROOT_DIR = res
    fig.BIOM_DIR = "biom_paired"
    sys.argv = ["module_network_figures", "--layouts", "2x2"]
    fig.main()   # writes *_fs{S} names only (SUFFIX mechanism)

    os.makedirs(pub, exist_ok=True)
    src_dir = os.path.join(res, "module_networks")
    n_pub = 0
    for mode in MODES:
        for ext in ("png", "pdf"):
            src = os.path.join(
                src_dir, f"module_networks_2x2_{mode}{SUFFIX}.{ext}")
            if os.path.exists(src):
                dst = os.path.join(
                    pub, f"paired_module_networks_2x2_{mode}{SUFFIX}.{ext}")
                shutil.copy2(src, dst)
                print(f"[pub] {os.path.relpath(dst, ROOT)}")
                n_pub += 1
            else:
                print(f"[warn] missing {src}")
    print(f"\n[done] {n_pub} files -> {pub} "
          "(standard 2x2 and 1x4 outputs untouched)")


if __name__ == "__main__":
    main()