#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_paired_2x2.py — 2x2 arrangement of the PAIRED module networks.

Same mechanism as run_paired_figures.py: imports module_network_figures,
patches its ROOT_DIR / BIOM_DIR to the paired results
(result_paired/{backbone_v2, modularity, biom_paired}) and runs it with
--layouts 2x2, so ONLY the new 2x2 combined figures are generated
(the original 1x4 deliverables are NOT regenerated or overwritten).

Outputs published into result_paired/backbone_v2/paired_figures/:
  paired_module_networks_2x2_{plain,font}.{png,pdf}

Usage (from network/):
  python script_hc_pd_dc_oscc/make_paired_2x2.py
"""

import os
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

MODES = ("plain", "font")


def main():
    res = os.path.join(ROOT, "result_paired")
    pub = os.path.join(ROOT, "result_paired", "backbone_v2",
                       "paired_figures")

    import module_network_figures as fig
    fig.ROT_DEG = 0.0    # pin: the approved 2x2 layout stays UNROTATED
                         # (module default ROT 45° applies to the 1x4 only)
    fig.ROOT_DIR = res
    fig.BIOM_DIR = "biom_paired"
    sys.argv = ["module_network_figures", "--layouts", "2x2"]
    fig.main()                       # writes module_networks_2x2_* only

    os.makedirs(pub, exist_ok=True)
    src_dir = os.path.join(res, "module_networks")
    n_pub = 0
    for mode in MODES:
        for ext in ("png", "pdf"):
            src = os.path.join(src_dir, f"module_networks_2x2_{mode}.{ext}")
            if os.path.exists(src):
                dst = os.path.join(
                    pub, f"paired_module_networks_2x2_{mode}.{ext}")
                shutil.copy2(src, dst)
                print(f"[pub] {os.path.relpath(dst, ROOT)}")
                n_pub += 1
            else:
                print(f"[warn] missing {src}")
    print(f"\n[done] {n_pub} files -> {pub} (1x4 outputs untouched)")


if __name__ == "__main__":
    main()