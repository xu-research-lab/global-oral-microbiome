#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_paired_figures.py — Draw the paired-network figures and collect them
under result_paired/backbone_v2/paired_figures/ (a NEW subdir; the
original figures of the >30-threshold pipeline stay untouched).

Original scripts are not modified:
  - module_network_figures.main() is imported and its ROOT_DIR / BIOM_DIR
    constants are patched to the paired result root, so it reads
    result_paired/{backbone_v2, modularity, biom_paired} and writes to
    result_paired/module_networks/ first.
  - plot_jaccard_facets_paired.R (a one-line OUT_DIR variant of the
    original plot_jaccard_facets.R, also kept unmodified) draws the facets
    from result_paired/jaccard/.

Everything is then COPIED (never moved/deleted) to
  result_paired/backbone_v2/paired_figures/
  module_networks_1x4.png|.pdf      the 1x4 combined deliverable
  module_network_<grp>_genus.png    standalone panels
  jaccard_facets.png|.pdf           module-overlap facets
  <grp>_backbone.png etc.           the pooling QC plots (copied + renamed)

Usage (from network/):
  python script_hc_pd_dc_oscc/run_paired_figures.py
"""

import argparse
import os
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)                      # network/
sys.path.insert(0, SCRIPT_DIR)

RSERIALIZE = os.environ.get("RSCRIPT",
                     "/home/zhangyuchen/.conda/envs/r_latest/bin/Rscript")

# (src rel to ROOT, dst name) figure files to publish
# the generator (module_network_figures.py) now emits two emphasis versions
# of every module-network figure: "_font" (key genera = bold enlarged label,
# no ring) and "_plain" (no magenta anywhere). The un-suffixed mixed version
# was retired; keep the suffix list in sync with module_network_figures.
_HL_MODES = ("font", "plain")
MODULE_FIGS = (
    [f"module_networks_1x4_{m}.{ext}" for m in _HL_MODES
     for ext in ("png", "pdf")]
    + [f"module_network_{g}_genus_{m}.png" for g in ("HC", "PD", "DC", "OSCC")
       for m in _HL_MODES])
POOLING_FIGS = ["backbone.png", "reproducibility.png",
                "pooled_corr_distribution.png"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result-dir", default=os.path.join(ROOT, "result_paired"))
    ap.add_argument("--biom-dir", default="biom_paired")
    ap.add_argument("--publish-dir",
                    default=os.path.join(ROOT, "result_paired",
                                         "backbone_v2", "paired_figures"))
    args = ap.parse_args()

    res = os.path.abspath(args.result_dir)
    pub = os.path.abspath(args.publish_dir)

    # ── 1x4 module networks (reads result_paired via patched constants) ──
    import module_network_figures as fig
    fig.ROOT_DIR = res
    fig.BIOM_DIR = args.biom_dir
    sys.argv = ["module_network_figures"]
    fig.main()      # writes into <res>/module_networks/

    # ── jaccard facets (unmodified-style copy with paired OUT_DIR) ──
    rs_rscript = RSERIALIZE
    os.system(
        f"{rs_rscript} "
        f"{SCRIPT_DIR}/plot_jaccard_facets_paired.R")

    # ── collect (copy) everything into the publish dir ──
    os.makedirs(pub, exist_ok=True)
    src_net = os.path.join(res, "module_networks")
    for fn in MODULE_FIGS:
        src = os.path.join(src_net, fn)
        if os.path.exists(src):
            dst = os.path.join(pub, fn.replace("module_networks_1x4",
                                               "paired_module_networks_1x4"))
            shutil.copy2(src, dst)
            print(f"[pub] {fn} -> {os.path.relpath(dst, ROOT)}")
        else:
            print(f"[warn] missing {src}")

    src_jac = os.path.join(res, "jaccard")
    for ext in ("png", "pdf"):
        src = os.path.join(src_jac, f"jaccard_facets_paired.{ext}")
        if not os.path.exists(src):            # R-script default name fallback
            src = os.path.join(src_jac, f"jaccard_facets.{ext}")
        if os.path.exists(src):
            dst = os.path.join(pub, f"paired_jaccard_facets.{ext}")
            shutil.copy2(src, dst)
            print(f"[pub] jaccard_facets.{ext} -> {os.path.relpath(dst, ROOT)}")

    for grp in ("HC", "PD", "DC", "OSCC"):
        d = os.path.join(res, "backbone_v2", f"{grp}_genus_uw")
        for fn in POOLING_FIGS:
            src = os.path.join(d, fn)
            if os.path.exists(src):
                dst = os.path.join(pub, f"paired_{grp}_{fn}")
                shutil.copy2(src, dst)
                print(f"[pub] {grp}/{fn} -> {os.path.relpath(dst, ROOT)}")

    print(f"\n[done] figures -> {pub}")


if __name__ == "__main__":
    main()