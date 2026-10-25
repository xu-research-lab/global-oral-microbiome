#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_paired_downstream.py — Run the HC/PD/DC/OSCC downstream chain on
result_paired/ WITHOUT touching the original scripts or original results.

The original scripts (backbone_pooling, modularity, extract_module_roles,
combine_module_members_table, jaccard_heatmap) stay byte-identical; this
runner imports them and patches their module-level path constants to the
paired result root before calling their main(), so every output lands under
  result_paired/{backbone_v2, modularity, module_networks, jaccard}/

Chain (matching the >30-threshold pipeline):
  1. backbone_pooling.main()            (skipped if backbone edges exist)
  2. modularity.main()                  Louvain modules + Q/null diagnostics
  3. extract_module_roles.main()        module hubs / main connectors
  4. combine_module_members_table.main()  members table (md + csv)
  5. jaccard_heatmap.main()             module overlap vs HC

Because the original scripts hard-code the ">30-threshold / DC only 5
cohorts / HC 113 cohorts" wording in their md outputs (wrong for the paired
run), the last step
rewrites those known strings to paired-run wording. The original
result_paired/ md files are never re-generated, so they keep the old
wording untouched.

Usage (from network/):
  python script_hc_pd_dc_oscc/run_paired_downstream.py
  python script_hc_pd_dc_oscc/run_paired_downstream.py --skip-pooling
"""

import argparse
import importlib
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SCRIPT_DIR)                      # network/
sys.path.insert(0, SCRIPT_DIR)


def stage(name, modname, patches, argv=()):
    """Import a module, patch its path constants, run its main()."""
    print(f"\n===== {name} =====", flush=True)
    mod = importlib.import_module(modname)
    for k, v in patches.items():
        setattr(mod, k, v)
    sys.argv = [modname, *argv]
    mod.main()


def fixup_paired_text(base):
    """Rewrite original-run wording inside the NEW md outputs under
    result_paired/ (the originals under result_paired/ are untouched)."""
    fixes = [
        # combine_module_members_table.md
        ("result_paired/module_networks/all_module_members_table.md", [
            ("# HC / PD / OSCC / DC networks: consolidated module member table",
             "# HC / PD / OSCC / DC networks (paired, 37 projects): "
             "consolidated module member table"),
            ("**Note**: DC has only 5 cohorts; its network has limited "
             "robustness and the module partition is indicative only.",
             "**Note**: paired scope (HC 35 / PD 18 / DC 9 / OSCC 10 cohorts, "
             "each project containing both disease and healthy samples); "
             "networks with fewer cohorts (especially DC) have limited "
             "partition robustness, values are indicative only."),
        ]),
        # jaccard_heatmap summary
        ("result_paired/jaccard/jaccard_summary.md", [
            ("Module Jaccard similarity, HC / ",
             "Module Jaccard similarity (paired, 37 projects), HC / "),
            ("Healthy network = all-site HC network from this pipeline; "
             "disease networks = ",
             "Healthy network = paired-scope HC network (HC units of the 37 "
             "paired projects); disease networks = "),
            ("- HC 113 cohorts vs PD 11 / DC 5 / OSCC 10 cohorts: the "
             "difference in genus-space coverage systematically depresses\n"
             "  Jaccard (the HC module genus pool is far larger).",
             "- Cohorts: HC 35 / PD 18 / DC 9 / OSCC 10; differences in "
             "genus-space coverage affect the Jaccard scale."),
            ("- DC has only 5 cohorts; its network has limited robustness "
             "and the overlap values are indicative only.",
             "- DC has only 9 cohorts (paired); its network has limited "
             "robustness and the overlap values are indicative only."),
        ]),
    ]
    rel_md = os.path.relpath(base, ROOT)
    for rel, pairs in fixes:
        fp = os.path.join(ROOT, rel)
        if not os.path.exists(fp):
            continue
        with open(fp, encoding="utf-8") as fh:
            txt = fh.read()
        for old, new in pairs:
            txt = txt.replace(old, new)
        with open(fp, "w", encoding="utf-8") as fh:
            fh.write(txt)
        print(f"[text-fix] {rel_md}->{rel} applied")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result-dir", default=os.path.join(ROOT, "result_paired"))
    ap.add_argument("--biom-dir", default="biom_paired",
                    help="per-group genus count subdir under --result-dir")
    ap.add_argument("--skip-pooling", action="store_true",
                    help="backbone pooling already done")
    args = ap.parse_args()

    res = os.path.abspath(args.result_dir)
    biom = os.path.join(res, args.biom_dir)
    backbone = os.path.join(res, "backbone_v2", "HC_genus_uw",
                            "backbone_edges.csv")
    mod_dir = os.path.join(res, "modularity")
    out_net = os.path.join(res, "module_networks")

    # 1. NetMoss pooling (only if not already done — it is deterministic;
    #    backbone_pooling.main() reads units.csv from its RESULT global)
    if os.path.exists(backbone) or args.skip_pooling:
        print(f"[pool] backbone edges present ({backbone}), skip pooling")
    else:
        stage("pooling", "backbone_pooling", dict(
            RESULT=res, SPARCC_BASE=os.path.join(res, "per_site_sparcc"),
            OUT_BASE=os.path.join(res, "backbone_v2"),
            SUM_DIR=os.path.join(res, "backbone_summary")))

    # 2. Louvain modularity
    stage("modularity", "modularity", dict(
        BACKBONE_DIR=os.path.join(res, "backbone_v2"), OUT_DIR=mod_dir))

    # 3. module roles — patch module_network_figures BEFORE importing
    #    extract_module_roles (its `from ... import ROOT_DIR` binds at import)
    print("\n===== module roles =====", flush=True)
    import module_network_figures as fig
    fig.ROOT_DIR = res
    fig.BIOM_DIR = args.biom_dir
    roles = importlib.import_module("extract_module_roles")
    roles.OUT_DIR = out_net
    sys.argv = ["extract_module_roles"]
    roles.main()

    # 4. members table (needs roles)
    stage("members table", "combine_module_members_table", dict(
        ABU_DIR=biom, MOD_DIR=mod_dir, OUT_DIR=out_net,
        ROLES_PATH=os.path.join(out_net, "module_roles.csv")))

    # 5. jaccard vs HC
    stage("jaccard", "jaccard_heatmap", dict(
        MOD_DIR=mod_dir, OUT_DIR=os.path.join(res, "jaccard")))

    fixup_paired_text(res)
    print("\n[done] paired downstream outputs ->", res)


if __name__ == "__main__":
    main()