#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
combine_module_members_table.py — Merge the module member tables of the four
single-condition all-site networks (HC / PD / OSCC / DC, >30-threshold
pipeline) into one long-format table, mirroring
result_disease/modularity/all_module_members_table.md:

  module | health_status | member (degree, mean rel. abund.)

Each module is one row; the member cell is a compact `genus(degree,abund%)` list
sorted by degree ↓ then mean relative abundance ↓. Abundance is the per-sample
mean relative abundance across ALL cohorts of that network (six decimals),
computed exactly like plot_stable_core_abundance.py.

Content order in the table (per user request): HC, PD, OSCC, DC.

Inputs:
  result_paired/modularity/modules_<grp>.csv   (taxon_id, module, degree, ...)
  result_paired/biom_paired/<grp>_genus/*.tsv      (genus counts)
  result_paired/module_networks/module_roles.csv (module_hub, main_connector)
Outputs (under result_paired/module_networks/):
  all_module_members_table.md
  all_module_members_table.csv

Usage (from network/):
  python script_hc_pd_dc_oscc/combine_module_members_table.py
"""

import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ABU_DIR = os.path.join(ROOT, "result_paired", "biom_paired")
MOD_DIR = os.path.join(ROOT, "result_paired", "modularity")
OUT_DIR = os.path.join(ROOT, "result_paired", "module_networks")
ROLES_PATH = os.path.join(OUT_DIR, "module_roles.csv")

# (group, status label, module-id prefix) — content order HC, PD, OSCC, DC
GROUPS = [
    ("HC", "healthy", "HC"),
    ("PD", "PD", "PD"),
    ("OSCC", "OSCC", "OSCC"),
    ("DC", "DC", "DC"),
]


def strip_g(t):
    return str(t)[3:] if str(t).startswith("g__") else str(t)


def load_abundance(grp):
    """Mean per-sample rel. abundance per genus across all cohorts of grp."""
    d = os.path.join(ABU_DIR, f"{grp}_genus")
    sums, n_samp = {}, 0
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".tsv"):
            continue
        df = pd.read_csv(os.path.join(d, fn), sep="\t", index_col=0)
        obs = df.to_numpy(dtype=np.float64)
        s = obs.sum(axis=0)
        s[s == 0] = np.nan
        rel = obs / s
        for i, tax in enumerate(df.index):
            key = strip_g(tax)
            v = rel[i]
            v = v[np.isfinite(v)]
            if len(v):
                sums[key] = sums.get(key, 0.0) + v.sum()
        n_samp += df.shape[1]
    return {g: v / n_samp for g, v in sums.items()}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = []
    for grp, status, prefix in GROUPS:
        abu = load_abundance(grp)
        mod = pd.read_csv(os.path.join(MOD_DIR, f"modules_{grp}.csv"),
                          dtype={"taxon_id": str})
        for m in sorted(mod.module.unique()):
            sub = mod.loc[mod.module == m]
            for _, r in sub.iterrows():
                genus = strip_g(r.taxon_id)
                rows.append(dict(
                    module=f"{prefix}-M{int(m)}",
                    health_status=status,
                    member=genus,
                    degree=int(r.degree),
                    mean_rel_abund=abu.get(genus, 0.0) * 100.0))

    # module roles (hub / main connector) — merged per module, blank if empty
    roles = pd.read_csv(ROLES_PATH, dtype=str).fillna("")
    roles = roles.rename(columns={"module_hub": "module hub",
                                  "main_connector": "main connector"})

    df = pd.DataFrame(rows)
    # order: HC, PD, OSCC, DC; modules numeric; members degree ↓ then abund ↓
    cat_order = {status: i for i, (_, status, _) in enumerate(GROUPS)}
    df["_cat"] = df.health_status.map(cat_order)
    df["_mid"] = df.module.str.split("-").str[-1].str[1:].astype(int)
    df = df.sort_values(["_cat", "_mid", "degree", "mean_rel_abund"],
                        ascending=[True, True, False, False])
    df = df.drop(columns=["_cat", "_mid"])

    # attach hub / main connector (one value per module)
    df = df.merge(roles[["module", "module hub", "main connector"]],
                  on="module", how="left")
    df["module hub"] = df["module hub"].fillna("")
    df["main connector"] = df["main connector"].fillna("")

    df_out = df.copy()
    df_out["mean_rel_abund"] = df_out["mean_rel_abund"].map(
        lambda v: f"{v:.6f}%")
    df_out.to_csv(os.path.join(OUT_DIR, "all_module_members_table.csv"),
                  index=False)

    md = ["# HC / PD / OSCC / DC networks: consolidated module member table",
          "",
          "module numbering: HC = all-site healthy network, PD = periodontal "
          "disease, OSCC = oral squamous cell carcinoma, DC = dental caries",
          "(all-site single-condition networks from this pipeline); one row "
          "per module, the member cell is a compact `genus(degree,mean rel. "
          "abund.)` list ordered by degree then abundance, descending;",
          "abundance is the genus-level mean relative abundance across ALL "
          "cohorts of that network (six decimals).",
          "",
          "**Note**: DC has only 5 cohorts; its network has limited "
          "robustness and the module partition is indicative only.",
          "",
          "| module | health_status | module hub | main connector | member (degree, mean rel. abund.) |",
          "|---|---|---|---|---|"]
    for (module, status), sub in df.groupby(["module", "health_status"],
                                            sort=False):
        cell = "; ".join(f"{r.member}({r.degree},{r.mean_rel_abund:.6f}%)"
                         for _, r in sub.iterrows())
        hub = sub["module hub"].iloc[0]
        conn = sub["main connector"].iloc[0]
        md.append(f"| {module} | {status} | {hub} | {conn} | {cell} |")

    md_out = os.path.join(OUT_DIR, "all_module_members_table.md")
    with open(md_out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    print(f"[done] {df.module.nunique()} modules / {len(df)} members "
          f"-> {md_out}; long csv -> all_module_members_table.csv")


if __name__ == "__main__":
    main()
