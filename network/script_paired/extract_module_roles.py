#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_module_roles.py — Export the module hubs and main connectors shown in
module_networks_1x4.png (the cartography role-marked nodes) to a CSV table.

Roles come from the SAME cartography rule as the figure
(module_network_figures.py, Guimerà-Amaral):
  module hub     = max within-module degree z-score of the module
                   (only assigned when k_in varies within the module)
  main connector = max participation coefficient P of the module
                   (only assigned when cross-module links exist, P > 0)
  a node that is both is listed in BOTH columns (drawn purple in the figure).

Columns: module, health_status, module_hub, main_connector
(empty cells when a module has no such node). Order: HC, PD, OSCC, DC.

Inputs:
  result_paired/backbone_v2/<grp>_genus_uw/backbone_edges.csv
  result_paired/modularity/modules_<grp>.csv
  result_paired/biom_paired/<grp>_genus/*.tsv
Output:
  result_paired/module_networks/module_roles.csv

Usage (from network/):
  python script_hc_pd_dc_oscc/extract_module_roles.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from module_network_figures import (  # noqa: E402
    ROOT_DIR, MODULARITY_DIR, BACKBONE_DIR, GRP_LABEL, ROLE_BLUE,
    ROLE_PURPLE, assemble, cartography_role_nodes)

import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, "result_paired", "module_networks")

# content order: HC, PD, OSCC, DC (as in the module members table)
ORDER = ["HC", "PD", "OSCC", "DC"]


def main():
    rows = []
    for grp in ORDER:
        mods_path = os.path.join(ROOT, ROOT_DIR, MODULARITY_DIR,
                                 f"modules_{grp}.csv")
        edges_path = os.path.join(ROOT, ROOT_DIR, BACKBONE_DIR,
                                  f"{grp}_genus_uw", "backbone_edges.csv")
        if not (os.path.exists(mods_path) and os.path.exists(edges_path)):
            print(f"[skip] {grp}: missing backbone/modules output")
            continue
        G = assemble(grp, os.path.join(ROOT, ROOT_DIR))
        roles = cartography_role_nodes(G)
        all_mods = sorted({G.nodes[n]["rank"] for n in G.nodes()})

        # group role nodes by module, keeping the figure's colour semantics
        hub_of, con_of = {}, {}
        for n, colour in roles.items():
            m = G.nodes[n]["rank"]
            if colour == ROLE_BLUE:          # connector-only node
                con_of[m] = n
            elif colour == ROLE_PURPLE:      # hub AND connector on one node
                hub_of[m] = n
                con_of[m] = n
            else:                            # hub-only node (yellow)
                hub_of[m] = n

        for m in all_mods:      # every module gets a row; blank when no role
            rows.append(dict(
                module=f"{grp}-M{m}",
                health_status=GRP_LABEL.get(grp, grp),
                module_hub=hub_of.get(m, ""),
                main_connector=con_of.get(m, "")))
        print(f"[roles] {grp}: {len(hub_of)} hubs, {len(con_of)} connectors")

    df = pd.DataFrame(rows, columns=["module", "health_status",
                                     "module_hub", "main_connector"])
    out = os.path.join(OUT_DIR, "module_roles.csv")
    df.to_csv(out, index=False)
    print(f"[done] {len(df)} rows -> {out}")


if __name__ == "__main__":
    main()