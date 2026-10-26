"""Shared network data, abundance and cartography roles."""
from pathlib import Path
import glob
import os

import numpy as np
import pandas as pd
import networkx as nx

MODULE_DIR = Path(__file__).resolve().parents[1]
RESULT_DIR = Path(os.environ.get("NETWORK_RESULT_DIR", MODULE_DIR / "result_paired"))
GROUPS = ["HC", "PD", "DC", "OSCC"]
GRP_LABEL = {"HC": "Healthy control", "PD": "Periodontal disease", "DC": "Dental caries", "OSCC": "OSCC"}
BACKBONE_DIR, MODULARITY_DIR, BIOM_DIR = "backbone_v2", "modularity", "biom_paired"
ROLE_RED, ROLE_BLUE, ROLE_PURPLE = "#ff8c00", "#2ca02c", "#7b2fa0"

def load_mean_abundance(grp, root):
    """Mean per-sample relative abundance of every genus across ALL units of
    the group (the pooled size basis for this group's figure)."""
    sums, n_samp = {}, 0
    pattern = os.path.join(root, BIOM_DIR, f"{grp}_genus", "*.tsv")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"Missing cohort count tables: {pattern}")
    for fn in files:
        df = pd.read_csv(fn, sep="\t", index_col=0)
        obs = df.to_numpy(dtype=float)         # genus x sample (counts)
        s = obs.sum(axis=0)
        s[s == 0] = np.nan
        rel = obs / s                          # per-sample relative abundance
        for i, tax in enumerate(df.index):
            v = rel[i]
            v = v[np.isfinite(v)]
            if len(v):
                sums[str(tax)] = sums.get(str(tax), 0.0) + v.sum()
        n_samp += obs.shape[1]
    if not sums or not n_samp:
        return {}
    return {g: v / n_samp for g, v in sums.items()}

def assemble(grp, root):
    """Graph + module assignment for one group network."""
    d = os.path.join(root, BACKBONE_DIR, f"{grp}_genus_uw")
    edges = pd.read_csv(os.path.join(d, "backbone_edges.csv"))
    mods = pd.read_csv(os.path.join(root, MODULARITY_DIR, f"modules_{grp}.csv"),
                       dtype={"taxon_id": str})
    mod_of = dict(zip(mods.taxon_id, mods.module))
    edge_nodes = set(edges.taxon_i.astype(str)) | set(edges.taxon_j.astype(str))
    if mods.taxon_id.duplicated().any() or edge_nodes - set(mod_of):
        raise ValueError(f"Invalid or incomplete module assignments for {grp}")
    abu = load_mean_abundance(grp, root)

    G = nx.Graph()
    for _, e in edges.iterrows():
        G.add_edge(str(e.taxon_i), str(e.taxon_j),
                   r=float(e.pooled_obs), w=abs(float(e.pooled_obs)))
    for n in G.nodes():
        G.nodes[n]["abu"] = abu.get(n, np.nan)
        G.nodes[n]["mod"] = mod_of.get(n, 0)
        G.nodes[n]["rank"] = G.nodes[n]["mod"]

    hubs = {}
    for n, d0 in G.nodes(data=True):
        m = d0["rank"]
        if m not in hubs or G.degree(n) > G.degree(hubs[m]):
            hubs[m] = n
    G.graph["hubs"] = hubs
    G.graph["site"] = grp
    return G

def cartography_df(G):
    """z_in (within-module degree z-score, sd=0 modules -> 0) and the
    participation coefficient P = 1 - sum_s (k_is / k_i)^2."""
    rows = []
    for m in sorted({d0["rank"] for _, d0 in G.nodes(data=True)}):
        members = [n for n in G.nodes() if G.nodes[n]["rank"] == m]
        kin = {n: sum(1 for v in G[n] if G.nodes[v]["rank"] == m)
               for n in members}
        arr = np.array([kin[n] for n in members], float)
        sd = float(arr.std())
        for n in members:
            z = (kin[n] - arr.mean()) / sd if sd > 1e-12 else 0.0
            kout = G.degree(n) - kin[n]
            rows.append(dict(taxon=n, module=m, k_in=kin[n], z_in=z,
                             k_out=kout, degree=kin[n] + kout))
    df = pd.DataFrame(rows)

    def _P(name, degree):
        if degree == 0:
            return 0.0
        shares = {}
        for v in G[name]:
            s = G.nodes[v]["rank"]
            shares[s] = shares.get(s, 0) + 1
        return 1.0 - sum((c / degree) ** 2 for c in shares.values())

    df["P"] = [_P(n, d) for n, d in zip(df.taxon, df.degree)]

    def _role(r):
        if r.z_in >= 2.5:
            if r.P < 0.3:
                return "provincial hub"
            return "connector hub" if r.P < 0.75 else "kinless hub"
        if r.P < 0.05:
            return "ultra-peripheral"
        return "peripheral" if r.P < 0.62 else "satellite connector"
    df["role"] = df.apply(_role, axis=1)
    return df

def cartography_role_nodes(G):
    """node -> outline/label colour. module hub = max within-module z (only when
    k_in varies, yellow), main connector = max P (only when cross-module links
    exist, green); hub and connector on the same node -> purple."""
    dfc = cartography_df(G)
    role_nodes = {}
    for m, sub in dfc.groupby("module"):
        sd = float(sub.k_in.std(ddof=0))
        hub = (sub.sort_values(["z_in", "degree"]).iloc[-1]
               if sd > 1e-12 else None)
        con = sub.sort_values(["P", "z_in"]).iloc[-1]
        con = con if con.P > 1e-12 else None
        if hub is not None and con is not None and \
                hub["taxon"] == con["taxon"]:
            role_nodes[hub["taxon"]] = ROLE_PURPLE
        else:
            if hub is not None:
                role_nodes[hub["taxon"]] = ROLE_RED
            if con is not None:
                role_nodes[con["taxon"]] = ROLE_BLUE
    return role_nodes

def prepare_graphs(root=RESULT_DIR):
    graphs = {group: assemble(group, root) for group in GROUPS}
    for graph in graphs.values():
        graph.graph["role_nodes"] = cartography_role_nodes(graph)
    return graphs
