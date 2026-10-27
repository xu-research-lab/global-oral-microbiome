#!/usr/bin/env python3
"""
pipeline_v2 Step 12g: per-phylotype profile table for ALL phylotypes (26,311).

One row per phylotype, ranked by max_subject_prevalence descending:
  - phylo_n_asvs, phylo_max_gg2_id, phylo_min_gg2_id, phylo_max_homd_id,
    phylo_max_subject_prevalence, phylo_n_bioprojects, phylo_n_continents,
    vsearch_phylum (rep / max-identity member, GG2 2024.09) and
    tree_phylum (MAJORITY phylum of the placement-edge descendant reference tips,
    GG2 2022.10 tree taxonomy — consistent with the tree), species_in_homd,
    contamination flag

Output: 12g_phylotype_profile.tsv
"""
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
TAXONOMY = Path("/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv")
TAXONOMY22 = PIPE_V2 / "gappa" / "gg2_2022.10_tree_taxonomy.tsv"
PLACEMENTS = BASE / "pipeline/05_sepp/placements_export/placements.json"

log = lambda s: print(f"[12g_profile] {s}", flush=True)


def main():
    log("Loading ASV -> phylotype map...")
    mp = pd.read_csv(PIPE_V2 / "10_asv_phylotype_map.tsv", sep="\t", low_memory=False,
                     usecols=["derep_id", "phylotype_id", "edge_num"])

    log("Loading per-ASV GG2 annotation...")
    ann = pd.read_csv(PIPE_V2 / "10_gg2_closedref_annotation.tsv", sep="\t",
                      low_memory=False,
                      usecols=["derep_id", "best_identity", "best_target"])

    log("Loading per-ASV HOMD identity (step 01)...")
    h01 = pd.read_csv(PIPE_V2 / "01_chimera_filtered.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "identity"])
    h01 = h01.rename(columns={"identity": "homd_identity"})

    log("Loading phylotype-level table (step 11)...")
    phy = pd.read_csv(PIPE_V2 / "11_phylotypes_homd_novelty.tsv", sep="\t",
                      low_memory=False)
    phy = phy[["phylotype_id", "n_asvs", "max_subject_prevalence",
               "n_bioprojects_total", "n_continents", "clean_phylotype",
               "novelty_class", "species_in_homd", "rep_derep_id"]]

    # ---- per-phylotype aggregates ----
    m = mp.merge(ann, on="derep_id", how="left").merge(h01, on="derep_id", how="left")
    m["best_identity"] = pd.to_numeric(m["best_identity"], errors="coerce")
    m["homd_identity"] = pd.to_numeric(m["homd_identity"], errors="coerce")

    agg = m.groupby("phylotype_id").agg(
        phylo_max_gg2_id=("best_identity", "max"),
        phylo_min_gg2_id=("best_identity", "min"),
        phylo_max_homd_id=("homd_identity", "max"),
    ).reset_index()
    phy = phy.merge(agg, on="phylotype_id", how="left")

    # ---- phylum of the rep ASV (fallback: max-identity member) ----
    log("Loading GG2 2024.09 taxonomy (phylum lookup)...")
    tax_phylum = {}
    with open(TAXONOMY) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t", 1)
            if len(parts) < 2:
                continue
            phylum = ""
            for tok in parts[1].split(";"):
                tok = tok.strip()
                if tok.startswith("p__"):
                    phylum = tok[3:]
                    break
            tax_phylum[parts[0]] = phylum

    def phylum_of(target):
        if not isinstance(target, str) or target not in tax_phylum:
            return ""
        return tax_phylum[target]

    m["phylum"] = m["best_target"].map(phylum_of)
    # rep phylum (by rep_derep_id), else max-identity member's phylum
    rep_phylum = m.dropna(subset=["best_identity"]).sort_values(
        "best_identity", ascending=False).drop_duplicates("phylotype_id")
    rep_phylum = rep_phylum.set_index("phylotype_id")["phylum"]
    phy["phylum"] = phy["phylotype_id"].map(rep_phylum).fillna("")
    phy = phy.rename(columns={"phylum": "vsearch_phylum"})

    # ---- TREE-BASED phylum (placement edge descendant-tip majority, 2022.10) ----
    log("Loading 2022.10 placement tree + tree taxonomy (tree_phylum)...")
    sys.setrecursionlimit(1000000)
    with open(PLACEMENTS) as f:
        s = json.load(f)["tree"]
    par, kids, edge, tip = [], [], [], []

    def mk(p=None, e=None, t=None):
        idx = len(par); par.append(p); kids.append([]); edge.append(e); tip.append(t)
        if p is not None:
            kids[p].append(idx)
        return idx

    i = 0
    root = mk()
    stack = [root]
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append(mk(stack[-1])); i += 1
        elif c == ")":
            closed = stack.pop(); i += 1
            while i < len(s) and s[i] not in ":[],();":
                i += 1
            if i < len(s) and s[i] == ":":
                i += 1
                while i < len(s) and s[i] not in "[],();":
                    i += 1
            if i < len(s) and s[i] == "[":
                k = s.index("]", i); edge[closed] = int(s[i + 1:k]); i = k + 1
        elif c == ",":
            i += 1
        elif c == ";":
            break
        else:
            j = i
            while j < len(s) and s[j] not in ":,()[];":
                j += 1
            nd = mk(stack[-1], t=s[i:j]); i = j
            if i < len(s) and s[i] == ":":
                i += 1
                while i < len(s) and s[i] not in "[],();":
                    i += 1
            if i < len(s) and s[i] == "[":
                k = s.index("]", i); edge[nd] = int(s[i + 1:k]); i = k + 1
    edge_to_node = {e: idx for idx, e in enumerate(edge) if e is not None}

    def desc_tips(n):
        if tip[n] is not None:
            return [tip[n]]
        out = []
        st = list(kids[n])
        while st:
            k = st.pop()
            if tip[k] is not None:
                out.append(tip[k])
            else:
                st.extend(kids[k])
        return out

    tax22 = {}
    with open(TAXONOMY22) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            ph = ""
            for tok in parts[1].split(";"):
                tok = tok.strip()
                if tok.startswith("p__"):
                    ph = tok[3:]
                    break
            tax22[parts[0]] = ph

    def resolve_phylum(node_idx):
        """Majority phylum of the descendant clade; walk UP to the nearest
        ancestor whose clade has annotated tips (deepest resolvable phylum)."""
        n = node_idx
        seen = set()
        while n is not None and n not in seen:
            seen.add(n)
            dts = desc_tips(n)
            phys = [tax22.get(t, "") for t in dts if tax22.get(t, "")]
            if phys:
                return Counter(phys).most_common(1)[0][0]
            n = par[n]
        return ""

    tree_phylum = {}
    edge_map = mp.drop_duplicates("phylotype_id")[["phylotype_id", "edge_num"]]
    for _, r in edge_map.iterrows():
        e = r["edge_num"]
        if pd.isna(e) or int(e) not in edge_to_node:
            tree_phylum[r["phylotype_id"]] = ""
            continue
        tree_phylum[r["phylotype_id"]] = resolve_phylum(edge_to_node[int(e)])
    phy["tree_phylum"] = phy["phylotype_id"].map(tree_phylum).fillna("")
    log(f"  tree_phylum assigned for {int((phy['tree_phylum'] != '').sum()):,} phylotypes")

    phy = phy.sort_values("max_subject_prevalence", ascending=False).reset_index(drop=True)
    phy["rank"] = phy.index + 1

    # categorical tracks (already in phy): clean, novelty, species_in_homd
    for c in ["clean_phylotype", "species_in_homd"]:
        phy[c] = phy[c].astype(str).str.lower().isin(["true", "1"])

    out = PIPE_V2 / "12g_phylotype_profile.tsv"
    phy.to_csv(out, sep="\t", index=False)
    log(f"saved {out} ({len(phy):,} rows)")
    log(f"  vsearch phyla: {phy['vsearch_phylum'].nunique()} | tree phyla: "
        f"{phy['tree_phylum'].nunique()} | tree phylum missing: "
        f"{int((phy['tree_phylum'] == '').sum()):,}")
    log("Done!")


if __name__ == "__main__":
    main()
