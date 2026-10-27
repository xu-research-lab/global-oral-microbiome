#!/usr/bin/env python3
"""Supplementary table: the phylotypes called differential in head and neck
cancer and in periodontal disease.

The phylotypes and their statistics come from MaAsLin3
(17_maaslin3_significant.tsv, restricted to the two disease cohorts), the
taxonomy and support from 14_supported_cluster_annotation.tsv, and the member
ASVs, their sequences and their distance to the annotated node from
07_master_asv_table.tsv.  The representative ASV of a multi-ASV phylotype is
the one the pipeline used to aggregate it (rep_derep_id,
08_supported_clusters.tsv).

Two tables are written:

  features  one row per cohort x phylotype -- what was tested, how it is
            annotated, how far its representative ASV sits from that annotation,
            the effect and the q value.
  members   one row per member ASV -- every sequence in a 40- or 112-member
            phylotype is listed here, not hidden behind the representative.

A phylotype whose lineage stops above the species rank carries no species name;
its taxon and rank in the table are the deepest named rank of the lineage.

Run from pipeline_v3.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.styles import Font

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/"
            "HOMD_uncapture_analysis/pipeline_v3")
OUTDIR = BASE / "analysis_report" / "tables"

COHORTS = ["Head and Neck Cancer", "Periodontal Disease"]
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
BASIS = {
    "cluster_refs": "reference sequences within the phylotype",
    "cluster_ancestor": "ancestor node of the phylotype",
    "singleton_ancestor": "nearest reference-bearing ancestor",
}

# ---- inputs -----------------------------------------------------------------
sig = pd.read_csv(BASE / "17_maaslin3_significant.tsv", sep="\t")
sig = sig[sig["cohort"].isin(COHORTS)].copy()

ann = pd.read_csv(BASE / "14_supported_cluster_annotation.tsv", sep="\t",
                  dtype={"cluster_id": str})
sup = pd.read_csv(BASE / "08_supported_clusters.tsv", sep="\t",
                  dtype={"cluster_id": str, "rep_derep_id": str})
mt = pd.read_csv(BASE / "07_master_asv_table.tsv", sep="\t",
                 dtype={"asv_id": str, "cluster_id": str})

features = sorted(sig["feature"].unique())
missing = set(features) - set(ann["cluster_id"]) - set(mt["cluster_id"])
if missing:
    raise SystemExit(f"features missing from the tables: {sorted(missing)}")

ann = ann.set_index("cluster_id")
sup = sup.set_index("cluster_id")
mem = mt[mt["cluster_id"].isin(features)].copy()

# ---- taxonomy ---------------------------------------------------------------
def split_lineage(s):
    """(taxonomy string, deepest named taxon, its rank) from a d__..s__ lineage."""
    parts = [p.strip() for p in str(s).split(";")] if isinstance(s, str) else []
    parts = (parts + [""] * len(RANKS))[:len(RANKS)]
    named = [i for i, p in enumerate(parts) if p]
    if not named:
        return "", "", ""
    deepest = named[-1]
    # an empty rank between two named ones has no name in the reference tree
    shown = [p if p else "unclassified" for p in parts[:deepest + 1]]
    return "; ".join(shown), parts[deepest], RANKS[deepest]


# ---- distance of every member to the node its annotation comes from ----------
mem["annotation_distance"] = pd.to_numeric(mem["annotation_distance"], errors="coerce")
dist = mem.groupby("cluster_id")["annotation_distance"].agg(
    dist_min="min", dist_median="median", dist_max="max")

# ---- representatives --------------------------------------------------------
rep_of = {}
for f in features:
    rep = sup.at[f, "rep_derep_id"] if f in sup.index else None
    ids = mem.loc[mem["cluster_id"] == f, "asv_id"]
    if rep not in set(ids):
        rep = ids.iloc[0]
    rep_of[f] = rep

rep_rows = mem[mem.apply(lambda r: rep_of[r["cluster_id"]] == r["asv_id"], axis=1)]
rep_rows = rep_rows.drop_duplicates("cluster_id").set_index("cluster_id")
if len(rep_rows) != len(features):
    raise SystemExit("could not find one representative row per feature")

# ---- feature table ----------------------------------------------------------
rows = []
for r in sig.itertuples(index=False):
    f = r.feature
    lin, taxon, rank = split_lineage(ann.at[f, "annotation_lineage"])
    rep = rep_rows.loc[f]
    rows.append({
        "cohort": r.cohort,
        "phylotype": f,
        "phylotype_type": "multi-ASV" if f.startswith("c") else "single ASV",
        "n_members": int(ann.at[f, "n_asvs"]),
        "taxonomy": lin,
        "assigned_taxon": taxon,
        "taxon_rank": rank,
        "phylum": ann.at[f, "phylum"] if f in ann.index else "",
        "annotation_source": BASIS[ann.at[f, "annotation_basis"]],
        "annotation_distance_rep": float(rep["annotation_distance"]),
        "annotation_distance_min": float(dist.at[f, "dist_min"]),
        "annotation_distance_median": float(dist.at[f, "dist_median"]),
        "annotation_distance_max": float(dist.at[f, "dist_max"]),
        "homd_identity_pct": float(rep["homd_identity"]),
        "nearest_homd_ref": rep["nearest_homd_ref"],
        "nearest_homd_ref_species": rep["nearest_homd_ref_species"],
        "nearest_homd_ref_distance": float(rep["nearest_homd_ref_dist"]),
        "max_n_subjects": int(ann.at[f, "max_subject_prevalence"]),
        "n_bioprojects": int(ann.at[f, "n_bioprojects_total"]),
        "n_continents": int(ann.at[f, "n_continents"]),
        "total_reads": int(ann.at[f, "total_reads_sum"]),
        "representative_asv": rep["asv_id"],
        "representative_sequence": rep["sequence"],
        "representative_sequence_length": len(rep["sequence"]),
        "member_asvs": ";".join(sorted(mem.loc[mem["cluster_id"] == f, "asv_id"])),
        "coef": float(r.coef),
        "higher_in": r.direction,
        "model": r.model,
        "qval_individual": float(r.qval_individual),
        "qval_joint": float(r.qval_joint),
        "n_samples": int(r.N),
        "n_samples_nonzero": int(r.N_not_zero),
    })

feat = pd.DataFrame(rows)
order_cohort = {c: i for i, c in enumerate(COHORTS)}
order_type = {"multi-ASV": 0, "single ASV": 1}
feat = feat.sort_values(["cohort", "phylotype_type", "qval_joint"],
                        key=lambda s: s.map(order_cohort) if s.name == "cohort"
                        else s.map(order_type) if s.name == "phylotype_type"
                        else s).reset_index(drop=True)

# ---- member table -----------------------------------------------------------
cohorts_of = (feat.groupby("phylotype")["cohort"]
              .apply(lambda s: ";".join(sorted(set(s), key=COHORTS.index))))

mem["is_representative"] = [rep_of[c] == a for c, a in zip(mem["cluster_id"], mem["asv_id"])]
mem_out = pd.DataFrame({
    "phylotype": mem["cluster_id"],
    "asv_id": mem["asv_id"],
    "is_representative": mem["is_representative"],
    "sequence": mem["sequence"],
    "sequence_length": mem["sequence"].str.len(),
    "annotation_distance": mem["annotation_distance"],
    "homd_identity_pct": pd.to_numeric(mem["homd_identity"]),
    "nearest_homd_ref": mem["nearest_homd_ref"],
    "nearest_homd_ref_species": mem["nearest_homd_ref_species"],
    "nearest_homd_ref_distance": pd.to_numeric(mem["nearest_homd_ref_dist"]),
    "n_subjects": pd.to_numeric(mem["subject_prevalence"]),
    "n_bioprojects": pd.to_numeric(mem["n_bioprojects"]),
    "total_reads": pd.to_numeric(mem["total_reads"]),
    "cohorts": mem["cluster_id"].map(cohorts_of),
})
mem_out = mem_out.sort_values(
    ["phylotype", "is_representative", "asv_id"],
    ascending=[True, False, True]).reset_index(drop=True)

# ---- column definitions -----------------------------------------------------
DEFS = [
    ("cohort", "Disease cohort the MaAsLin3 model was fitted in; head and neck "
               "cancer and periodontal disease were modelled separately."),
    ("phylotype", "Phylotype id: c<number> for a phylotype of two or more ASVs, "
                  "s<number> for a single ASV."),
    ("phylotype_type", "multi-ASV or single ASV."),
    ("n_members", "Number of ASVs in the phylotype (1 for a single ASV)."),
    ("taxonomy", "Lineage of the node the annotation comes from, to its deepest "
                 "named rank."),
    ("assigned_taxon", "Name of that deepest named rank; the taxon the phylotype "
                       "is annotated as."),
    ("taxon_rank", "Rank of assigned_taxon."),
    ("phylum", "Phylum from the lineage."),
    ("annotation_source", "Where the annotation comes from: reference sequences "
                          "within the phylotype, the ancestor node of a phylotype "
                          "with no internal reference, or the nearest "
                          "reference-bearing ancestor of a single ASV."),
    ("annotation_distance_rep", "Tree distance (substitutions per site) from the "
                                "representative ASV tip to the reference tip or "
                                "ancestor node that carries the annotation."),
    ("annotation_distance_min", "Smallest annotation distance among the member ASVs."),
    ("annotation_distance_median", "Median annotation distance among the member ASVs."),
    ("annotation_distance_max", "Largest annotation distance among the member ASVs."),
    ("homd_identity_pct", "Identity (%) of the representative ASV to its nearest "
                          "eHOMD reference."),
    ("nearest_homd_ref", "Nearest eHOMD reference sequence of the representative ASV."),
    ("nearest_homd_ref_species", "eHOMD species of that reference."),
    ("nearest_homd_ref_distance", "Tree distance from the representative ASV to that "
                                  "reference."),
    ("max_n_subjects", "Number of subjects in which the phylotype was detected "
                       "(highest among its member ASVs)."),
    ("n_bioprojects", "Number of bioprojects in which the phylotype was detected."),
    ("n_continents", "Number of continents in which the phylotype was detected."),
    ("total_reads", "Reads summed over the member ASVs."),
    ("representative_asv", "ASV id of the representative (the derep sequence the "
                           "phylotype was aggregated on; the ASV itself for a "
                           "single-ASV phylotype)."),
    ("representative_sequence", "16S sequence of the representative ASV."),
    ("representative_sequence_length", "Length of that sequence (bp)."),
    ("member_asvs", "Semicolon-separated ids of all member ASVs; their sequences are "
                    "in the member table."),
    ("coef", "MaAsLin3 coefficient for the disease term."),
    ("higher_in", "Group in which the phylotype is higher: Disease or HC."),
    ("model", "MaAsLin3 model the phylotype was called in: prevalence or abundance."),
    ("qval_individual", "MaAsLin3 individual q value (BH within the cohort)."),
    ("qval_joint", "MaAsLin3 joint q value across both models."),
    ("n_samples", "Samples in the model."),
    ("n_samples_nonzero", "Samples in which the feature was non-zero."),
    ("asv_id", "Member ASV id (member table)."),
    ("is_representative", "TRUE for the representative ASV of its phylotype (member table)."),
    ("sequence", "Full 16S sequence of the member ASV (member table)."),
    ("sequence_length", "Length of that sequence in bp (member table)."),
    ("annotation_distance", "Tree distance from this ASV to the node its annotation "
                            "comes from (member table)."),
    ("n_subjects", "Number of subjects in which this ASV was detected "
                   "(member table)."),
    ("cohorts", "Cohorts in which the phylotype was differential (member table)."),
]

# ---- write ------------------------------------------------------------------
OUTDIR.mkdir(parents=True, exist_ok=True)
f_tsv = OUTDIR / "TableS_differential_phylotypes_HNC_PD.tsv"
m_tsv = OUTDIR / "TableS_differential_phylotypes_HNC_PD_members.tsv"
xlsx = OUTDIR / "TableS_differential_phylotypes_HNC_PD.xlsx"

feat.to_csv(f_tsv, sep="\t", index=False, float_format="%.6g")
mem_out.to_csv(m_tsv, sep="\t", index=False, float_format="%.6g")

with pd.ExcelWriter(xlsx, engine="openpyxl") as xl:
    feat.to_excel(xl, sheet_name="features", index=False)
    mem_out.to_excel(xl, sheet_name="members", index=False)
    pd.DataFrame(DEFS, columns=["column", "definition"]).to_excel(
        xl, sheet_name="column_definitions", index=False)
    for name in ("features", "members", "column_definitions"):
        ws = xl.sheets[name]
        ws.freeze_panes = "A2"
        for cell in ws[1]:
            cell.font = Font(bold=True)
        for col in ws.columns:
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:40])
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)

print(f"[18] features {feat.shape} -> {f_tsv}")
print(f"[18] members  {mem_out.shape} -> {m_tsv}")
print(f"[18] workbook -> {xlsx}")
print(feat.groupby(["cohort", "phylotype_type"]).size().to_string())
print("\nmulti-ASV phylotypes:")
print(feat[feat["phylotype_type"] != "single ASV"][
    ["cohort", "phylotype", "n_members", "assigned_taxon",
     "annotation_distance_rep", "qval_joint"]].to_string(index=False))
