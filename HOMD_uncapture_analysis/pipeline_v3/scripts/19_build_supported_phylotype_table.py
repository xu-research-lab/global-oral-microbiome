#!/usr/bin/env python3
"""Supplementary table: the supported catalog, one row per phylotype.

The source is 08_supported_clusters.tsv, the cross-continent supported catalog
(subject prevalence > 50, at least two bioprojects, at least two continents,
no contamination-flagged member, no organelle member).  Two kinds of change are
made to it:

  * "cluster" becomes "phylotype" throughout, and the columns are renamed to
    say the same thing in words rather than in pipeline keys;
  * columns that carry no information in this catalog are dropped -- every row
    is clean and free of organelles by construction, `headline` is the same
    boolean as "eHOMD-uncaptured", `candidate_novel` adds nothing to it, and
    `max_n_bioprojects` is the per-member maximum where the table keeps the
    union.

The taxonomy comes from the reference-tree lineage, cut at its deepest named
rank; this names the 69 phylotypes whose `annotation` is empty because the
ancestor node carried no label at the rank the pipeline stopped at (their
`annotation_level` then overstates what is known).

Row order is the source order: support descending (subjects, then bioprojects,
then amplicon regions).

Run from pipeline_v3.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.styles import Font

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/"
            "HOMD_uncapture_analysis/pipeline_v3")
OUTDIR = BASE / "analysis_report" / "tables"

RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
BASIS = {
    "cluster_refs": "reference sequences within the phylotype",
    "cluster_ancestor": "ancestor node of the phylotype",
    "singleton_ancestor": "nearest reference-bearing ancestor",
}
TIER = {
    "<=0.01 (~99%)": "≥99%",
    "0.01-0.03 (97-99%)": "97–99%",
    ">0.03 (<97%)": "<97%",
}
NOVELTY = {
    "homd_captured": "eHOMD-captured",
    "homd_uncaptured": "eHOMD-uncaptured",
}

d = pd.read_csv(BASE / "08_supported_clusters.tsv", sep="\t", dtype=str)

# ---- taxonomy ---------------------------------------------------------------
def split_lineage(s):
    parts = [p.strip() for p in str(s).split(";")] if isinstance(s, str) else []
    parts = (parts + [""] * len(RANKS))[:len(RANKS)]
    named = [i for i, p in enumerate(parts) if p]
    if not named:
        return "", "", "", ""
    deepest = named[-1]
    shown = [p if p else "unclassified" for p in parts[:deepest + 1]]
    return "; ".join(shown), parts[deepest], RANKS[deepest], parts[1]

lin = [split_lineage(s) for s in d["annotation_lineage"]]
taxonomy = [x[0] for x in lin]
deepest = [x[1] for x in lin]
rank = [x[2] for x in lin]
phylum = [x[3] for x in lin]
assigned = [a if isinstance(a, str) and a else b
            for a, b in zip(d["annotation"], deepest)]

num = lambda c: pd.to_numeric(d[c], errors="coerce").astype("Int64")
flt = lambda c: pd.to_numeric(d[c], errors="coerce")

out = pd.DataFrame({
    "phylotype": d["cluster_id"],
    "phylotype_type": np.where(d["n_asvs"].astype(int) >= 2,
                               "multi-ASV", "single ASV"),
    "n_members": num("n_asvs"),
    "n_references": num("n_refs_in_cluster"),
    "representative_asv": d["rep_derep_id"],
    "member_asvs": d["member_derep_ids"],
    "taxonomy": taxonomy,
    "assigned_taxon": assigned,
    "taxon_rank": rank,
    "phylum": phylum,
    "annotation_source": d["annotation_basis"].map(BASIS),
    "annotation_distance_min": flt("annotation_distance_min"),
    "annotation_distance_max": flt("annotation_distance_max"),
    "novelty_class": d["cluster_novelty_class"].map(NOVELTY),
    "n_members_annotation_in_homd": num("n_annotation_homd"),
    "genus_in_homd": d["any_annotation_genus_homd"].map(
        {"True": True, "False": False}),
    "max_n_subjects": num("max_subject_prevalence"),
    "max_n_samples": num("max_prevalence"),
    "n_bioprojects": num("n_bioprojects_total"),
    "n_continents": num("n_continents"),
    "continents": d["continents"],
    "n_amplicon_regions": num("n_regions"),
    "amplicon_regions": d["regions"],
    "total_reads": num("total_reads_sum"),
    "homd_identity_pct": flt("rep_homd_identity"),
    "nearest_homd_target": d["rep_homd_target"],
    "nearest_reference_distance": flt("rep_nearest_ref_dist"),
    "nearest_reference_identity_tier": d["tree_tier"].map(TIER),
    "support_score": num("support_score"),
})

if out["phylotype"].duplicated().any():
    raise SystemExit("duplicated phylotype ids")
if out["annotation_source"].isna().any() or out["novelty_class"].isna().any():
    raise SystemExit("unmapped categorical value")
if out["taxonomy"].eq("").any():
    raise SystemExit("a phylotype has no named rank at all")

DEFS = [
    ("phylotype", "Phylotype id. c<number> for a phylotype of two or more ASVs; "
                  "s<number> for a single ASV that clustered with nothing else; a "
                  "plain number for a single-ASV phylotype that contains reference "
                  "sequences (its treecluster id is positive, so there is no sign "
                  "to replace)."),
    ("phylotype_type", "multi-ASV or single ASV."),
    ("n_members", "Number of member ASVs."),
    ("n_references", "Number of reference sequences inside the phylotype."),
    ("representative_asv", "Member ASV the phylotype was aggregated on."),
    ("member_asvs", "Semicolon-separated ids of all member ASVs."),
    ("taxonomy", "Reference-tree lineage, cut at its deepest named rank."),
    ("assigned_taxon", "Name of that deepest named rank; the taxon the phylotype "
                       "is annotated as. Filled from the pipeline annotation where "
                       "one was made."),
    ("taxon_rank", "Rank of assigned_taxon."),
    ("phylum", "Phylum from the lineage."),
    ("annotation_source", "Where the annotation comes from: reference sequences "
                          "within the phylotype, the ancestor node of a phylotype "
                          "with no internal reference, or the nearest "
                          "reference-bearing ancestor of a single ASV."),
    ("annotation_distance_min", "Smallest tree distance (substitutions per site) "
                                "from a member ASV to the node carrying the "
                                "annotation."),
    ("annotation_distance_max", "Largest such distance."),
    ("novelty_class", "eHOMD-captured when a member is annotated as an eHOMD "
                      "species, eHOMD-uncaptured otherwise."),
    ("n_members_annotation_in_homd", "Number of member ASVs annotated as an eHOMD "
                                     "species."),
    ("genus_in_homd", "TRUE when a member's annotated genus occurs in eHOMD."),
    ("max_n_subjects", "Highest number of subjects in which a member ASV was "
                       "detected."),
    ("max_n_samples", "Highest number of samples in which a member ASV was "
                      "detected."),
    ("n_bioprojects", "Number of distinct bioprojects the members occur in."),
    ("n_continents", "Number of continents."),
    ("continents", "Semicolon-separated continents."),
    ("n_amplicon_regions", "Number of 16S amplicon regions."),
    ("amplicon_regions", "Semicolon-separated amplicon regions."),
    ("total_reads", "Reads summed over the member ASVs."),
    ("homd_identity_pct", "Identity (%) of the representative ASV to its nearest "
                          "eHOMD reference."),
    ("nearest_homd_target", "HMT target of that eHOMD reference."),
    ("nearest_reference_distance", "Tree distance from the representative ASV to "
                                   "the nearest reference sequence."),
    ("nearest_reference_identity_tier", "That distance as the identity band used "
                                        "in the text (≥99%, 97–99%, <97%)."),
    ("support_score", "max_n_subjects x n_bioprojects x n_amplicon_regions."),
]

# ---- write ------------------------------------------------------------------
OUTDIR.mkdir(parents=True, exist_ok=True)
tsv = OUTDIR / "TableS_supported_phylotypes.tsv"
xlsx = OUTDIR / "TableS_supported_phylotypes.xlsx"

out.to_csv(tsv, sep="\t", index=False, float_format="%.6g")

with pd.ExcelWriter(xlsx, engine="openpyxl") as xl:
    out.to_excel(xl, sheet_name="phylotypes", index=False)
    pd.DataFrame(DEFS, columns=["column", "definition"]).to_excel(
        xl, sheet_name="column_definitions", index=False)
    for name in ("phylotypes", "column_definitions"):
        ws = xl.sheets[name]
        ws.freeze_panes = "A2"
        for cell in ws[1]:
            cell.font = Font(bold=True)
        for col in ws.columns:
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:60])
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60)

print(f"[19] {out.shape} -> {tsv}")
print(f"[19] workbook -> {xlsx}")
print(out["phylotype_type"].value_counts().to_string())
print(out["novelty_class"].value_counts().to_string())
print(out["annotation_source"].value_counts().to_string())
print(out["taxon_rank"].value_counts().to_string())
print("\nfilled from the lineage (pipeline annotation empty):",
      int((d["annotation"].isna()).sum()))
