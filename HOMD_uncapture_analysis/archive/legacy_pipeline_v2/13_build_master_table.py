#!/usr/bin/env python3
"""
pipeline_v2 Step 13: master per-ASV supplementary table (CLEAN, new-scheme only).

One row per dereplicated ASV (54,290). Column blocks:
  A. sequence + support      : derep_id, sequence, n_asvs, n_bioprojects,
                               total_reads, prevalence_samples, subject_prevalence
  B. GG2 2024.09 annotation  : best_identity, tier (species/sgb/none),
                               best_qcov, ambiguous, best_target, n_valid_hits,
                               best_species, best_species_token, best_genus,
                               sgb_placeholder
  C. HOMD status             : homd_identity, homd_target, has_homd_hit,
                               uncaptured_reason
  D. geography / regions     : bioprojects, amplicon_regions, sequencing_regions,
                               continents, countries
  E. phylotype (edge x label): phylotype_id, edge_num, unit_label + phylo_*:
                               member/support stats, novelty_class, oral_clade,
                               headline_candidate, supported_clean
  F. contamination flags     : contamination_flag, flag_type

Legacy columns from the old scheme (gg2_genus/name/identity, uchime_*,
is_chimera, homd_uncaptured/alnlen/qcov, phylo_gg2_anchored,
phylo_clean_candidate_novel, phylo_credible_candidate) are intentionally
REMOVED — the new annotation (best_*) supersedes them.

Output: 13_master_asv_table.tsv
"""
import json
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"

log = lambda s: print(f"[13_master] {s}", flush=True)


def main():
    log("Loading step 04 (support floor metadata)...")
    base = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False)
    base = base[["derep_id", "sequence", "n_bioprojects", "total_reads",
                 "prevalence_samples", "contamination_flag", "flag_type",
                 "bioprojects", "amplicon_regions", "sequencing_regions",
                 "continents", "countries"]]
    log(f"  {len(base):,} ASVs")

    log("Loading step 01b (subject prevalence)...")
    subj = pd.read_csv(PIPE_V2 / "01b_derep_subject_prevalence.tsv", sep="\t",
                       low_memory=False, usecols=["derep_id", "subject_prevalence"])
    base = base.merge(subj, on="derep_id", how="left")

    log("Loading step 01 (HOMD mapping + per-derep ASV count)...")
    h01 = pd.read_csv(PIPE_V2 / "01_chimera_filtered.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "n_asvs", "homd_target", "identity",
                               "has_homd_hit", "uncaptured_reason"])
    h01 = h01.rename(columns={"identity": "homd_identity"})
    base = base.merge(h01, on="derep_id", how="left")

    log("Loading step 10 ASV -> phylotype map...")
    mp = pd.read_csv(PIPE_V2 / "10_asv_phylotype_map.tsv", sep="\t", low_memory=False)
    base = base.merge(mp, on="derep_id", how="left")

    log("Loading step 10 GG2 2024.09 annotation...")
    ann = pd.read_csv(PIPE_V2 / "10_gg2_closedref_annotation.tsv", sep="\t",
                      low_memory=False)
    base = base.merge(ann, on="derep_id", how="left")
    base["tier"] = base["tier"].fillna("none")

    log("Loading step 11 phylotype-level evidence...")
    phy = pd.read_csv(PIPE_V2 / "11_phylotypes_homd_novelty.tsv", sep="\t",
                      low_memory=False)
    keep = ["phylotype_id", "n_asvs", "rep_derep_id", "rep_identity", "rep_tier",
            "rep_species", "rep_species_token", "rep_genus", "total_reads_sum",
            "max_subject_prevalence", "n_bioprojects_total", "n_regions",
            "n_continents", "regions", "n_contaminated_asvs", "any_contamination",
            "clean_phylotype", "novelty_class", "oral_clade",
            "renamed_in_homd", "species_in_homd", "headline_candidate",
            "supported_clean"]
    phy = phy[keep].rename(columns={
        "n_asvs": "phylo_n_asvs", "rep_derep_id": "phylo_rep_derep_id",
        "rep_identity": "phylo_rep_identity", "rep_tier": "phylo_rep_tier",
        "rep_species": "phylo_rep_species",
        "rep_species_token": "phylo_rep_species_token",
        "rep_genus": "phylo_rep_genus", "total_reads_sum": "phylo_total_reads",
        "max_subject_prevalence": "phylo_max_subject_prevalence",
        "n_bioprojects_total": "phylo_n_bioprojects",
        "n_regions": "phylo_n_regions", "n_continents": "phylo_n_continents",
        "regions": "phylo_regions",
        "n_contaminated_asvs": "phylo_n_contaminated_asvs",
        "any_contamination": "phylo_any_contamination",
        "clean_phylotype": "phylo_clean", "novelty_class": "phylo_novelty_class",
        "oral_clade": "phylo_oral_clade",
        "renamed_in_homd": "phylo_renamed_in_homd",
        "species_in_homd": "phylo_species_in_homd",
        "headline_candidate": "phylo_headline_candidate",
        "supported_clean": "phylo_supported_clean"})
    base = base.merge(phy, on="phylotype_id", how="left")

    # ---- column order (blocks A-F) ----
    first = ["derep_id", "sequence",
             "n_asvs", "n_bioprojects", "total_reads", "prevalence_samples",
             "subject_prevalence",
             "best_identity", "tier", "best_qcov", "ambiguous", "best_target",
             "n_valid_hits", "best_species", "best_species_token", "best_genus",
             "sgb_placeholder",
             "homd_identity", "homd_target", "has_homd_hit", "uncaptured_reason",
             "bioprojects", "amplicon_regions", "sequencing_regions",
             "continents", "countries",
             "phylotype_id", "edge_num", "unit_label",
             "phylo_n_asvs", "phylo_rep_derep_id", "phylo_rep_identity",
             "phylo_rep_tier", "phylo_rep_species", "phylo_rep_species_token",
             "phylo_rep_genus", "phylo_total_reads", "phylo_max_subject_prevalence",
             "phylo_n_bioprojects", "phylo_n_regions", "phylo_n_continents",
             "phylo_regions", "phylo_n_contaminated_asvs", "phylo_any_contamination",
             "phylo_clean", "phylo_novelty_class", "phylo_oral_clade",
             "phylo_renamed_in_homd", "phylo_species_in_homd",
             "phylo_headline_candidate", "phylo_supported_clean",
             "contamination_flag", "flag_type"]
    rest = [c for c in base.columns if c not in first]
    base = base[first + rest]

    out = PIPE_V2 / "13_master_asv_table.tsv"
    base.to_csv(out, sep="\t", index=False)
    log(f"saved {out} ({len(base):,} rows x {len(base.columns)} cols, "
        f"{out.stat().st_size / 1e6:.1f} MB)")

    summary = {
        "n_rows": int(len(base)),
        "n_columns": int(len(base.columns)),
        "n_with_phylotype": int(base["phylotype_id"].notna().sum()),
        "n_with_gg2_annotation": int(base["best_identity"].notna().sum()),
        "tier": base["tier"].value_counts().to_dict(),
        "n_contamination_flagged": int(base["contamination_flag"].astype(str)
                                       .str.lower().isin(["true", "1"]).sum()),
        "n_ambiguous": int(base["ambiguous"].fillna(False).astype(bool).sum()),
        "n_in_headline_phylotypes": int(
            base["phylo_headline_candidate"].fillna(False).astype(bool).sum()),
        "n_in_supported_clean": int(
            base["phylo_supported_clean"].fillna(False).astype(bool).sum()),
    }
    with open(PIPE_V2 / "13_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(f"summary: {summary}")
    log("Done!")


if __name__ == "__main__":
    main()
