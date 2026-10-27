#!/usr/bin/env python3
"""
pipeline_v2 Step 11: HOMD novelty re-check (new phylotype units) + STRICT tier.

The new phylotypes (Step 10, edge x species-label) replaced the old scheme, so
the HOMD novelty classification is recomputed per new phylotype rep:
  * homd_best_id    : rep ASV's best HOMD identity (cached in 01_chimera_filtered)
  * species_in_homd : rep's GG2 2024.09 species name present in the HOMD catalog
  * renamed_in_homd : GTDB genus != HOMD genus but the epithet matches exactly
                      one HOMD species (e.g. Eubacterium yurii = Peptoanaerobacter
                      yurii) -> treated as captured
  * oral_clade      : rep's genus in the HOMD genus set
  * novelty_class   : homd_captured (species in HOMD, >=97% HOMD hit, or renamed)
                      vs homd_uncaptured (candidate missing from HOMD)

Result layers (2026-08-18, user decision):
  * SUPPORT GATE (single standard, user 2026-08-18): subjects > 50 &
    n_bioprojects_total >= 2 & n_continents >= 2  (cross-CONTINENT evidence,
    NOT amplicon regions)  (no extra BP filter layer)
  * SUPPORTED CATALOG (user 2026-08-18): support gate AND clean_phylotype
    (no contamination-flagged member) -> 186
  * HEADLINE (the scientific claim): supported catalog + GG2 2024.09 identity
    >= 99% AND homd_uncaptured (species absent from HOMD) -> 53

Outputs (pipeline_v2/):
  - 11_phylotypes_homd_novelty.tsv   : all new phylotypes + HOMD novelty columns
  - 11_supported_clean.tsv           : SUPPORTED CATALOG (clean, 186)
  - 11_headline_candidates.tsv       : HEADLINE layer (53)
  - 11_summary.json
"""
import json
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
HOMD_FASTA = Path("/beegfs/huangxiaochang/global_oral/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta")

SUPPORT = {"subject": 50, "bioprojects": 2, "continents": 2}   # cross-continent
ID_SPECIES = 99.0   # user's GG2 closed-reference annotation/filter standard

log = lambda s: print(f"[11_homd_novelty] {s}", flush=True)


def load_homd_catalog():
    homd_genus = set()
    homd_species = set()
    homd_epithets = {}
    n = 0
    with open(HOMD_FASTA) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.split()
                if len(parts) >= 3:
                    n += 1
                    homd_genus.add(parts[2].strip())
                    if len(parts) >= 4:
                        sp = f"{parts[2].strip()} {parts[3].strip()}"
                        homd_species.add(sp)
                        homd_epithets.setdefault(parts[3].strip(), set()).add(sp)
    log(f"HOMD refs: {n:,} | genera: {len(homd_genus):,} | species: {len(homd_species):,}")
    return homd_genus, homd_species, homd_epithets


def clean_species(s):
    if not isinstance(s, str) or not s:
        return None
    parts = s.split(" ")
    if len(parts) < 2:
        return None
    genus = parts[0].split("_")[0]
    sp = " ".join(parts[1:]).strip()
    if sp.startswith("sp") and sp[2:].isdigit():
        return None
    return f"{genus} {sp}"


def main():
    log("Loading new phylotypes (step 10)...")
    phylo = pd.read_csv(PIPE_V2 / "10_phylotypes_edge_label.tsv", sep="\t", low_memory=False)
    log(f"  {len(phylo):,} phylotypes")

    log("Loading per-ASV HOMD mapping (step 01)...")
    homd01 = pd.read_csv(PIPE_V2 / "01_chimera_filtered.tsv", sep="\t", low_memory=False,
                         usecols=["derep_id", "homd_target", "identity", "homd_uncaptured"])
    homd01["homd_uncaptured"] = homd01["homd_uncaptured"].astype(str).str.lower().isin(["true", "1", "t"])
    homd01 = homd01.rename(columns={"identity": "homd_best_id_01"})
    phylo = phylo.merge(homd01, left_on="rep_derep_id", right_on="derep_id", how="left")

    log("Computing per-phylotype continent count...")
    mp = pd.read_csv(PIPE_V2 / "10_asv_phylotype_map.tsv", sep="\t", low_memory=False,
                     usecols=["derep_id", "phylotype_id"])
    s04 = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "continents"])

    def parse_set(val):
        if pd.isna(val) or val == "":
            return set()
        s = str(val).strip("{}")
        for sep in ["; ", ";", ", ", ","]:
            if sep in s:
                return set(x.strip().strip("'\"") for x in s.split(sep))
        return {s.strip().strip("'\"")}

    cont_count = {}
    for pid, g in mp.groupby("phylotype_id"):
        c = set()
        for v in s04.loc[s04["derep_id"].isin(g["derep_id"]), "continents"].dropna():
            c |= parse_set(v)
        cont_count[pid] = len(c)
    phylo["n_continents"] = phylo["phylotype_id"].map(cont_count).fillna(0).astype(int)
    log(f"  n_continents merged (mean {phylo['n_continents'].mean():.1f})")

    log("Loading HOMD catalog (species/genus sets)...")
    homd_genus, homd_species, homd_epithets = load_homd_catalog()

    def species_renamed_in_homd(sp):
        if not isinstance(sp, str):
            return False
        epithet = sp.split(" ")[-1].strip()
        return len(homd_epithets.get(epithet, ())) == 1

    rows = []
    for _, r in phylo.iterrows():
        gg2_sp = clean_species(r["rep_species"])          # e.g. "Neisseria subflava"
        hid = pd.to_numeric(r["homd_best_id_01"], errors="coerce")
        hid = float(hid) if pd.notna(hid) else 0.0
        rep_genus = r["rep_genus"]
        if pd.isna(rep_genus) or not str(rep_genus).strip():
            rep_genus = gg2_sp.split(" ")[0] if gg2_sp else None

        if gg2_sp is not None and gg2_sp in homd_species:
            cls, renamed = "homd_captured", False
        elif hid >= 97.0:
            cls, renamed = "homd_captured", False
        elif gg2_sp is not None and species_renamed_in_homd(gg2_sp):
            cls, renamed = "homd_captured", True
        else:
            cls, renamed = "homd_uncaptured", False

        rows.append({
            "homd_best_id": round(hid, 2),
            "homd_target": r["homd_target"],
            "species_in_homd": bool(gg2_sp is not None and gg2_sp in homd_species),
            "renamed_in_homd": renamed,
            "oral_clade": bool(rep_genus is not None and rep_genus in homd_genus),
            "novelty_class": cls,
        })
    nov = pd.DataFrame(rows, index=phylo.index)
    phylo = pd.concat([phylo, nov], axis=1)
    phylo["homd_uncaptured_new"] = phylo["novelty_class"] == "homd_uncaptured"

    log("Novelty distribution:")
    print(phylo["novelty_class"].value_counts().to_string())
    print(pd.crosstab(phylo["novelty_class"], phylo["oral_clade"]).to_string())

    # ---------------- support gate (single standard) ----------------
    support_gate = ((phylo["max_subject_prevalence"] > SUPPORT["subject"])
                    & (phylo["n_bioprojects_total"] >= SUPPORT["bioprojects"])
                    & (phylo["n_continents"] >= SUPPORT["continents"]))
    log(f"\nSUPPORT GATE (S>50 & BP>=2 & continents>=2): {int(support_gate.sum()):,}")
    supported_clean = support_gate & phylo["clean_phylotype"]
    log(f"SUPPORTED CATALOG (support gate + clean): {int(supported_clean.sum()):,}")

    # ---------------- HEADLINE layer ----------------
    headline = (supported_clean
                & (phylo["rep_identity"] >= ID_SPECIES)
                & phylo["homd_uncaptured_new"])
    log(f"HEADLINE (support gate + GG2>=99% + HOMD-uncaptured): {int(headline.sum()):,}")

    # rep sequence for the paper table
    seq4 = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                       usecols=["derep_id", "sequence"])
    phylo = phylo.merge(seq4.rename(columns={"derep_id": "rep_derep_id"}),
                        on="rep_derep_id", how="left")
    phylo = phylo.rename(columns={"sequence": "rep_sequence"})

    # member lists
    mp = pd.read_csv(PIPE_V2 / "10_asv_phylotype_map.tsv", sep="\t", low_memory=False)
    members = mp.groupby("phylotype_id")["derep_id"].apply(lambda x: ";".join(x)).rename("member_derep_ids")
    phylo = phylo.merge(members, on="phylotype_id", how="left")

    phylo["headline_candidate"] = headline
    phylo["supported_clean"] = supported_clean
    phylo.to_csv(PIPE_V2 / "11_phylotypes_homd_novelty.tsv", sep="\t", index=False)
    log(f"saved 11_phylotypes_homd_novelty.tsv ({len(phylo):,} rows)")

    sup_clean = phylo[supported_clean].sort_values(
        ["max_subject_prevalence", "n_bioprojects_total", "n_regions"],
        ascending=False)
    sup_clean.to_csv(PIPE_V2 / "11_supported_clean.tsv", sep="\t", index=False)
    log(f"saved 11_supported_clean.tsv ({len(sup_clean):,} rows)")

    cand = phylo[headline].sort_values(
        ["max_subject_prevalence", "n_bioprojects_total", "n_regions"],
        ascending=False)
    cand.to_csv(PIPE_V2 / "11_headline_candidates.tsv", sep="\t", index=False)
    log(f"saved 11_headline_candidates.tsv (HEADLINE, {len(cand):,} rows)")

    summary = {
        "n_phylotypes": int(len(phylo)),
        "novelty_class": {c: int((phylo["novelty_class"] == c).sum())
                          for c in ["homd_captured", "homd_uncaptured"]},
        "support_cutoffs": SUPPORT,
        "n_supported": int(support_gate.sum()),
        "n_supported_clean": int(supported_clean.sum()),
        "headline": int(headline.sum()),
        "headline_by_oral_clade": {
            "in_homd_genus": int((headline & phylo["oral_clade"]).sum()),
            "outside_homd_genus": int((headline & ~phylo["oral_clade"]).sum()),
        },
    }
    with open(PIPE_V2 / "11_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(f"summary: {summary}")
    log("Done!")


if __name__ == "__main__":
    main()
