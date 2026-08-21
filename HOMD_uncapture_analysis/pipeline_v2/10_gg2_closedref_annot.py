#!/usr/bin/env python3
"""
pipeline_v2 Step 10: GG2 closed-reference annotation + edge-based phylotypes.

Design (2026-08-18, agreed with user; GG2 2024.09 reference 2026-08-18):
  * Phylotype = ASVs placed on the SAME SEPP edge (no D-value splitting).
    Exception: if members of one edge annotate to DIFFERENT species/SGB at
    identity >= 97%, split them into EDGE_<edge>__<species-token> units so that
    support statistics are never pooled across species.
  * Species annotation = vsearch --usearch_global vs the GG2 2024.09 backbone
    full-length 16S, --iddef 2, --strand both, qcov >= 0.9:
        id >= 99%   -> species-level annotation (named species or GTDB SGB)
        97 <= id < 99 -> SGB-level annotation only (placeholder/stable unit)
        id < 97%    -> NOT annotated (GG2-uncaptured; deprioritized)
  * Ambiguity: if the 2nd-best valid hit is within 0.5% of the best and has a
    different species, the ASV is flagged ambiguous (no species claim).

Inputs:
  - 04_support_floor_passed.tsv            : ASV sequences + support metadata
  - 05_asv_phylotype_assignments.tsv       : SEPP placement (edge_num, phylotype)
  - 01b_derep_subject_prevalence.tsv       : subject-level prevalence
  - 06_novelty/phylotype_reps.fasta        : (unused here; kept for context)
  - /beegfs/dongbiao/greengene2/backbone_fasta/dna-sequences.fasta : GG2 2024.09
  - /beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv     : GG2 2024.09

Outputs (pipeline_v2/):
  - gg2_ref/gg2_2024.09_full_length.fasta         : deduped backbone
  - 10_gg2_hits_iddef2.tsv                        : raw vsearch hits (all)
  - 10_gg2_closedref_annotation.tsv               : per-ASV annotation
  - 10_phylotypes_edge_label.tsv                  : new edge-based phylotypes
  - 10_credible_candidates.tsv                    : GG2-anchored clean candidates
  - 10_summary.json
"""
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

BASE = Path("/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis")
PIPE_V2 = BASE / "pipeline_v2"
OUT = PIPE_V2
REF_DIR = PIPE_V2 / "gg2_ref"
GG2_2024_09_FASTA = Path("/beegfs/dongbiao/greengene2/backbone_fasta/dna-sequences.fasta")
TAXONOMY_TSV = Path("/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv")
VSEARCH = "/home/huangxc/.conda/envs/rachis-qiime2-2026.4/bin/vsearch"
THREADS = int(os.environ.get("SLURM_CPUS_PER_TASK", "28"))

ID_SPECIES = 99.0
ID_SGB = 97.0
QCOV_MIN = 100.0  # closed-reference requires FULL query coverage (vsearch qcov = 0-100)
AMBIG_DELTA = 0.5  # % identity window for ambiguity

log = lambda s: print(f"[10_gg2_annot] {s}", flush=True)


# ============================================================
# 0. Reference preparation (dedupe the GG2 2024.09 backbone)
# ============================================================
def prep_reference():
    fasta = REF_DIR / "gg2_2024.09_full_length.fasta"
    if fasta.exists():
        n = sum(1 for _ in open(fasta) if _.startswith(">"))
        log(f"reference already prepared: {fasta} ({n:,} seqs)")
        return fasta
    REF_DIR.mkdir(parents=True, exist_ok=True)
    log(f"preparing reference from {GG2_2024_09_FASTA}")
    seen = set()
    kept = 0
    dup = 0
    with open(GG2_2024_09_FASTA) as fin, open(fasta, "w") as fout:
        cur_id = None
        cur_seq = []
        def flush():
            nonlocal kept, dup
            if cur_id is None:
                return
            seq = "".join(cur_seq).upper()
            if not seq:
                return
            if seq in seen:
                dup += 1
                return
            seen.add(seq)
            fout.write(f">{cur_id}\n{seq}\n")
            kept += 1
        for line in fin:
            if line.startswith(">"):
                flush()
                cur_id = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        flush()
    log(f"reference ready: {kept:,} unique sequences ({dup:,} duplicate records dropped)")
    return fasta


# ============================================================
# 1. Write ASV fasta (support-floor set)
# ============================================================
def write_asv_fasta(limit=None):
    asv = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                      usecols=["derep_id", "sequence"])
    if limit:
        asv = asv.head(limit)
    fa = OUT / "10_asvs_support_floor.fasta"
    with open(fa, "w") as f:
        for _, r in asv.iterrows():
            f.write(f">{r['derep_id']}\n{r['sequence']}\n")
    log(f"ASV fasta written: {len(asv):,} seqs -> {fa}")
    return fa


# ============================================================
# 2. vsearch closed-reference (cached)
# ============================================================
def run_vsearch(query_fa, ref_fa):
    hits = OUT / "10_gg2_hits_iddef2.tsv"
    if hits.exists() and hits.stat().st_size > 0:
        log(f"vsearch hits cached: {hits}")
        return hits
    log(f"running vsearch (iddef 2, strand both, maxaccepts 8)...")
    subprocess.run([
        VSEARCH, "--usearch_global", str(query_fa),
        "--db", str(ref_fa),
        "--id", "0.0", "--iddef", "2",
        "--strand", "both",
        "--maxaccepts", "8", "--maxrejects", "64",
        "--userout", str(hits),
        "--userfields", "query+target+id+alnlen+mism+opens+qcov+qlo+qhi",
        "--threads", str(THREADS), "--quiet",
    ], check=True)
    log(f"vsearch done: {hits}")
    return hits


# ============================================================
# 3. Parse hits + annotate
# ============================================================
def load_taxonomy():
    RANKS = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
             "f": "family", "g": "genus", "s": "species"}
    tax = {}
    with open(TAXONOMY_TSV) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            if parts[0] == "Feature ID":
                continue
            lineage = {}
            for tok in parts[1].split(";"):
                tok = tok.strip()
                if len(tok) >= 4 and tok[1:3] == "__" and tok[3:]:
                    lineage[RANKS.get(tok[0])] = tok[3:]
            tax[parts[0]] = lineage
    return tax


def parse_hits(hits_path):
    """Return dict: query -> list of (target, id, qcov) sorted by id desc."""
    out = {}
    with open(hits_path) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 7:
                continue
            q, t, ident = p[0], p[1], float(p[2])
            # userfields: query target id alnlen mism opens qcov qlo qhi qlen tlen
            qcov = float(p[6])
            if qcov < QCOV_MIN:
                continue
            out.setdefault(q, []).append((t, ident, qcov))
    for q in out:
        out[q].sort(key=lambda x: -x[1])
    return out


PLACEHOLDER_SP = re.compile(r"^sp\d+$")


def clean_species(s):
    """'Treponema_B_986142 vincentii' -> 'Treponema vincentii'; spNNN -> None."""
    if not isinstance(s, str) or not s:
        return None
    parts = s.split(" ")
    if len(parts) < 2:
        return None
    genus = parts[0].split("_")[0]
    sp = " ".join(parts[1:]).strip()
    if PLACEHOLDER_SP.fullmatch(sp):
        return None
    return f"{genus} {sp}"


def clean_genus(s):
    if not isinstance(s, str) or not s:
        return None
    return s.split("_")[0]


def annotate_asv(hits, tax):
    rows = []
    for q, hl in hits.items():
        best_t, best_id, best_qcov = hl[0]
        lineage = tax.get(best_t, {})
        best_species_raw = lineage.get("species")
        best_genus = clean_genus(lineage.get("genus"))
        best_species = clean_species(best_species_raw)
        is_placeholder = bool(best_species_raw and PLACEHOLDER_SP.fullmatch(best_species_raw.split(" ")[-1]))

        if best_id >= ID_SPECIES:
            tier = "species"
        elif best_id >= ID_SGB:
            tier = "sgb"
        else:
            tier = "none"

        # ambiguity: 2nd-best valid hit within AMBIG_DELTA of best, different species
        ambiguous = False
        if len(hl) > 1:
            t2, id2, _ = hl[1]
            if best_id - id2 <= AMBIG_DELTA:
                sp2 = clean_species(tax.get(t2, {}).get("species"))
                if sp2 != best_species or (best_species is None and sp2 is not None):
                    ambiguous = True

        rows.append({
            "derep_id": q,
            "best_identity": round(best_id, 2),
            "best_target": best_t,
            "best_qcov": round(best_qcov, 3),
            "n_valid_hits": len(hl),
            "tier": tier,
            "ambiguous": ambiguous,
            "best_genus": best_genus,
            "best_species_token": best_species_raw,
            "best_species": best_species,
            "sgb_placeholder": is_placeholder,
        })
    return pd.DataFrame(rows)


# ============================================================
# 4. Phylotype = edge (split by species label when mixed)
# ============================================================
def build_phylotypes(ann):
    placed = pd.read_csv(PIPE_V2 / "05_asv_phylotype_assignments.tsv", sep="\t",
                         low_memory=False,
                         usecols=["derep_id", "edge_num", "distal_length", "pendant_length"])
    meta = pd.read_csv(PIPE_V2 / "04_support_floor_passed.tsv", sep="\t", low_memory=False,
                       usecols=["derep_id", "contamination_flag", "n_bioprojects",
                                "total_reads", "amplicon_regions"])
    subj = pd.read_csv(PIPE_V2 / "01b_derep_subject_prevalence.tsv", sep="\t",
                       low_memory=False, usecols=["derep_id", "subject_prevalence"])
    df = placed.merge(meta, on="derep_id", how="left").merge(subj, on="derep_id", how="left")
    df = df.merge(ann, on="derep_id", how="left")
    df["edge_num"] = df["edge_num"].astype("Int64").astype(str)
    df["contamination_flag"] = df["contamination_flag"].astype(str).str.lower().isin(["true", "1", "t"])
    df["subject_prevalence"] = pd.to_numeric(df["subject_prevalence"], errors="coerce").fillna(0).astype(int)
    df["total_reads"] = pd.to_numeric(df["total_reads"], errors="coerce").fillna(0)

    # label used for splitting: >=97% species/SGB token (or genus/target if the
    # species field is empty). Unannotated ASVs get a reserved UNRESOLVED label
    # so they never pool support with annotated members of the same edge.
    UNRESOLVED = "__UNRESOLVED__"

    def unit_label(r):
        if pd.isna(r["tier"]) or r["tier"] == "none":
            return UNRESOLVED
        tok = None
        for cand in (r["best_species_token"], r["best_genus"], r["best_target"]):
            if isinstance(cand, str) and cand.strip():
                tok = cand
                break
        if tok is None:
            return UNRESOLVED
        return str(tok).replace(" ", "_")

    df["unit_label"] = df.apply(unit_label, axis=1)

    # primary phylotype = edge; split when a single edge has >1 distinct label
    edge_labels = df.groupby("edge_num")["unit_label"].nunique()
    split_edges = set(edge_labels[edge_labels > 1].index)

    def phy_id(r):
        e = r["edge_num"]
        if e in split_edges:
            lab = r["unit_label"]
            if lab == UNRESOLVED:
                lab = "unresolved"
            return f"EDGE_{e}__{lab}"
        return f"EDGE_{e}"
    df["phylotype_id"] = df.apply(phy_id, axis=1)
    log(f"edges split by species label: {len(split_edges):,} of "
        f"{df['edge_num'].nunique():,} with ASVs")
    df[["derep_id", "phylotype_id", "edge_num", "unit_label"]].to_csv(
        OUT / "10_asv_phylotype_map.tsv", sep="\t", index=False)
    log("saved 10_asv_phylotype_map.tsv (ASV -> phylotype_id)")

    def parse_set(val):
        if pd.isna(val) or val == "":
            return set()
        s = str(val).strip("{}")
        items = []
        for sep in ["; ", ";", ", ", ","]:
            if sep in s:
                items = [x.strip().strip("'\"") for x in s.split(sep)]
                break
        if not items:
            items = [s.strip().strip("'\"")]
        return set(i for i in items if i)

    df = df.sort_values(["subject_prevalence", "total_reads"], ascending=False)
    phylo = df.groupby("phylotype_id").agg(
        n_asvs=("derep_id", "count"),
        rep_derep_id=("derep_id", "first"),
        edge_num=("edge_num", "first"),
        rep_identity=("best_identity", "first"),
        rep_tier=("tier", "first"),
        rep_species=("best_species", "first"),
        rep_species_token=("best_species_token", "first"),
        rep_genus=("best_genus", "first"),
        total_reads_sum=("total_reads", "sum"),
        max_subject_prevalence=("subject_prevalence", "max"),
        n_contaminated_asvs=("contamination_flag", "sum"),
        any_contamination=("contamination_flag", "any"),
    ).reset_index()
    phylo["n_bioprojects_total"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["n_bioprojects"].apply(
            lambda x: len(set(r for v in x.dropna() for r in parse_set(v)))))
    phylo["regions"] = phylo["phylotype_id"].map(
        df.groupby("phylotype_id")["amplicon_regions"].apply(
            lambda x: ";".join(sorted(set(
                r for v in x.dropna() for r in parse_set(v))))))
    phylo["n_regions"] = phylo["regions"].apply(
        lambda s: len(s.split(";")) if s and s.strip() else 0)
    phylo["clean_phylotype"] = phylo["n_contaminated_asvs"] == 0
    phylo["gg2_anchored"] = phylo["rep_identity"] >= ID_SGB
    return df, phylo


# ============================================================
# 5. Credible candidates (GG2-anchored, HOMD-uncaptured, clean)
# ============================================================
def credible_candidates(phylo):
    # All support-floor ASVs are HOMD-uncaptured by construction (step 01);
    # phylotype-level gate = clean + GG2-anchored (>=97%) + subject support.
    phylo["homd_uncaptured"] = True
    phylo["clean_candidate_novel"] = (
        phylo["homd_uncaptured"] & phylo["clean_phylotype"] & phylo["gg2_anchored"])
    # default support gate: >=10 subjects (tunable downstream)
    phylo["credible_candidate"] = (
        phylo["clean_candidate_novel"] & (phylo["max_subject_prevalence"] >= 10))
    return phylo


def main():
    limit = None
    if len(sys.argv) > 1 and sys.argv[1] == "--limit":
        limit = int(sys.argv[2])
    ref = prep_reference()
    qfa = write_asv_fasta(limit)
    hits = run_vsearch(qfa, ref)
    tax = load_taxonomy()
    ann = annotate_asv(parse_hits(hits), tax)
    ann.to_csv(OUT / "10_gg2_closedref_annotation.tsv", sep="\t", index=False)
    log(f"annotation saved: {len(ann):,} ASVs "
        f"| species: {(ann['tier']=='species').sum():,} "
        f"| sgb: {(ann['tier']=='sgb').sum():,} "
        f"| none: {(ann['tier']=='none').sum():,}")

    df, phylo = build_phylotypes(ann)
    phylo = credible_candidates(phylo)
    phylo.to_csv(OUT / "10_phylotypes_edge_label.tsv", sep="\t", index=False)
    log(f"phylotypes: {len(phylo):,} | clean: {int(phylo['clean_phylotype'].sum()):,} "
        f"| gg2_anchored: {int(phylo['gg2_anchored'].sum()):,} "
        f"| credible (>=10 subjects): {int(phylo['credible_candidate'].sum()):,}")

    credible = phylo[phylo["credible_candidate"]].sort_values(
        ["max_subject_prevalence", "n_bioprojects_total", "n_regions"],
        ascending=False)
    credible.to_csv(OUT / "10_credible_candidates.tsv", sep="\t", index=False)

    summary = {
        "n_asvs_annotated": int(len(ann)),
        "tier_species": int((ann["tier"] == "species").sum()),
        "tier_sgb": int((ann["tier"] == "sgb").sum()),
        "tier_none": int((ann["tier"] == "none").sum()),
        "ambiguous_species": int(ann["ambiguous"].sum()),
        "n_phylotypes": int(len(phylo)),
        "n_phylotypes_clean": int(phylo["clean_phylotype"].sum()),
        "n_phylotypes_gg2_anchored": int(phylo["gg2_anchored"].sum()),
        "n_credible_candidates_subj10": int(phylo["credible_candidate"].sum()),
        "cutoffs": {"species_id": ID_SPECIES, "sgb_id": ID_SGB, "qcov_min": QCOV_MIN,
                    "ambiguity_delta": AMBIG_DELTA},
    }
    with open(OUT / "10_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log(f"summary: {summary}")
    log("Done!")


if __name__ == "__main__":
    main()
