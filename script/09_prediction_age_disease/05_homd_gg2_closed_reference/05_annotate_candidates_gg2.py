#!/usr/bin/env python3
"""Annotate reproducible HOMD-low candidates with Greengenes2 taxonomy."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from collections import Counter, defaultdict

import pandas as pd

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
DEFAULT_GG2_REF = "/beegfs/dongbiao/greengene2/2024.09.seqs.fna"
DEFAULT_GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"


def log(msg):
    print(f"[gg2] {msg}", flush=True)


def write_fasta(df, path):
    with open(path, "w") as handle:
        for row in df.itertuples(index=False):
            handle.write(f">{row.derep_id}\n")
            seq = row.sequence
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")


def load_taxonomy(path):
    tax = {}
    with open(path) as handle:
        header = next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                tax[parts[0]] = parts[1]
    return tax


def parse_rank(taxonomy, prefix):
    if pd.isna(taxonomy):
        return "Unassigned"
    for part in str(taxonomy).split(";"):
        part = part.strip()
        if part.startswith(prefix):
            value = part.split("__", 1)[1].strip()
            return value if value else "Unassigned"
    return "Unassigned"


def classify_qc(taxonomy):
    t = str(taxonomy)
    levels = [x.strip() for x in t.split(";")]
    if not t or t == "nan" or t == "Unassigned":
        return "unassigned"
    if any("Chloroplast" in x for x in levels):
        return "chloroplast"
    if any("Mitochondria" in x for x in levels):
        return "mitochondria"
    if levels and levels[0].startswith("d__Bacteria"):
        return "bacterial"
    if levels and levels[0].startswith("d__Archaea"):
        return "archaea"
    if levels and levels[0].startswith("d__Eukaryota"):
        return "eukaryota"
    return "other_non_bacterial"


def split_values(value):
    if pd.isna(value):
        return ["NA"]
    vals = [v for v in str(value).split(";") if v and v != "nan"]
    return vals or ["NA"]


def summarize_group(df, group_col, tax_col, out_path):
    rows = []
    for _, row in df.iterrows():
        for group in split_values(row[group_col]):
            rows.append({
                "grouping": group_col,
                "level": group,
                "taxon_rank": tax_col,
                "taxon": row[tax_col],
                "derep_sequences": 1,
                "source_asvs": row["n_asvs"],
                "total_reads": row["total_reads"],
                "prevalence_samples": row["prevalence_samples"],
                "n_bioprojects_sum": row["n_bioprojects"],
            })
    if not rows:
        pd.DataFrame().to_csv(out_path, sep="\t", index=False)
        return
    ex = pd.DataFrame(rows)
    agg = ex.groupby(["grouping", "level", "taxon_rank", "taxon"], dropna=False).agg(
        derep_sequences=("derep_sequences", "sum"),
        source_asvs=("source_asvs", "sum"),
        total_reads=("total_reads", "sum"),
        prevalence_samples=("prevalence_samples", "sum"),
        n_bioprojects_sum=("n_bioprojects_sum", "sum"),
    ).reset_index()
    agg = agg.sort_values(["grouping", "level", "total_reads", "derep_sequences"], ascending=[True, True, False, False])
    agg.to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--vsearch", default=DEFAULT_VSEARCH)
    parser.add_argument("--gg2-ref", default=DEFAULT_GG2_REF)
    parser.add_argument("--gg2-tax", default=DEFAULT_GG2_TAX)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--min-id", type=float, default=0.7)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cand_path = os.path.join(args.outdir, "qcov_sensitivity_reproducible_abundant_not_high_identity_homd.tsv")
    fasta_path = os.path.join(args.outdir, "candidates_not_homd_id99_qcov99.fa")
    hit_path = os.path.join(args.outdir, "candidates_vs_gg2_best_hits.tsv")
    annotated_path = os.path.join(args.outdir, "candidates_not_homd_id99_qcov99_gg2_annotated.tsv")
    qc_summary_path = os.path.join(args.outdir, "candidate_gg2_qc_summary.tsv")
    tax_summary_path = os.path.join(args.outdir, "candidate_gg2_taxonomy_summary.tsv")
    group_tax_summary_path = os.path.join(args.outdir, "candidate_gg2_group_taxonomy_summary.tsv")
    stats_path = os.path.join(args.outdir, "candidate_gg2_annotation_stats.json")

    log("loading candidates")
    cand = pd.read_csv(cand_path, sep="\t")
    write_fasta(cand, fasta_path)
    log(f"wrote FASTA: {fasta_path}")

    if args.force or not (os.path.exists(hit_path) and os.path.getsize(hit_path) > 0):
        cmd = [
            args.vsearch, "--usearch_global", fasta_path,
            "--db", args.gg2_ref,
            "--id", str(args.min_id),
            "--maxaccepts", "1", "--maxrejects", "0", "--strand", "both",
            "--userout", hit_path,
            "--userfields", "query+target+id+alnlen+qcov+tcov+ql+tl",
            "--threads", str(args.threads), "--quiet",
        ]
        log("running GG2 alignment")
        log(" ".join(cmd))
        subprocess.run(cmd, check=True)
    else:
        log(f"using existing GG2 alignment: {hit_path}")

    hit_cols = ["derep_id", "gg2_target", "gg2_identity", "gg2_alnlen", "gg2_qcov", "gg2_tcov", "query_length", "target_length"]
    if os.path.getsize(hit_path) > 0:
        hits = pd.read_csv(hit_path, sep="\t", names=hit_cols)
        hits = hits.sort_values(["derep_id", "gg2_identity", "gg2_qcov", "gg2_alnlen"], ascending=[True, False, False, False])
        hits = hits.drop_duplicates("derep_id", keep="first")
    else:
        hits = pd.DataFrame(columns=hit_cols)

    log("loading taxonomy")
    tax = load_taxonomy(args.gg2_tax)
    hits["gg2_taxonomy"] = hits["gg2_target"].map(tax).fillna("Unassigned")
    ann = cand.merge(hits, on="derep_id", how="left")
    ann["gg2_taxonomy"] = ann["gg2_taxonomy"].fillna("Unassigned")
    ann["gg2_qc_class"] = ann["gg2_taxonomy"].map(classify_qc)
    for col, prefix in [
        ("gg2_domain", "d__"), ("gg2_phylum", "p__"), ("gg2_class", "c__"),
        ("gg2_order", "o__"), ("gg2_family", "f__"), ("gg2_genus", "g__"), ("gg2_species", "s__"),
    ]:
        ann[col] = ann["gg2_taxonomy"].map(lambda x, p=prefix: parse_rank(x, p))

    ann.to_csv(annotated_path, sep="\t", index=False)

    qc = ann.groupby("gg2_qc_class", dropna=False).agg(
        derep_sequences=("derep_id", "count"),
        source_asvs=("n_asvs", "sum"),
        total_reads=("total_reads", "sum"),
        prevalence_samples=("prevalence_samples", "sum"),
    ).reset_index().sort_values("total_reads", ascending=False)
    qc.to_csv(qc_summary_path, sep="\t", index=False)

    tax_rows = []
    for rank in ["gg2_phylum", "gg2_class", "gg2_order", "gg2_family", "gg2_genus", "gg2_species"]:
        tmp = ann.groupby(rank, dropna=False).agg(
            derep_sequences=("derep_id", "count"),
            source_asvs=("n_asvs", "sum"),
            total_reads=("total_reads", "sum"),
            prevalence_samples=("prevalence_samples", "sum"),
            max_bioprojects=("n_bioprojects", "max"),
        ).reset_index().rename(columns={rank: "taxon"})
        tmp.insert(0, "rank", rank.replace("gg2_", ""))
        tax_rows.append(tmp)
    tax_summary = pd.concat(tax_rows, ignore_index=True).sort_values(["rank", "total_reads"], ascending=[True, False])
    tax_summary.to_csv(tax_summary_path, sep="\t", index=False)

    bacteria = ann[ann["gg2_qc_class"].eq("bacterial")].copy()
    group_parts = []
    for group_col in ["continents", "countries", "disease_studies", "statuses", "sequencing_regions", "sites"]:
        rows = []
        for _, row in bacteria.iterrows():
            for group in split_values(row[group_col]):
                rows.append({
                    "grouping": group_col,
                    "level": group,
                    "taxon_rank": "genus",
                    "taxon": row["gg2_genus"],
                    "derep_sequences": 1,
                    "source_asvs": row["n_asvs"],
                    "total_reads": row["total_reads"],
                    "prevalence_samples": row["prevalence_samples"],
                    "n_bioprojects_sum": row["n_bioprojects"],
                })
        if rows:
            group_parts.append(pd.DataFrame(rows))
    if group_parts:
        group_tax = pd.concat(group_parts, ignore_index=True).groupby(
            ["grouping", "level", "taxon_rank", "taxon"], dropna=False
        ).agg(
            derep_sequences=("derep_sequences", "sum"),
            source_asvs=("source_asvs", "sum"),
            total_reads=("total_reads", "sum"),
            prevalence_samples=("prevalence_samples", "sum"),
            n_bioprojects_sum=("n_bioprojects_sum", "sum"),
        ).reset_index().sort_values(["grouping", "level", "total_reads"], ascending=[True, True, False])
        group_tax.to_csv(group_tax_summary_path, sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(group_tax_summary_path, sep="\t", index=False)

    stats = {
        "candidate_sequences": int(len(cand)),
        "gg2_hits": int(ann["gg2_target"].notna().sum()),
        "qc_class_counts": qc.to_dict(orient="records"),
        "outputs": {
            "candidate_fasta": fasta_path,
            "gg2_hits": hit_path,
            "annotated_candidates": annotated_path,
            "qc_summary": qc_summary_path,
            "taxonomy_summary": tax_summary_path,
            "group_taxonomy_summary": group_tax_summary_path,
        },
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    log(f"complete: {len(cand):,} candidates; GG2 hits={stats['gg2_hits']:,}")


if __name__ == "__main__":
    main()
