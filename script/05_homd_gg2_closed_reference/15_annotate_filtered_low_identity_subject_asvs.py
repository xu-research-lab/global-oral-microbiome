#!/usr/bin/env python3
"""Annotate strictly filtered low-HOMD-identity ASVs with GG2 and plot bacterial GG2 identity distribution."""

from __future__ import annotations

import argparse
import json
import os
import subprocess

import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_OUTDIR = "/beegfs/huangxiaochang/global_oral/ASV_vs_HOMD_analysis/region_aware_pipeline"
DEFAULT_VSEARCH = "/home/xuhuan/miniconda3/envs/vsearch_env/bin/vsearch"
DEFAULT_GG2_REF = "/beegfs/dongbiao/greengene2/2024.09.seqs.fna"
DEFAULT_GG2_TAX = "/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv"

sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150, "savefig.dpi": 300})


def log(msg):
    print(f"[filtered-gg2] {msg}", flush=True)


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
        next(handle, None)
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
            val = part.split("__", 1)[1].strip()
            return val if val else "Unassigned"
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--vsearch", default=DEFAULT_VSEARCH)
    parser.add_argument("--gg2-ref", default=DEFAULT_GG2_REF)
    parser.add_argument("--gg2-tax", default=DEFAULT_GG2_TAX)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    in_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered.tsv")
    fasta_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered.fa")
    hit_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_vs_gg2_best_hits.tsv")
    annotated_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_annotated.tsv")
    bacterial_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_only.tsv")
    qc_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_qc_summary.tsv")
    tax_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_bacterial_taxonomy_summary.tsv")
    stats_path = os.path.join(args.outdir, "low_identity_lt97_subject_prevalence_filtered_gg2_stats.json")
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig_base = os.path.join(figdir, "fig9_filtered_low_identity_bacterial_gg2_identity_distribution")

    log("loading filtered ASVs")
    df = pd.read_csv(in_path, sep="\t")
    write_fasta(df, fasta_path)

    if args.force or not (os.path.exists(hit_path) and os.path.getsize(hit_path) > 0):
        cmd = [args.vsearch, "--usearch_global", fasta_path, "--db", args.gg2_ref, "--id", "0.7", "--maxaccepts", "1", "--maxrejects", "0", "--strand", "both", "--userout", hit_path, "--userfields", "query+target+id+alnlen+qcov+tcov+ql+tl", "--threads", str(args.threads), "--quiet"]
        log("running GG2 alignment")
        log(" ".join(cmd))
        subprocess.run(cmd, check=True)
    else:
        log(f"using existing GG2 hits: {hit_path}")

    hit_cols = ["derep_id", "gg2_target", "gg2_identity", "gg2_alnlen", "gg2_qcov", "gg2_tcov", "query_length", "target_length"]
    hits = pd.read_csv(hit_path, sep="\t", names=hit_cols) if os.path.getsize(hit_path) else pd.DataFrame(columns=hit_cols)
    tax = load_taxonomy(args.gg2_tax)
    if len(hits):
        hits["gg2_taxonomy"] = hits["gg2_target"].map(tax).fillna("Unassigned")
        hits["has_taxonomy"] = hits["gg2_taxonomy"].ne("Unassigned")

        # GG2 FASTA contains some IDs absent from the exported taxonomy table.
        # Keep the absolute best alignment for audit, but use the best taxonomy-resolved
        # hit for QC filtering so annotated bacterial ASVs are not discarded incorrectly.
        sort_cols = ["derep_id", "gg2_identity", "gg2_qcov", "gg2_alnlen"]
        best_overall = (
            hits.sort_values(sort_cols, ascending=[True, False, False, False])
            .drop_duplicates("derep_id")
            .rename(columns={
                "gg2_target": "gg2_best_target",
                "gg2_identity": "gg2_best_identity",
                "gg2_alnlen": "gg2_best_alnlen",
                "gg2_qcov": "gg2_best_qcov",
                "gg2_tcov": "gg2_best_tcov",
                "query_length": "gg2_best_query_length",
                "target_length": "gg2_best_target_length",
                "gg2_taxonomy": "gg2_best_taxonomy",
                "has_taxonomy": "gg2_best_has_taxonomy",
            })
        )
        best_tax = (
            hits[hits["has_taxonomy"]]
            .sort_values(sort_cols, ascending=[True, False, False, False])
            .drop_duplicates("derep_id")
        )
        hits = best_overall.merge(
            best_tax[["derep_id", "gg2_target", "gg2_identity", "gg2_alnlen", "gg2_qcov", "gg2_tcov", "query_length", "target_length", "gg2_taxonomy"]],
            on="derep_id",
            how="left",
        )
    else:
        hits["gg2_taxonomy"] = []

    ann = df.merge(hits, on="derep_id", how="left")
    ann["gg2_taxonomy"] = ann["gg2_taxonomy"].fillna("Unassigned")
    ann["gg2_qc_class"] = ann["gg2_taxonomy"].map(classify_qc)
    for col, prefix in [("gg2_domain", "d__"), ("gg2_phylum", "p__"), ("gg2_class", "c__"), ("gg2_order", "o__"), ("gg2_family", "f__"), ("gg2_genus", "g__"), ("gg2_species", "s__")]:
        ann[col] = ann["gg2_taxonomy"].map(lambda x, p=prefix: parse_rank(x, p))
    ann.to_csv(annotated_path, sep="\t", index=False)

    bacterial = ann[ann["gg2_qc_class"].eq("bacterial")].copy()
    bacterial.to_csv(bacterial_path, sep="\t", index=False)

    qc = ann.groupby("gg2_qc_class", dropna=False).agg(
        derep_sequences=("derep_id", "count"),
        source_asvs=("n_asvs", "sum"),
        total_reads=("total_reads", "sum"),
        subject_prevalence=("subject_prevalence", "sum"),
        max_subject_prevalence=("subject_prevalence", "max"),
        max_bioprojects=("n_bioprojects", "max"),
    ).reset_index().sort_values("subject_prevalence", ascending=False)
    qc.to_csv(qc_path, sep="\t", index=False)

    tax_rows = []
    for rank in ["gg2_phylum", "gg2_class", "gg2_order", "gg2_family", "gg2_genus", "gg2_species"]:
        tmp = bacterial.groupby(rank, dropna=False).agg(
            derep_sequences=("derep_id", "count"),
            source_asvs=("n_asvs", "sum"),
            total_reads=("total_reads", "sum"),
            subject_prevalence=("subject_prevalence", "sum"),
            max_subject_prevalence=("subject_prevalence", "max"),
            max_bioprojects=("n_bioprojects", "max"),
        ).reset_index().rename(columns={rank: "taxon"})
        tmp.insert(0, "rank", rank.replace("gg2_", ""))
        tax_rows.append(tmp)
    pd.concat(tax_rows, ignore_index=True).sort_values(["rank", "subject_prevalence"], ascending=[True, False]).to_csv(tax_path, sep="\t", index=False)

    log("plotting GG2 identity distribution")
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    sns.histplot(data=bacterial, x="gg2_identity", bins=30, color="#356A8A", edgecolor="white", ax=ax)
    ax.axvline(97, color="#B22222", linestyle="--", linewidth=1.2, label="97%")
    ax.axvline(99, color="#333333", linestyle=":", linewidth=1.2, label="99%")
    ax.set_xlabel("Best Greengenes2 identity (%)")
    ax.set_ylabel("Bacterial ASV count")
    ax.set_title("GG2 identity distribution of filtered low-HOMD-identity bacterial ASVs")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_base + ".pdf", bbox_inches="tight")
    fig.savefig(fig_base + ".png", bbox_inches="tight", dpi=300)
    plt.close(fig)

    stats = {
        "filtered_input_asvs": int(len(df)),
        "gg2_hits": int(ann["gg2_best_target"].notna().sum()) if "gg2_best_target" in ann.columns else int(ann["gg2_target"].notna().sum()),
        "gg2_taxonomy_resolved_hits": int(ann["gg2_target"].notna().sum()),
        "bacterial_retained_asvs": int(len(bacterial)),
        "qc_summary": qc.to_dict(orient="records"),
        "outputs": {
            "annotated_all": annotated_path,
            "bacterial_only": bacterial_path,
            "qc_summary": qc_path,
            "taxonomy_summary": tax_path,
            "identity_distribution_pdf": fig_base + ".pdf",
            "identity_distribution_png": fig_base + ".png",
        },
    }
    with open(stats_path, "w") as handle:
        json.dump(stats, handle, indent=2)
    log(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
