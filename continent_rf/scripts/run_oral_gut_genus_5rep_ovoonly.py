#!/usr/bin/env python3
"""Genus-level 5-repeat OVO-only mirror of the oral/gut continent RF.

Targets the two-panel upper-triangle OVO figure: it needs Saliva+Gut OVO
only, so OVA (the expensive one-vs-all pass for all 6 continents) is skipped.
CV, RF settings and fold scheme match run_oral_gut_project_cv_repeat5.py;
features are aggregated to GENUS (taxonomy.tsv) then renormalized to unit sum.

Outputs (original result/ and result_genus_1rep/ are untouched):
  result_genus_5rep/oral/  region_ovo_{summary,fold_metrics,predictions}.tsv
  result_genus_5rep/gut/   same
At the end writes result_genus_5rep/comparison_genus_vs_feature.csv with
original (feature, 5-repeat) vs genus (5-repeat) OVO mean AUCs side by side.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_oral_gut_project_cv_repeat5 as rr

rr.N_REPEATS = 5
OUTBASE = Path(__file__).resolve().parents[1] / "result_genus_5rep"
OUTBASE.mkdir(parents=True, exist_ok=True)

TAXONOMY_TSV = Path("/beegfs/huangxiaochang/global_oral/taxonomy.tsv")


def parse_genus(taxon_str):
    """Extract the genus label (g__...) from a semicolon-separated taxon."""
    if not isinstance(taxon_str, str):
        return None
    for part in taxon_str.split(";"):
        part = part.strip()
        if part.startswith("g__"):
            label = part[len("g__"):].strip()
            return label or None
    return None


def load_genus_map():
    import pandas as pd
    tax = pd.read_csv(TAXONOMY_TSV, sep="\t")
    tax["genus"] = tax["Taxon"].apply(parse_genus)
    tax = tax.dropna(subset=["genus"])
    return dict(zip(tax["Feature ID"].astype(str), tax["genus"]))


def main():
    import biom
    from scipy import sparse
    import numpy as np
    import pandas as pd

    def ids_from_table(path, wanted, norm=True, empty=True):
        table = biom.load_table(str(path))
        shared = [sid for sid in table.ids(axis="sample") if str(sid) in wanted]
        if not shared:
            raise ValueError(f"No sample overlap for {path}")
        table = table.filter(shared, axis="sample", inplace=False)
        if empty:
            table.remove_empty(inplace=True)
        if norm:
            table.norm(axis="sample", inplace=False)
        return table

    def to_genus(table, genus_map, tag):
        obs = np.asarray(table.ids(axis="observation"), dtype=str)
        col_genus = np.array([genus_map.get(f) for f in obs], dtype=object)
        known = col_genus != None  # noqa: E711
        print(f"[{tag}] features kept->genus: {known.sum():,}/{len(obs):,} "
              f"({len(set(col_genus[known])):,} genera)", flush=True)
        sub = table.matrix_data.T.tocsr()[:, known]
        genera = np.array(sorted(set(col_genus[known])), dtype=object)
        gcol = {g: k for k, g in enumerate(genera)}
        mapping = np.fromiter((gcol[col_genus[i]] for i in np.nonzero(known)[0]),
                              dtype=np.int64, count=known.sum())
        mat = sub @ sparse.csr_matrix(
            (np.ones(len(mapping)), (np.arange(len(mapping)), mapping)),
            shape=(len(mapping), len(genera)))
        # renormalize each sample to unit sum
        row_sum = np.asarray(mat.sum(axis=1)).ravel()
        nz = row_sum > 0
        s = np.ones_like(row_sum); s[nz] = row_sum[nz]
        from scipy.sparse import diags
        mat = diags(1.0 / s) @ mat
        return mat.tocsr(), genera

    # ---- oral ----
    meta_o, _mat = rr.load_oral()
    table = ids_from_table(rr.ORAL_BIOM, set(meta_o["sample_id"].astype(str)))
    genus_map = load_genus_map()
    mat_o, genera_o = to_genus(table, genus_map, "oral")
    sids = np.asarray(table.ids(axis="sample"), dtype=str)
    assert np.array_equal(sids, meta_o["sample_id"].to_numpy()), \
        "oral sample order drift between loaders"
    print(f"Oral genus matrix: {mat_o.shape}", flush=True)

    # ---- gut ----
    meta_g, _matg = rr.load_gut()
    table_g = ids_from_table(rr.GUT_BIOM, set(meta_g["sample_id"].astype(str)))
    mat_g, genera_g = to_genus(table_g, genus_map, "gut")
    sids_g = np.asarray(table_g.ids(axis="sample"), dtype=str)
    assert np.array_equal(sids_g, meta_g["sample_id"].to_numpy()), \
        "gut sample order drift between loaders"
    print(f"Gut genus matrix: {mat_g.shape}", flush=True)

    # ---- OVO only (OVA skipped: not needed for the two-panel figure) ----
    outdirs = {"oral": OUTBASE / "oral", "gut": OUTBASE / "gut"}
    ovo_summaries = []
    for dataset, (meta, matrix) in (
        ("oral", (meta_o, mat_o)),
        ("gut", (meta_g, mat_g)),
    ):
        out = outdirs[dataset]
        out.mkdir(parents=True, exist_ok=True)
        pred_o = out / "region_ovo_predictions.tsv"
        if pred_o.exists():
            pred_o.unlink()
        folds_o, sums_o = rr.run_ovo(meta, matrix, dataset, pred_o)
        pd.DataFrame(folds_o).to_csv(out / "region_ovo_fold_metrics.tsv",
                                     sep="\t", index=False)
        pd.DataFrame(sums_o).to_csv(out / "region_ovo_summary.tsv",
                                    sep="\t", index=False)
        ovo_summaries += sums_o

    # ---- side-by-side comparison with the original 5-repeat results -------
    # result/oral/region_ovo_summary.tsv already contains BOTH datasets'
    # rows (oral/gut files are full merged copies); use one and normalize.
    orig_ovo = pd.read_csv(rr.ROOT / "continent_rf" / "result" /
                           "oral/region_ovo_summary.tsv", sep="\t")
    orig_ovo["dataset"] = orig_ovo["dataset"].str.lower()
    orig_key = orig_ovo.assign(
        site=orig_ovo["site"].astype(str),
        pair=orig_ovo["region_1"] + " vs " + orig_ovo["region_2"])
    genus_key = pd.DataFrame(ovo_summaries).assign(
        dataset=lambda d: d["dataset"].str.lower(),
        site=lambda d: d["site"].astype(str),
        pair=lambda d: d["region_1"] + " vs " + d["region_2"])
    comp = (orig_key[["dataset", "site", "pair", "mean_auc"]]
            .merge(genus_key[["dataset", "site", "pair", "mean_auc"]],
                   on=["dataset", "site", "pair"],
                   suffixes=("_feature_5rep", "_genus_5rep")))
    comp["delta_auc"] = comp["mean_auc_genus_5rep"] - comp["mean_auc_feature_5rep"]
    comp = comp.sort_values("delta_auc")
    comp.to_csv(OUTBASE / "comparison_genus_vs_feature.csv", index=False)
    print("\n===== OVO mean AUC: genus(5rep) vs feature(5rep) =====")
    print(comp.to_string(index=False), flush=True)

    from datetime import datetime
    (OUTBASE / "run_config.json").write_text(
        __import__("json").dumps({
            "script": "run_oral_gut_genus_5rep_ovoonly.py",
            "generated": datetime.now().isoformat(),
            "repeats": 5, "trees": rr.N_TREES, "seed": rr.SEED,
            "scope": "OVO only (OVA skipped — two-panel OVO figure needs none)",
            "level": "genus (taxonomy.tsv, g__ aggregation, unit-sum renorm)",
            "taxonomy": str(TAXONOMY_TSV),
        }, indent=2, sort_keys=True) + "\n")
    print(f"[done] outputs in {OUTBASE}", flush=True)


if __name__ == "__main__":
    main()