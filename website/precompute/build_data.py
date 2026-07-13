#!/usr/bin/env python3
"""Precompute static JSON aggregates for the Global Oral Microbiome Atlas website.

Reads the harmonized QC metadata (plain TSV) and the genus-level relative-abundance
BIOM (HDF5, extracted from the QIIME2 .qza) and writes compact JSON consumed by the
static front-end in website/site/. No server or database is required at serve time.

Outputs (website/site/data/):
  overview.json  - global totals + counts by continent / site / amplicon / disease
  studies.json   - one row per BioProject (geography, sites, amplicons, disease, N)
  taxa.json      - per-genus prevalence & mean relabund, overall + by site/continent/status
  refgap.json    - GG2 closed-reference mapping rate by continent + amplicon-adjusted model

Run:
  <venv>/bin/python website/precompute/build_data.py \
      --meta-dir data/metadata-and-biom-v2-qc \
      --genus-biom /path/to/feature-table.biom \
      --out website/site/data
"""
from __future__ import annotations
import argparse, csv, json, math, os, sys
from collections import defaultdict

# The six oral sites used throughout the resource (site_cleaned values).
ORAL_SITES = ["saliva", "tooth", "tongue", "oral mucosa", "oral swab", "gingiva"]
CONTINENTS = ["North America", "Asia", "Europe", "Africa", "South America", "Oceania"]

# Amplicon-adjusted continent effect on GG2 closed-reference mapping, study-level
# weighted logit regression (ref = North America), from
# results/geographic_reference_bias.md section 4.2. These are validated results;
# the per-continent descriptive summary below is recomputed from the raw table.
REFGAP_REGRESSION = [
    {"continent": "Asia",          "coef": -0.735, "ci": [-1.03, -0.44], "or": 0.48, "p": 1e-6,  "sig": True},
    {"continent": "Africa",        "coef": -1.166, "ci": [-1.82, -0.51], "or": 0.31, "p": 5e-4,  "sig": True},
    {"continent": "Oceania",       "coef": -0.617, "ci": [-1.43,  0.20], "or": 0.54, "p": 0.14,  "sig": False},
    {"continent": "South America", "coef": -0.267, "ci": [-1.06,  0.52], "or": 0.77, "p": 0.51,  "sig": False},
    {"continent": "Europe",        "coef": -0.215, "ci": [-0.54,  0.11], "or": 0.81, "p": 0.20,  "sig": False},
]

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def read_meta(path):
    """Yield metadata rows as dicts from a TSV."""
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def clean(v):
    if v is None:
        return ""
    v = v.strip()
    return "" if v.lower() in ("", "na", "nan", "none", "null") else v


def status_group(status_cleaned):
    """Map status_cleaned -> ('healthy'|'disease'|'', disease_name_or_None).

    status_cleaned is 'HC' for healthy controls, 'NA'/'' for unknown, and the
    disease name itself (e.g. 'Dental Caries') for cases.
    """
    s = clean(status_cleaned)
    if s == "":
        return "", None
    if s.upper() == "HC":
        return "healthy", None
    return "disease", s


def round_sig(x, n=4):
    if x == 0 or x is None:
        return 0
    return round(x, n)


# ---------------------------------------------------------------------------
# 1. Metadata-driven aggregates (overview + studies), from the oral QC table.
# ---------------------------------------------------------------------------
def build_meta_aggregates(oral_tsv):
    studies = defaultdict(lambda: {
        "n_samples": 0, "subjects": set(), "sites": defaultdict(int),
        "continents": set(), "amplicons": set(), "countries": set(),
        "diseases": set(), "n_case": 0, "n_control": 0,
        "platform": "", "public_time": "", "title": "", "disease_study": False,
    })
    by_continent = defaultdict(int)
    by_site = defaultdict(int)
    by_amplicon = defaultdict(int)
    by_country = defaultdict(int)
    disease_counts = defaultdict(int)          # sample-level case counts by disease label
    disease_projects = defaultdict(set)
    total_samples = 0
    all_subjects = set()

    for r in read_meta(oral_tsv):
        sid = clean(r.get("sample_id"))
        bp = clean(r.get("bioproject_id"))
        if not sid or not bp:
            continue
        site = clean(r.get("site_cleaned")).lower()
        cont = clean(r.get("continent_group"))
        amp = clean(r.get("sequencing_region_cleaned"))
        country = clean(r.get("country_cleaned"))
        subj = clean(r.get("subject_id"))
        sgroup, disease_lbl = status_group(r.get("status_cleaned"))
        is_disease_study = clean(r.get("disease_study_cleaned")).lower() == "yes"

        total_samples += 1
        if cont:
            by_continent[cont] += 1
        if site:
            by_site[site] += 1
        if amp:
            by_amplicon[amp] += 1
        if country:
            by_country[country] += 1
        if subj:
            all_subjects.add((bp, subj))

        s = studies[bp]
        s["n_samples"] += 1
        if subj:
            s["subjects"].add(subj)
        if site:
            s["sites"][site] += 1
        if cont:
            s["continents"].add(cont)
        if amp:
            s["amplicons"].add(amp)
        if country:
            s["countries"].add(country)
        if clean(r.get("platform_cleaned")):
            s["platform"] = clean(r.get("platform_cleaned"))
        if clean(r.get("public_time")):
            s["public_time"] = clean(r.get("public_time"))
        if clean(r.get("title")) and not s["title"]:
            s["title"] = clean(r.get("title"))[:300]
        if is_disease_study:
            s["disease_study"] = True
        if sgroup == "disease":
            s["n_case"] += 1
            s["diseases"].add(disease_lbl)
            disease_counts[disease_lbl] += 1
            disease_projects[disease_lbl].add(bp)
        elif sgroup == "healthy":
            s["n_control"] += 1

    # Serialize study rows.
    study_rows = []
    for bp, s in sorted(studies.items(), key=lambda kv: -kv[1]["n_samples"]):
        study_rows.append({
            "bioproject": bp,
            "n_samples": s["n_samples"],
            "n_subjects": len(s["subjects"]),
            "sites": sorted(s["sites"].keys()),
            "top_site": max(s["sites"].items(), key=lambda kv: kv[1])[0] if s["sites"] else "",
            "continents": sorted(s["continents"]),
            "amplicons": sorted(s["amplicons"]),
            "countries": sorted(s["countries"]),
            "diseases": sorted(s["diseases"]),
            "is_disease_study": s["disease_study"],
            "n_case": s["n_case"],
            "n_control": s["n_control"],
            "platform": s["platform"],
            "public_time": s["public_time"],
            "title": s["title"],
        })

    overview = {
        "totals": {
            "samples": total_samples,
            "studies": len(studies),
            "subjects": len(all_subjects),
            "sites": len([k for k in by_site if k]),
            "continents": len([k for k in by_continent if k]),
            "genera": None,  # filled in by taxa step
        },
        "by_continent": dict(sorted(by_continent.items(), key=lambda kv: -kv[1])),
        "by_site": dict(sorted(by_site.items(), key=lambda kv: -kv[1])),
        "by_amplicon": dict(sorted(by_amplicon.items(), key=lambda kv: -kv[1])),
        "by_country": dict(sorted(by_country.items(), key=lambda kv: -kv[1])),
        "disease_case_counts": dict(sorted(disease_counts.items(), key=lambda kv: -kv[1])),
        "disease_n_projects": {k: len(v) for k, v in sorted(
            disease_projects.items(), key=lambda kv: -len(kv[1]))},
    }
    return overview, study_rows


# ---------------------------------------------------------------------------
# 2. Taxon aggregates from the genus relative-abundance BIOM (HDF5).
# ---------------------------------------------------------------------------
def build_taxa(genus_biom, oral_tsv):
    import h5py
    import numpy as np

    # sample_id -> (site, continent, status_group) for oral samples only.
    sample_meta = {}
    for r in read_meta(oral_tsv):
        sid = clean(r.get("sample_id"))
        if not sid:
            continue
        site = clean(r.get("site_cleaned")).lower()
        cont = clean(r.get("continent_group"))
        sg, _ = status_group(r.get("status_cleaned"))
        sample_meta[sid] = (site if site in ORAL_SITES else "", cont, sg)

    f = h5py.File(genus_biom, "r")
    obs_ids = [x.decode() if isinstance(x, bytes) else x for x in f["observation/ids"][:]]
    samp_ids = [x.decode() if isinstance(x, bytes) else x for x in f["sample/ids"][:]]
    n_obs = len(obs_ids)
    n_samp = len(samp_ids)

    # Column (sample) -> group indices. -1 means "skip / not oral".
    col_site = np.full(n_samp, -1, dtype=np.int8)
    col_cont = np.full(n_samp, -1, dtype=np.int8)
    col_status = np.full(n_samp, -1, dtype=np.int8)
    site_idx = {s: i for i, s in enumerate(ORAL_SITES)}
    cont_idx = {c: i for i, c in enumerate(CONTINENTS)}
    status_idx = {"healthy": 0, "disease": 1}
    oral_mask = np.zeros(n_samp, dtype=bool)
    for j, sid in enumerate(samp_ids):
        m = sample_meta.get(sid)
        if not m:
            continue
        site, cont, sg = m
        if site == "" and cont == "":
            continue
        oral_mask[j] = True
        if site in site_idx:
            col_site[j] = site_idx[site]
        if cont in cont_idx:
            col_cont[j] = cont_idx[cont]
        if sg in status_idx:
            col_status[j] = status_idx[sg]

    n_oral = int(oral_mask.sum())
    # Denominators: number of oral samples per group (for prevalence & mean incl. zeros).
    site_totals = np.array([(col_site == i).sum() for i in range(len(ORAL_SITES))], dtype=np.int64)
    cont_totals = np.array([(col_cont == i).sum() for i in range(len(CONTINENTS))], dtype=np.int64)
    status_totals = np.array([(col_status == i).sum() for i in range(2)], dtype=np.int64)

    # observation/matrix is CSR: for each genus (row), the nonzero sample columns.
    data = f["observation/matrix/data"]
    indices = f["observation/matrix/indices"]
    indptr = f["observation/matrix/indptr"][:]

    taxa = []
    for i in range(n_obs):
        a, b = indptr[i], indptr[i + 1]
        if a == b:
            continue
        cols = indices[a:b][:]
        vals = data[a:b][:]
        # Restrict to oral samples.
        keep = oral_mask[cols]
        if not keep.any():
            continue
        cols = cols[keep]
        vals = vals[keep].astype(np.float64)

        overall_sum = float(vals.sum())
        overall_present = int((vals > 0).sum())

        def group_stats(col_group, totals):
            g = col_group[cols]
            k = totals.shape[0]
            sums = np.zeros(k); pres = np.zeros(k, dtype=np.int64)
            valid = g >= 0
            np.add.at(sums, g[valid], vals[valid])
            np.add.at(pres, g[valid], (vals[valid] > 0).astype(np.int64))
            out = {}
            for gi in range(k):
                if totals[gi] > 0 and pres[gi] > 0:
                    out[gi] = [round_sig(sums[gi] / totals[gi], 5),      # mean relabund (incl zeros)
                               round_sig(pres[gi] / totals[gi], 4)]      # prevalence
            return out

        # Parse taxonomy string into levels.
        tax = obs_ids[i]
        levels = [p.strip() for p in tax.split(";")]
        def lvl(prefix):
            for p in levels:
                if p.startswith(prefix):
                    val = p[len(prefix):]
                    return val if val else ""
            return ""
        genus = lvl("g__") or lvl("f__") or lvl("o__") or "(unclassified)"

        taxa.append({
            "id": i,
            "genus": genus,
            "phylum": lvl("p__"),
            "family": lvl("f__"),
            "tax": tax,
            "mean": round_sig(overall_sum / n_oral, 6),
            "prev": round_sig(overall_present / n_oral, 4),
            "by_site": group_stats(col_site, site_totals),
            "by_continent": group_stats(col_cont, cont_totals),
            "by_status": group_stats(col_status, status_totals),
        })

    taxa.sort(key=lambda t: -t["prev"])
    meta = {
        "n_oral_samples": n_oral,
        "sites": ORAL_SITES,
        "continents": CONTINENTS,
        "status": ["healthy", "disease"],
        "site_totals": site_totals.tolist(),
        "continent_totals": cont_totals.tolist(),
        "status_totals": status_totals.tolist(),
    }
    return {"meta": meta, "taxa": taxa}, len(taxa)


# ---------------------------------------------------------------------------
# 3. Reference-gap summary, recomputed from the raw (pre-filter) merged table.
# ---------------------------------------------------------------------------
def build_refgap(raw_tsv):
    per = defaultdict(lambda: {"rates": [], "studies": defaultdict(list)})
    for r in read_meta(raw_tsv):
        site = clean(r.get("site_cleaned")).lower()
        if site not in ORAL_SITES:
            continue
        cont = clean(r.get("continent_group"))
        rate = clean(r.get("sample_cr_mapping_rate"))
        bp = clean(r.get("bioproject_id"))
        amp = clean(r.get("sequencing_region_cleaned"))
        if not cont or not rate:
            continue
        try:
            rate = float(rate)
        except ValueError:
            continue
        for subset in ("all", "v4"):
            if subset == "v4" and amp != "V4":
                continue
            key = (subset, cont)
            per[key]["rates"].append(rate)
            per[key]["studies"][bp].append(rate)

    def median(xs):
        xs = sorted(xs); n = len(xs)
        if n == 0:
            return None
        return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2

    rows = []
    for (subset, cont), d in per.items():
        rates = d["rates"]
        if not rates:
            continue
        below = sum(1 for x in rates if x < 0.30) / len(rates) * 100
        rows.append({
            "subset": "V4 only" if subset == "v4" else "All amplicons",
            "continent": cont,
            "n_samples": len(rates),
            "n_studies": len(d["studies"]),
            "median_mapping": round_sig(median(rates), 3),
            "mean_unmapped_pct": round_sig((1 - sum(rates) / len(rates)) * 100, 2),
            "pct_below_floor": round_sig(below, 2),
        })
    order = {c: i for i, c in enumerate(CONTINENTS)}
    rows.sort(key=lambda r: (0 if r["subset"] == "All amplicons" else 1,
                             order.get(r["continent"], 99)))
    return {"summary": rows, "regression": REFGAP_REGRESSION,
            "regression_ref": "North America", "floor": 0.30}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-dir", default="data/metadata-and-biom-v2-qc")
    ap.add_argument("--genus-biom", required=True)
    ap.add_argument("--out", default="website/site/data")
    args = ap.parse_args()

    oral_tsv = os.path.join(args.meta_dir, "metadata_oral_qc.tsv")
    raw_tsv = os.path.join(args.meta_dir, "metadata_all_for_qc_raw.tsv")
    os.makedirs(args.out, exist_ok=True)

    print("[1/4] metadata aggregates (overview + studies) ...", flush=True)
    overview, studies = build_meta_aggregates(oral_tsv)

    print("[2/4] taxon aggregates from genus BIOM ...", flush=True)
    taxa, n_taxa = build_taxa(args.genus_biom, oral_tsv)
    overview["totals"]["genera"] = n_taxa

    print("[3/4] reference-gap summary ...", flush=True)
    refgap = build_refgap(raw_tsv)

    print("[4/4] writing JSON ...", flush=True)
    for name, obj in [("overview", overview), ("studies", {"studies": studies}),
                      ("taxa", taxa), ("refgap", refgap)]:
        p = os.path.join(args.out, f"{name}.json")
        with open(p, "w") as fh:
            json.dump(obj, fh, separators=(",", ":"), ensure_ascii=False)
        print(f"   wrote {p}  ({os.path.getsize(p)/1e6:.2f} MB)")

    print("done. totals:", overview["totals"])


if __name__ == "__main__":
    main()
