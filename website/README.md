# Global Oral Microbiome Atlas — companion website

A static, self-contained data portal for the Global Oral Microbiome Resource. It is the
interactive companion to the resource manuscript: a study/cohort browser, a genus-level
**taxon explorer**, and an interactive view of the **geographic reference gap** in GG2
closed-reference mapping.

No server, database, or build toolchain is required at serve time — the front-end reads
pre-computed static JSON. Deploy the `site/` directory as-is to Cloudflare Pages, GitHub
Pages, Netlify, or any static host.

## Layout

```
website/
  precompute/
    build_data.py     # metadata TSVs + genus BIOM  ->  static JSON aggregates
  site/               # <-- deploy this directory
    index.html
    assets/{style.css, app.js}
    data/{overview,studies,taxa,refgap}.json   # generated; committed for convenience
```

## Regenerating the data

Requires `h5py` and `pandas` (only for reading the HDF5 BIOM; the rest is pure Python).
The genus table ships as a QIIME2 `.qza`; extract its inner BIOM first.

```bash
python3 -m venv .venv && .venv/bin/pip install h5py pandas

# extract the genus feature table from the .qza (a zip archive)
unzip -o data/metadata-and-biom-v2-qc/genus_table_relabund_all_qc.qza -d /tmp/genus
BIOM=$(find /tmp/genus -name feature-table.biom)

.venv/bin/python website/precompute/build_data.py \
    --meta-dir data/metadata-and-biom-v2-qc \
    --genus-biom "$BIOM" \
    --out website/site/data
```

Runtime is a few seconds. Outputs (compact JSON):

| File | Content | Approx size |
|---|---|---|
| `overview.json` | Global totals; counts by continent / site / amplicon / disease | 4 KB |
| `studies.json`  | One row per BioProject (geography, sites, amplicons, disease, N) | 120 KB |
| `taxa.json`     | 5,220 genera: prevalence & mean relabund, overall + by site/continent/status | 2 MB |
| `refgap.json`   | GG2 mapping rate by continent + amplicon-adjusted regression | 3 KB |

## Serving locally

Browsers block `fetch()` of local files over `file://`, so serve over HTTP:

```bash
cd website/site && python3 -m http.server 8000   # then open http://localhost:8000
```

## Notes

- **Heavy data belongs on a DOI archive, not here.** The multi-GB BIOM/metadata tables
  should be deposited on Zenodo/figshare; update the "Data & download" tab links and the
  DOI badge before publication. This portal hosts only the small aggregates and the
  interactive views.
- The taxon explorer aggregates the genus-level **relative-abundance** table over the
  73,340 QC-passed oral samples. Prevalence = fraction of samples in a group with abundance
  > 0; mean = mean relative abundance including zeros.
- The reference-gap regression numbers are the validated study-level weighted logit model
  from `results/geographic_reference_bias.md` §4.2; the per-continent descriptive table is
  recomputed from the pre-filter merged metadata at build time.
