#!/usr/bin/env python
"""Per-phylum ASV counts in 1%-identity bins (70-100) vs best HOMD identity."""
import pandas as pd

E = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13f_homd_phylum_ecdf.tsv"
C = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13d_homd_phylum_coverage.tsv"
OUT = "/beegfs/huangxiaochang/global_oral/global-oral-microbiome/HOMD_uncapture_analysis/pipeline_v3/13n_homd_phylum_identity_bins.tsv"

e = pd.read_csv(E, sep="\t")
c = pd.read_csv(C, sep="\t")
nmap = dict(zip(c.phylum, c.n_asvs_total))

levels = [float(x) for x in range(70, 101)]
rows = []
bins = [(70, 75), (75, 80), (80, 85), (85, 90), (90, 95), (95, 100)]
for ph, g in e.groupby("phylum"):
    n = nmap[ph]
    g = g[g.identity.isin(levels)].sort_values("identity")
    fg = dict(zip(g.identity, g.frac_ge))
    for (a, b) in bins:
        if b == 100:
            cnt = round(fg[a] * n)          # [a,100] includes identity == 100
        else:
            cnt = round((fg[a] - fg[b]) * n)
        rows.append((ph, a, b, (a + b) / 2.0, cnt))

out = pd.DataFrame(rows, columns=["phylum", "bin_start", "bin_end", "bin_mid", "n_asvs"])
out.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT}: {len(out)} rows")
