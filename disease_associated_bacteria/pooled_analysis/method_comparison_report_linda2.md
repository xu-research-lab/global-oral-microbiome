# da_v2 genus-level differential abundance: MaAsLin3 / ANCOM-BC2 / LinDA2

**Working directory**: `/home/zhangcaiyi/project/global_oral_diseae_da/da_v2`  
**Generated**: 2026-07-21  

---

## 1. Analysis overview

This report reuses the existing **da_v2** MaAsLin3 and ANCOM-BC2 outputs, and reruns LinDA as **LinDA2** with `is.winsor = FALSE` into `linda2/`.

| Method | Design |
|--------|--------|
| **MaAsLin3** | Disease-level mixed model, `bioproject_id` random effect |
| **ANCOM-BC2** | Disease-level mixed model, `rand_formula = (1 \| bioproject_id)` |
| **LinDA2** | Disease-level mixed model with `(1\|bioproject_id)`, but `is.winsor = FALSE` |

**Rationale for LinDA2**: the original `linda/` run used `is.winsor = TRUE`, which introduced zero-library samples after winsorization in OSCC / Periodontal Disease / IBD and caused `NA/NaN/Inf in 'y'`. The new `linda2/` rerun disables winsorization and leaves the original `linda/` outputs unchanged.

## 2. MaAsLin3 run summary

Key call ([`run_maaslin3_disease.R`](../run_maaslin3_disease.R)):

```r
maaslin3(
  input_data = t(sub_abund), input_metadata = sub_meta,
  output = disease_dir, fixed_effects = fixed_effects,
  random_effects = c("bioproject_id")
)
```

| Disease | Cases | Controls | Matched_BioProjects | Fixed_Effects | Significant_q0.05 |
| --- | --- | --- | --- | --- | --- |
| Breast Cancer |  493 | 1768 |  2 | Disease_Status,sex_cleaned,sequencing_region_cleaned,continent_group |  78 |
| Dental Caries |  390 |  512 |  9 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,sequencing_region_cleaned,continent_group |   8 |
| Head and Neck Cancer |  763 |  885 |  5 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,continent_group | 114 |
| HIV |  303 |  142 |  6 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,sequencing_region_cleaned,continent_group |  36 |
| IBD |  399 | 2088 |  3 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,sequencing_region_cleaned,continent_group |  10 |
| Lung Cancer | 1316 | 3408 |  2 | Disease_Status,age_cleaned_group,sex_cleaned,sequencing_region_cleaned,continent_group |  72 |
| OSCC |  403 |  431 | 10 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,sequencing_region_cleaned,continent_group |  51 |
| Periodontal Disease | 1063 | 1138 | 18 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,sequencing_region_cleaned,continent_group |  99 |
| RA |  208 |  229 |  3 | Disease_Status,age_cleaned_group,sex_cleaned,site_cleaned,sequencing_region_cleaned,continent_group |  10 |

![MaAsLin3 UpSet](figures_linda2/maaslin3_upset.png)

## 3. ANCOM-BC2 run summary

Key call ([`run_ancombc2_disease.R`](../run_ancombc2_disease.R)):

```r
ancombc2(
  data = ps,
  fix_formula = ff,
  rand_formula = "(1 | bioproject_id)",
  p_adj_method = "BH",
  prv_cut = PRV_CUT
)
```

| Disease | Cases | Controls | Matched_BioProjects | Fixed_Effects | Significant_q0.05 |
| --- | --- | --- | --- | --- | --- |
| Breast Cancer |  493 | 1768 |  2 | Disease_Status | 10 |
| Dental Caries |  390 |  512 |  9 | Disease_Status |  0 |
| Head and Neck Cancer |  763 |  885 |  5 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + continent_group | 34 |
| HIV |  303 |  142 |  6 | Disease_Status |  1 |
| IBD |  396 | 2088 |  3 | Disease_Status |  4 |
| Lung Cancer | 1316 | 3408 |  2 | Disease_Status | 37 |
| OSCC |  403 |  431 | 10 | Disease_Status | 25 |
| Periodontal Disease | 1063 | 1137 | 18 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 62 |
| RA |  208 |  229 |  3 | Disease_Status |  4 |

![ANCOM-BC2 UpSet](figures_linda2/ancombc2_upset.png)

## 4. LinDA2 run summary

Key call ([`run_linda2_disease.R`](../run_linda2_disease.R)):

```r
linda(
  feature.dat = feature_mat,
  meta.dat = meta_df,
  formula = paste0("~", formula_str),
  feature.dat.type = "count",
  prev.filter = PRV_CUT,
  is.winsor = FALSE,
  outlier.pct = 0.03,
  p.adj.method = "BH"
)
```

| Disease | Cases | Controls | Matched_BioProjects | Fixed_Effects | Significant_q0.05 |
| --- | --- | --- | --- | --- | --- |
| Breast Cancer |  493 | 1768 |  2 | Disease_Status + sex_cleaned + sequencing_region_cleaned + continent_group | 54 |
| Dental Caries |  390 |  512 |  9 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 10 |
| Head and Neck Cancer |  763 |  885 |  5 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + continent_group | 67 |
| HIV |  303 |  142 |  6 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 31 |
| IBD |  396 | 2088 |  3 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 12 |
| Lung Cancer | 1316 | 3408 |  2 | Disease_Status + age_cleaned_group + sex_cleaned + sequencing_region_cleaned + continent_group | 64 |
| OSCC |  403 |  431 | 10 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 55 |
| Periodontal Disease | 1063 | 1137 | 18 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 77 |
| RA |  208 |  229 |  3 | Disease_Status + age_cleaned_group + sex_cleaned + site_cleaned + sequencing_region_cleaned + continent_group | 13 |

![LinDA2 UpSet](figures_linda2/linda2_upset.png)

## 5. Three-method comparison (Genus + Direction)

![Three-way Venn](figures_linda2/venn_3way_per_disease.png)

| Disease | n_maaslin3 | n_ancombc2 | n_linda | n_all3 | n_only_maaslin3 | n_only_ancombc2 | n_only_linda |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Periodontal Disease | 56 | 55 | 67 | 46 |  1 | 4 | 13 |
| Lung Cancer | 53 | 35 | 57 | 31 |  9 | 2 | 13 |
| Head and Neck Cancer | 74 | 31 | 64 | 24 | 34 | 2 | 23 |
| OSCC | 43 | 24 | 53 | 20 |  3 | 1 | 10 |
| Breast Cancer | 57 | 10 | 47 |  8 | 16 | 0 |  6 |
| RA | 10 |  4 | 13 |  4 |  1 | 0 |  4 |
| Dental Caries |  8 |  0 | 10 |  0 |  0 | 0 |  2 |
| HIV | 31 |  1 | 27 |  0 | 13 | 1 |  9 |
| IBD |  9 |  4 | 12 |  0 |  2 | 3 |  6 |

### Direction definition

| Direction | Meaning |
|-----------|----------|
| **`+`** | Enriched in **Disease (Case)** relative to Healthy Control |
| **`-`** | Enriched in **Healthy Control** |

### Genera significant in all three methods with the same Direction (all3)

| Disease | Genus | Direction |
| --- | --- | --- |
| Periodontal Disease | Abiotrophia | - |
| Periodontal Disease | Actinomyces | - |
| Periodontal Disease | AEES01000073 | - |
| Periodontal Disease | Alloprevotella | + |
| Periodontal Disease | Anaeroglobus | + |
| Periodontal Disease | Arachnia | - |
| Periodontal Disease | Bulleidia | + |
| Periodontal Disease | Capnocytophaga | - |
| Periodontal Disease | Cardiobacterium | - |
| Periodontal Disease | Corynebacterium | - |
| Periodontal Disease | Desulfobulbus | + |
| Periodontal Disease | Dialister | + |
| Periodontal Disease | Eggerthia | + |
| Periodontal Disease | Eubacterium | + |
| Periodontal Disease | Filifactor | + |
| Periodontal Disease | Flexilinea | + |
| Periodontal Disease | Fretibacterium | + |
| Periodontal Disease | Gallibacter | + |
| Periodontal Disease | Gemella | - |
| Periodontal Disease | Granulicatella | - |
| Periodontal Disease | Haemophilus | - |
| Periodontal Disease | Hornefia | + |
| Periodontal Disease | Kingella | - |
| Periodontal Disease | Lancefieldella | + |
| Periodontal Disease | Lautropia | - |
| Periodontal Disease | Leptotrichia | - |
| Periodontal Disease | Metamycoplasma | + |
| Periodontal Disease | Mogibacterium | + |
| Periodontal Disease | Nanoperiomorbus | + |
| Periodontal Disease | Nanosynbacter | - |
| Periodontal Disease | Neisseria | - |
| Periodontal Disease | Parvimonas | + |
| Periodontal Disease | Pauljensenia | - |
| Periodontal Disease | Peptostreptococcus | + |
| Periodontal Disease | Phocaeicola | + |
| Periodontal Disease | Pseudoleptotrichia | - |
| Periodontal Disease | Rothia | - |
| Periodontal Disease | Scardovia | + |
| Periodontal Disease | Selenomonas | + |
| Periodontal Disease | Shuttleworthia | + |
| Periodontal Disease | Streptococcus | - |
| Periodontal Disease | Tannerella | + |
| Periodontal Disease | Treponema | + |
| Periodontal Disease | Veillonella | - |
| Periodontal Disease | W11650 | + |
| Periodontal Disease | W5053 | + |
| Lung Cancer | Abiotrophia | - |
| Lung Cancer | Actinomyces | - |
| Lung Cancer | Aggregatibacter | - |
| Lung Cancer | Alloprevotella | - |
| Lung Cancer | Bifidobacterium | + |
| Lung Cancer | Bulleidia | - |
| Lung Cancer | Campylobacter | - |
| Lung Cancer | Capnocytophaga | - |
| Lung Cancer | Cardiobacterium | - |
| Lung Cancer | Catonella | - |
| Lung Cancer | Centipeda | - |
| Lung Cancer | Corynebacterium | - |
| Lung Cancer | Eikenella | - |
| Lung Cancer | Eubacterium | - |
| Lung Cancer | F0422 | - |
| Lung Cancer | Fusobacterium | - |
| Lung Cancer | Gemella | - |
| Lung Cancer | GU269551 | - |
| Lung Cancer | Haemophilus | - |
| Lung Cancer | HOT-345 | - |
| Lung Cancer | Lachnoanaerobaculum | - |
| Lung Cancer | Lactobacillus | + |
| Lung Cancer | Lautropia | - |
| Lung Cancer | Leptotrichia | - |
| Lung Cancer | Limosilactobacillus | + |
| Lung Cancer | Mogibacterium | - |
| Lung Cancer | Neisseria | - |
| Lung Cancer | Oribacterium | - |
| Lung Cancer | Peptidiphaga | - |
| Lung Cancer | Peptostreptococcus | - |
| Lung Cancer | Porphyromonas | - |
| Head and Neck Cancer | Abiotrophia | + |
| Head and Neck Cancer | Actinomyces | + |
| Head and Neck Cancer | Aggregatibacter | - |
| Head and Neck Cancer | Alloprevotella | - |
| Head and Neck Cancer | Centipeda | - |
| Head and Neck Cancer | Eikenella | + |
| Head and Neck Cancer | Enterococcus | + |
| Head and Neck Cancer | Eubacterium | - |
| Head and Neck Cancer | Filifactor | - |
| Head and Neck Cancer | Fusobacterium | - |
| Head and Neck Cancer | Gemella | + |
| Head and Neck Cancer | Lactobacillus | + |
| Head and Neck Cancer | Lautropia | + |
| Head and Neck Cancer | Leptotrichia | - |
| Head and Neck Cancer | Ligilactobacillus | + |
| Head and Neck Cancer | Limosilactobacillus | + |
| Head and Neck Cancer | Metamycoplasma | - |
| Head and Neck Cancer | Nanosynbacter | - |
| Head and Neck Cancer | Nanosyncoccus | - |
| Head and Neck Cancer | Peptidiphaga | + |
| Head and Neck Cancer | Porphyromonas | - |
| Head and Neck Cancer | Scardovia | + |
| Head and Neck Cancer | Streptococcus | + |
| Head and Neck Cancer | Tannerella | - |
| OSCC | Actinomyces | - |
| OSCC | Anaeroglobus | - |
| OSCC | Bacteroides | + |
| OSCC | Bifidobacterium | + |
| OSCC | Capnocytophaga | + |
| OSCC | Corynebacterium | - |
| OSCC | Dialister | + |
| OSCC | Eubacterium | - |
| OSCC | Haemophilus | - |
| OSCC | Lachnoanaerobaculum | - |
| OSCC | Leptotrichia | - |
| OSCC | Metamycoplasma | + |
| OSCC | Nanosynbacter | - |
| OSCC | Neisseria | - |
| OSCC | Oribacterium | - |
| OSCC | Parvimonas | + |
| OSCC | Pauljensenia | - |
| OSCC | Rothia | - |
| OSCC | Selenomonas | + |
| OSCC | Veillonella | - |
| Breast Cancer | Alloprevotella | - |
| Breast Cancer | Eubacterium | - |
| Breast Cancer | Filifactor | - |
| Breast Cancer | Fusobacterium | - |
| Breast Cancer | Gemella | - |
| Breast Cancer | Scardovia | + |
| Breast Cancer | Treponema | - |
| Breast Cancer | W11650 | - |
| RA | Lancefieldella | + |
| RA | Moryella | + |
| RA | Pauljensenia | + |
| RA | Sphingomonas | - |

### Consensus genus abundance heatmaps

#### Periodontal Disease (46 genera x 2201 samples)

![Periodontal Disease consensus heatmap](figures_linda2/consensus_heatmaps/consensus_heatmap_Periodontal_Disease.png)

#### Lung Cancer (31 genera x 4724 samples)

![Lung Cancer consensus heatmap](figures_linda2/consensus_heatmaps/consensus_heatmap_Lung_Cancer.png)

#### Head and Neck Cancer (24 genera x 1648 samples)

![Head and Neck Cancer consensus heatmap](figures_linda2/consensus_heatmaps/consensus_heatmap_Head_and_Neck_Cancer.png)

#### OSCC (20 genera x 834 samples)

![OSCC consensus heatmap](figures_linda2/consensus_heatmaps/consensus_heatmap_OSCC.png)

#### Breast Cancer (8 genera x 2261 samples)

![Breast Cancer consensus heatmap](figures_linda2/consensus_heatmaps/consensus_heatmap_Breast_Cancer.png)

#### RA (4 genera x 437 samples)

![RA consensus heatmap](figures_linda2/consensus_heatmaps/consensus_heatmap_RA.png)

### Direction conflicts

| Disease | Genus | Direction_MaAsLin3 | Direction_ANCOMBC2 | Direction_LinDA |
| --- | --- | --- | --- | --- |
| Head and Neck Cancer | Campylobacter | - | - | + |
| Head and Neck Cancer | Capnocytophaga | - |  | + |
| Head and Neck Cancer | Corynebacterium | - | + |  |
| Head and Neck Cancer | Granulicatella | - |  | + |
| Head and Neck Cancer | Haemophilus | - | + | - |
| Head and Neck Cancer | Kingella | - | + | + |
| Head and Neck Cancer | Mogibacterium | - |  | + |
| Head and Neck Cancer | Ottowia | - |  | + |
| Head and Neck Cancer | Peptoanaerobacter | - |  | + |
| Head and Neck Cancer | Phocaeicola | - |  | + |
| Head and Neck Cancer | Pyramidobacter | - |  | + |
| Head and Neck Cancer | Rothia | - |  | + |
| Head and Neck Cancer | Saccharimonas | - | + | + |
| Head and Neck Cancer | Slackia | - |  | + |
| Head and Neck Cancer | Veillonella | - |  | + |
| Head and Neck Cancer | W5053 | - |  | + |
| IBD | Pauljensenia | + | - | + |
| OSCC | Sutterella | + | - |  |
| Periodontal Disease | Aggregatibacter |  | - | + |
| Periodontal Disease | Campylobacter | + | - | + |
| Periodontal Disease | UBA1568 | - | + |  |
## 6. Output paths

- `maaslin3/`
- `ancombc2/`
- `linda2/`
- `method_comparison_linda2/`
- `reports/method_comparison_report_linda2.md`
- `reports/method_comparison_report_linda2.html`

## 7. Software versions

| Component | Version | Notes |
| --- | --- | --- |
| maaslin3 | WARNING: ignoring environment value of R_HOME | existing da_v2 MaAsLin3 runner |
 | maaslin3 | 0.99.16 | existing da_v2 MaAsLin3 runner |
| ANCOMBC | 2.4.0 | existing da_v2 ANCOM-BC2 runner |
| MicrobiomeStat | 1.4 | LinDA2 runner (`is.winsor = FALSE`) |
