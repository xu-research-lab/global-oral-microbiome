# Metadata Cleaning Mapping Summary

This document summarizes the age grouping and status cleaning mappings applied in `metadata_v2_csv_merge_clean.R`.

---

## 1. Age Grouping Mapping (`age_cleaned_group`)

### 1.1 Direct age value parsing (`age_cleaned_group_base`)

When `age` field contains a numeric value:

| Condition | `age_cleaned` | Group |
|-----------|---------------|-------|
| Contains "week" | Divide by 52 | — |
| `age_numeric < 0` or `> 150` | — | Anomaly (→ NA) |
| `age_numeric ≤ 0.5` | — | 0-6 months |
| `0.5 < age_numeric ≤ 3` | — | 6-36 months |
| `3 < age_numeric ≤ 12` | — | 3-12 years |
| `12 < age_numeric ≤ 18` | — | 12-18 years |
| `18 < age_numeric ≤ 60` | — | 18-60 years |
| `age_numeric > 60` | — | ＞60 years |

### 1.2 Age info text parsing (`age_info_cleaned`)

When direct age is unavailable, parse `age_info` text (case-insensitive):

| Text Pattern | Parsing Rule | Group |
|--------------|--------------|-------|
| `not provided`, `unspecified`, `age_group_1`, `toddlers 1-3 years; children 4-12 years; teenager-adults 13-59 years`, `full cohort median age 41.5 years...` | Explicit NA | NA |
| `15-40 years`, `11-33 years`, `18-50 years` (digit-digit years) | First number threshold | Range-based |
| `>=18`, `≥21`, `56+` (thresholds) | Number threshold | ≥60→＞60; >18→18-60; >12→12-18; >3→3-12; >0.5→6-36m |
| `above 18 years old`, `over 18 years` | Number threshold | Same as above |
| `median age: 61`, `average 48 years`, `mean 35 years` | Extract number | Range-based |
| `0-1.5`, `18-45` (pure numeric range) | Mean of range | Mean-based |
| `(35, 70)` (parenthesis range) | First number | Range-based |
| `6 months`, `12 months` (contains "month") | Month count | ≤6→0-6m; ≤36→6-36m |
| `average 48`, `mean 44` (no "years", no ±) | Extract number | Range-based |
| `20s`, `30s`, `40s` (decades) | Extract number | <12→3-12y; <18→12-18y; <60→18-60y |
| `55.39±7.80` (single ±) | First number | Range-based |
| `mean 35 +- 10` | First number | ≤18→12-18y; ≤60→18-60y |
| `between 32 and 56 years` | Mean of range | Mean-based |
| `60+`, `50+`, `above 60` (X+) | Extract number | ≥60→＞60y; >18→18-60y |
| `adult`, `adults` | Text match | 18-60 years |
| `infant`, `baby`, `newborn`, `preterm infant` | Text match | 0-6 months |
| `children`, `child`, `chilren` | Text match | 3-12 years |
| `elderly adults`, `elderly people and adults`, `elder` | Text match | ＞60 years |
| `teen`, `teenager` | Text match | 12-18 years |
| `young adults` | Text match | 18-60 years |
| `centenarians` | Text match | ＞60 years |
| `≥2 occurrences of ±` | Text match | 18-60 years |
| `child (12-year-olds)` (parenthesis format) | Extract number | Range-based |
| Any remaining text with numbers (catch-all) | First number | Range-based |

**Range-based group thresholds:** ≤0.5→0-6m; ≤3→6-36m; ≤12→3-12y; ≤18→12-18y; ≤60→18-60y; >60→＞60y

### 1.3 Final priority (`age_cleaned_group`)

1. `age_numeric` valid → use `age_cleaned_group_base`
2. `age_info_cleaned` available → use `age_info_cleaned`
3. `raw_host == "pregnant women"` → 18-60 years
4. `age_numeric` is Anomaly → Anomaly (final standardization converts to NA)
5. Otherwise → NA

---

## 2. Status Cleaning Mapping (`status_cleaned`)

### 2.1 Mapped to NA (non-disease / control / descriptive status)

| Original Value | Cleaned |
|----------------|---------|
| `mucositis healed`, `pre-conditioning`, `pre-cancer`, `atopic non-asthmatic controls` | NA |
| `allogeneic hematopoietic stem cell transplantation recipient`, `immunological ART responder` | NA |
| `dentally compromised`, `pulmonary nodule`, `thyroid nodules`, `black dental stain` | NA |
| `epithelial precursor lesion`, `oral benign tumors` | NA |
| `methamphetamine addiction`, `down syndrome`, `autism spectrum disorder` | NA |
| `chronic fatigue syndrome`, `chronic recurrent multifocal osteomyelitis`, `severe congenital neutropenia` | NA |
| `myasthenia gravis`, `primary biliary cholangitis`, `functional gastrointestinal complaints` | NA |
| `missing` | NA |
| `tonsillar hyperplasia`, `tonsillar hypertrophy` | NA |
| `migraine` | NA |
| `polycystic ovary syndrome` | NA |
| `pre-eclampsia` | NA |
| `breastfed...`, `never breastfed`, `orthodontic...`, `non-responder`, `high-risk for ra` | NA |
| `never smoker`, `current smoker` | NA |
| `dental implant` (PRJNA417321) | NA |
| `oral benign tumors` (PRJNA841771) | NA |
| `healthy` | HC |

### 2.2 Birth status mappings

| Original Value | Cleaned | `disease_study_cleaned` |
|----------------|---------|------------------------|
| `full-term` | **HC** | — |
| `preterm`, `preterm newborns` | **Preterm** | **yes** |

### 2.3 Mapped to "Other" (rare / non-target diseases)

| Original Value | Cleaned |
|----------------|---------|
| `halitosis`, `dysphagia`, `watery diarrhea NOS`, `truncal ataxia NOS` | Other |
| `cancer patient`, `Carcinoma in situ`, `Chronic mucosal inflammation`, `acute pancreatitis` | Other |
| `aortic thrombosis`, `tibial fracture`, `hip fracture`, `hypotension` | Other |
| `STEMI`, `NSTEMI`, `temporal lobe epilepsy` | Other |
| `Fanconi anemia with leukoplakia`, `Fanconi anemia with OSCC` | Other |
| `systemic lupus erythematosus with chronic periodontitis`, `systemic lupus erythematosus` | Other |
| `acute lymphoblastic leukemia`, `helminth infection`, `stomatitis` | Other |
| `early childhood caries and congenital heart disease`, `proliferative verrucous leukoplakia` | Other |
| `Refractory Celiac disease`, `IBD+Diabetes` | Other |
| `bradycardia NOS`, `juvenile idiopathic arthritis`, `chronic tonsillitis` | Other |
| `liver cirrhosis`, `chronic kidney disease`, `immunoglobulin A nephropathy` | Other |
| `membranous nephropathy`, `liver carcinoma`, `gliomas` | Other |

### 2.4 Exact disease mappings

| Original Value | Cleaned |
|----------------|---------|
| `Lung Cancer`, `lung cancer` | Lung Cancer |
| `Breast Cancer`, `breast cancer` | Breast Cancer |
| `Pancreatic Cancer`, `pancreatic cancer` | Pancreatic Cancer |
| `Oral Cancer`, `oral cancer` | Oral Cancer |
| `Digestive Tract Cancer`, `digestive tract cancer` | Digestive Tract Cancer |
| `Prostate Cancer`, `prostate cancer` | Prostate Cancer |
| `Esophageal Cancer`, `esophageal cancer` | Esophageal Cancer |
| `Colon Cancer`, `colon cancer` | Colon Cancer |
| `viremic ART-untreated`, `HIV`, `HIV-1`, `immunological ART non-responder` | HIV |
| `GERD`, `Reflux Esophagitis`, `reflux esophagitis`, `functional dyspepsia and gastroesophageal reflux disease` | Gastroesophageal Reflux Disease |
| `RAU` | Aphthous |
| `chronic periodontitis associated with type 2 diabetes mellitus` | Diabetes and Periodontitis |
| `CD and UC` | CD and UC |
| `caries`, `caries active`, `pre-caries`, `severe dental caries`, `cerebral palsy with mild-moderate/severe caries` | Dental Caries |
| `recurrent aphthous ulcer` | Aphthous |
| `RA` | RA |
| `CD` | CD |
| `UC` | UC |
| `T2DM`, `Type 1 Diabetes`, `type 2 diabetes mellitus` | Diabetes |
| `ESCC` | Esophageal Squamous Cell Carcinoma |
| `DTC` | Differentiated Thyroid Carcinoma |
| `axSpA` | Axial Spondyloarthritis |
| `C.diff` | Clostridium difficile infection |
| `IBD+C.diff` | IBD and Clostridium difficile infection |
| `Sjögren's syndrome` | Sjögren's syndrome |
| `obstructive sleep apnea syndrome` | Obstructive sleep apnea |
| `hand, foot and mouth disease symptomatic/asymptomatic` | Hand, foot and mouth disease |
| `healthy` | HC |
| `no-OM` | HC |
| `cerebral palsy caries-free` | HC |
| `certified saturation divers with valid health certificates` (PRJEB40804) | HC |
| `first-degree relatives` (PRJNA656387) | HC |
| `never smoker` / `current smoker` (PRJNA316469 & disease_study==no) | HC |

### 2.5 Regex-based keyword mappings

| Keyword Pattern | Cleaned |
|-----------------|---------|
| `periodontitis` / `gingivitis` / `periodontal` / `peri-implantitis` | Periodontal Disease |
| `caries` | Dental Caries |
| `coronavirus disease of 2019` / `covid` | COVID-19 |
| `chagas` | Chagas disease |
| `obese` / `obesity` | Obesity |
| `aphthous` | Aphthous |
| `type 2 diabetes mellitus` / `t2dm` / `type-2 diabetes` / `type 1 diabetes` / `diabetes` | Diabetes |
| `atherosclerosis` | Cardiovascular Disease |
| `oral mucositis` / `severe-OM` / `grade 1 OM` | Oral Mucositis |
| `olp` / `erosive OLP` / `reticular OLP` / `Non-erosive OLP` | Oral Lichen Planus |
| `helicobacter pylori` | Helicobacter pylori infection |
| `thyroid cancer` | Thyroid Cancer |
| `cancer` (fallback) | Cancer |
| `reflux esophagitis` | Gastroesophageal Reflux Disease |
| `axspa` | Axial Spondyloarthritis |
| `nasopharyngeal carcinoma` / `Nasopharyngeal carcinoma (NPC)` | Nasopharyngeal Carcinoma |
| `oscc` / `oral squamous cell carcinoma` | OSCC |
| `head and neck squamous cell carcinoma` / `hnscc` / `head and neck cancer` / `hnc` | Head and Neck Cancer |
| `cleft lip and palate` / `CLP` | Cleft Lip and Palate |
| `ad` / `alzheimer` | Alzheimer's disease |
| `parkinson` | Parkinson's disease |
| `post-traumatic stress disorder` / `ptsd` | PTSD |
| `celiac disease` | Celiac disease |
| `obstructive sleep apnea` | Obstructive sleep apnea |
| `copd` | COPD |
| `alcohol use disorder` | Alcohol Use Disorder |
| `atopic dermatitis` | Atopic Dermatitis |
| `atopic asthma` | Asthma |

### 2.6 Case correction after `tools::toTitleCase()`

| After TitleCase | Final |
|-----------------|-------|
| `Hiv` / `Hiv-1` | HIV |
| `Oscc` | OSCC |
| `Hc` | HC |
| `Cd And Uc` | CD and UC |
| `Covid-19` | COVID-19 |
| `Ra` | RA |
| `Cd` | CD |
| `Uc` | UC |
| `Ad` | AD |
| `Copd` | COPD |
| `Escc` | ESCC |
| `Dtc` | DTC |
| `Hnc` | HNC |
| `Hnscc` | HNSCC |
| `Olp` | OLP |
| `Rau` | RAU |
| `T2dm` | T2DM |
| `Axspa` | AxSpA |
| `Crc` | CRC |
| `Gerd` | GERD |
| `Ibd` | IBD |
| `Ptst` | PTSD |

### 2.7 Fallback HC assignment

When `status` is NA, the following bioprojects are forced to **HC**:

`PRJEB1799`, `PRJEB20889`, `PRJEB20890`, `PRJEB20894`, `PRJEB27304`, `PRJEB28088`, `PRJEB28769`, `PRJEB33119`, `PRJEB42056`, `PRJEB6518`, `PRJNA255332`, `PRJNA395436`, `PRJNA427101`, `PRJNA634162`, `PRJNA689848`, `PRJNA721325`, `PRJNA747639`, `PRJNA778070`, `PRJNA828039`, `PRJNA841067`

### 2.8 Preterm disease study override

When `status` is `preterm` or `preterm newborns`, `disease_study_cleaned` is forced to **yes** regardless of original value.
