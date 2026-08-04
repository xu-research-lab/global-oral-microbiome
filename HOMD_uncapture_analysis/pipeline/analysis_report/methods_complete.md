# 完整方法流程：从 ASV 到 56 个候选新型口腔 Phylotypes

## Phase 1 — 候选 ASV 生成

### Step 1: ASV 收集 + 精确去重
**脚本**: `01_collect_asvs_and_abundance.py`

从 334 个 BioProjects 读取 DADA2 产出的 `DADA2.biom`，提取所有 ASV 序列及其丰度，按**完全相同的序列**（100% identical）去重。同一条 16S 序列出现多次时只保留一个 `DEREP_XXX` 编号。

| 输入 | 输出 |
|------|------|
| `{bioproject}/results/DADA2.biom` (334 个) | `exact_derep_summary.tsv` |
| `metadata-v2-merge-oral.csv` | `exact_derep_representatives.fa` |
| | `asv_abundance_metadata.tsv` |

**输入**: 3,396,776 raw ASVs → **输出**: 2,844,174 derep sequences

### Step 2: HOMD 比对（identity-only）
**脚本**: `02_homd_mapping.sbatch`

vsearch `--usearch_global --iddef 1` 将全部 2,844,174 条 derep 序列比对到 HOMD V16.03 参考库（6,600 条 16S 序列）。

**HOMD-uncaptured 定义**（per 序列）：
- identity < 97%，**或**
- 无 HOMD hit

> qcov 不使用。global identity = (alnlen - mism - gaps) / query_length × 100

| 输入 | 输出 |
|------|------|
| `exact_derep_representatives.fa` | `exact_derep_homd_annotated.tsv` |
| HOMD V16.03 full (6,600 seqs) | `homd_mapping_stats.json` |

HOMD-captured (id ≥ 97%): 997,603 (35.1%)
**HOMD-uncaptured**: **1,843,301** (64.9%)

随后按 support floor 筛选: ≥2 BioProjects AND ≥50 reads → **69,820** candidates

### Step 3: 嵌合体去除
**脚本**: `03_chimera_removal.sbatch`

vsearch `--uchime_ref` 将 69,820 条 ASV 比对到合并参考库（HOMD 6,600 + GG2 337,506 + SILVA 510,495 ≈ 854,601 序列）。DADA2 只能去除单次 run 内的嵌合体，此步骤捕获**跨研究/跨分类群**的嵌合体。

| 输入 | 输出 |
|------|------|
| `exact_derep_representatives.fa` (69,820) | `exact_derep_chimera_annotated.tsv` |
| HOMD + GG2 + SILVA combined | `nonchimeras.fa` |

嵌合体: 15,566 (22.3%) → **非嵌合体**: **54,254**

---

## Phase 2 — 污染防火墙（两级标记，不删除）

### Step 4: 97% 聚类对照 ASV
**工具**: vsearch `--cluster_fast --id 0.97`

收集 20 个有对照样本的 BioProjects 中全部的对照 ASV（1,469 个对照样本），与 54,254 个 candidate ASV 合并做 97% OTU 聚类。如果 candidate ASV 与至少一个 HOMD-uncaptured 的对照 ASV 落在同一个 97% 簇内，则该 candidate 被打上 `cluster97_contaminant` 标记。

> HOMD-captured 对照 ASV（对照中出现的口腔菌）不被视为污染证据——它们更可能是交叉污染，而不是外源污染。

- 标记: 20,994 / 54,254 (38.7%)

### Step 5: 扩展属级参考库
**输入**: HOMD taxonomy（251 个已确认口腔属）+ 文献黑名单（94 属）+ 对照样本中的 ASV 属注释

找出在对照样本中出现、但**不在 HOMD 中**的属（891 个），加入扩展参考库。Candidate ASV 的 SILVA 属注释若命中了这个扩展参考库，则打上 `genus_contaminant` 标记。

- 标记: 2,754 / 54,254 (5.1%)
- 保护: 86 个口腔核心属 (never removed)

### Step 6-8: decontam + 文献黑名单 + 跨实验室污染目录 + oral-vs-blank
**脚本**: `flag_contaminants.py` → `update_integrated_annotations.py`

| 证据层 | 方法 |
|--------|------|
| **decontam (属级)** | decontam prevalence mode (threshold=0.5) 在 20 个对照研究的属级进行池化 |
| **跨实验室污染目录** | 在多个独立实验室的对照中出现的 taxa → 标记 `catalog_contaminant` |
| **文献黑名单** | `contaminant_reference_genera.tsv`（Salter 2014 + Weyrich 2019）→ auto_remove / flag_review / oral_core |
| **oral-vs-blank residence** | oral_prevalence ≥ 0.01 AND oral_ab > blank_ab AND blank_prevalence ≤ 0.10 |

**所有 54,254 条 ASV 都被保留**，每条 ASV 增加 16 列污染证据标注。

**最终标记分布**:
| 标记 | 数量 | 
|------|------|
| clean | 32,812 |
| soft_flag_cluster97 | 18,925 |
| hard_flag_both | 2,011 |
| soft_flag_genus | 420 |
| protected_oral_core | 86 |

**输出**: `candidates_flagged_updated.tsv` (54,254 rows × 16 cols)

---

## Phase 4 — 系统发育 Phylotype 分析

### Step 12: SEPP 系统发育插入
**脚本**: `05c_sepp_place.sbatch`

QIIME2 `q2-fragment-insertion sepp` 将 54,254 条 ASV 插入到 Greengenes2 (2022.10) 参考树骨架（~385K GG2 参考序列）中。

**SEPP 原理**:
1. 对每条 ASV，从 GG2 树的邻居分支构建 **HMM profile**（捕捉位点保守模式）
2. 通过 **pplacer**（最大似然）在树上搜索最佳插入位置
3. 输出 `edge_num`（插入的分支编号）和 `pendant_length`（插入长度）

54,254 submitted → 54,132 placed (99.8%)

**输出**: `sepp_insertion_tree.qza` → exported as `tree_export/tree.nwk` (20 MB, ~439K tips = 385K GG2 + 54K ASV)

### Step 13: Phylotype 定义 + Gates
**脚本**: `05d_parse_sepp_placements.py`

从 `placements.json` 提取每个 ASV 的 best placement（最高 `like_weight_ratio`）→ `edge_num` → phylotype ID = `EDGE_{edge_num}`。

**同一分支上的 ASV = 一个 phylotype**——这利用了 SEPP 的最大优势：不同扩增区域（V4 vs V3-V4）来的同源 ASV 会自动落在同一个分支上，因为它们在系统发育上是同源的。

54,132 ASVs → **23,088 phylotypes** (edge clusters)

随后应用 **5 道 gates**:
| Gate | 条件 | 通过 |
|------|------|------|
| 1 | HOMD-gap（最近 GG2 参考不在 HOMD） | 22,418 |
| 2 | ≥ 10 subjects | 10,980 |
| 3 | ≥ 2 BioProjects | 10,980 |
| 4 | ≥ 1 member passing oral-vs-blank residence | **56** |
| 5 | ≥ 1 non-chimeric member | 56 |

### Step 14: Novelty 多数据库鉴定
**脚本**: `07_novelty_sepp/`

vsearch 全局比对 56 个 phylotype 的 representative 序列 vs **HOMD + Greengenes2 + SILVA** 三个独立数据库。

| 分类 | 数量 | 定义 |
|------|------|------|
| **novel_all_refs** | **14** | < 97% 到所有三个数据库 → 候选新物种/谱系 |
| **known_not_in_homd** | **42** | ≥ 97% 到 GG2/SILVA 但不匹配 HOMD → 口腔目录新增 |

---

## 完整数据流

```
3,396,776 ASV raw (334 BPs)
    ↓ exact derep (100% identity)
2,844,174 derep sequences
    ↓ vsearch vs HOMD (id < 97% OR no hit)
1,843,301 HOMD-uncaptured
    ↓ support floor (≥2 BP, ≥50 reads)
69,820 candidates
    ↓ uchime_ref (vs HOMD+GG2+SILVA)
54,254 non-chimeric
    ↓ Contamination Firewall (FLAGS ONLY, no deletion)
    │  → 97% cluster with HOMD-uncaptured control ASVs
    │  → Genus in extended reference (not in HOMD)
    │  → Literature blacklist
    │  → Oral-vs-blank residence test
54,254 (flags added)
    ↓ SEPP placement on GG2 tree
54,132 placed (99.8%)
    ↓ Edge clustering
23,088 phylotypes
    ↓ 5 Gates (HOMD-gap, subjects≥10, BPs≥2, residence≥1, non-chimera)
**56 candidate novel oral phylotypes**
    ↓ Multi-DB novelty classification
14 novel_all_refs + 42 known_not_in_homd
```

## 关键脚本与输出文件

| Step | 脚本 | 核心输出 |
|------|------|----------|
| 1 | `01_collect_asvs_and_abundance.py` | `exact_derep_summary.tsv` (2,844,174 行) |
| 2 | `02_homd_mapping.sbatch` | `exact_derep_homd_annotated.tsv` |
| 3 | `03_chimera_removal.sbatch` | `exact_derep_chimera_annotated.tsv` (54,254 candidates) |
| 4-8 | `flag_contaminants.py` → `update_integrated_annotations.py` | `candidates_flagged_updated.tsv` (54,254 rows × 16 cols) |
| 12 | `05c_sepp_place.sbatch` | `tree_export/tree.nwk` (385K GG2 + 54K ASV tips) |
| 13 | `05d_parse_sepp_placements.py` | `sepp_placements_parsed.tsv` |
| 14 | `07_novelty_sepp/` | `phylotypes_novelty_classified.tsv` (56 phylotypes) |

## 最终 56 个候选新型口腔 Phylotypes

| 分类 | 数量 | 含义 |
|------|------|------|
| **novel_all_refs** | 14 | 所有数据库都找不到 → 候选新物种 |
| **known_not_in_homd** | 42 | GG2/SILVA 中有但 HOMD 无 → 口腔目录缺失 |

详见报告: `pipeline_report_comprehensive.md`

---

*完整方法文档 — 2026-07-29*
