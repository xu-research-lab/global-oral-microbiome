# 婴幼儿年龄差异与数据集来源核查

核查日期：2026-10-06。
对象：`figures_saliva_target_top10_globalpeak_maaslin2coef/figureB_abundance_saliva_2method_combined_zscore.png` 中归入 6–36 月龄的 10 个菌属。

当前证据支持这些菌属随婴幼儿发育增加的趋势，单纯的项目间混合不足以解释全部差异。但原始年龄精度和纵向样本的相关性限制了对“0–6 vs 6–36 月龄”及其显著性的解释。

## 图中方向的含义

- 气泡比较的是每个年龄组与其余所有有效年龄组，而不是 0–6 与 6–36 月龄直接比较。
- 这 10 个菌属本来就是按 6–36 月龄显著正效应、至少另一方法同向显著、且该组平均丰度为全组唯一峰值筛选的。因此这个区块不能代表全部菌属的变化。
- 原模型中，这 10 个菌属在 0–6 月龄比较中均为负，在 6–36 月龄比较中均为正，MaAsLin2 FDR 均 <0.05；见 `original_da_results_selected_genera.tsv`。
- 热图的 z-score 是对各年龄组原始平均相对丰度按菌属标准化；它不是批次校正后的效应值。相对丰度也不能直接解释为绝对菌量。

## 有多少个数据集

全图分析集：17,360 份样本，137 个 BioProject，均为当前元数据定义的 HC 唾液样本且有有效年龄组。

| 年龄组 | 样本数 | BioProject 数 |
|---|---:|---:|
| 0–6 months | 786 | 7 |
| 6–36 months | 543 | 8 |
| 3–12 years | 1,928 | 23 |
| 12–18 years | 1,316 | 13 |
| 18–60 years | 8,587 | 113 |
| >60 years | 4,200 | 36 |

前两个年龄组合计来自 12 个不同 BioProject；项目可以跨年龄组，不能直接把 7 和 8 相加。只有 3 个项目同时包含这两个年龄组。

| 共同项目 | 0–6 月龄样本数 | 6–36 月龄样本数 | 平均丰度增加的菌属数 |
|---|---:|---:|---:|
| PRJEB35824 | 360 | 272 | 10/10 |
| PRJNA997934 | 313 | 194 | 10/10 |
| PRJNA770071 | 2 | 23 | 10/10 |

前两个项目分别占低月龄组 85.6%、较大月龄组 85.8%。第三个项目低月龄样本只有 2 份，支持力度有限。BioProject 数也不自动等同于独立队列数或受试者数。

## 项目内部与排除项目后的结果

两个主要项目分别分析时，10 个菌属在较大月龄组的平均相对丰度全部更高，检出率也全部更高。下面给出部分平均相对丰度（百分比）：

| 菌属 | PRJEB35824：低月龄 → 较大月龄 | PRJNA997934：低月龄 → 较大月龄 |
|---|---:|---:|
| Neisseria_563205 | 1.548 → 10.509 | 9.226 → 11.678 |
| Actinomyces | 0.139 → 2.395 | 0.623 → 1.820 |
| Capnocytophaga | 0.017 → 1.623 | 0.118 → 1.251 |
| Lautropia | 0.00294 → 0.982 | 0.0649 → 0.915 |

分别剔除任一主要项目、同时剔除两个主要项目，以及仅剔除 PRJNA997934 的 AGE=0 样本后，10 个菌属的平均丰度方向均保持一致。同时剔除两个项目后剩余低月龄样本 113 份、较大月龄样本 77 份。

这是描述性敏感性检查，没有重新拟合原来的三种差异丰度模型；不能将方向一致说成“排除后仍全部显著”。剔除两个主要项目后的样本主要来自不同项目，残余项目混杂仍然存在。这些菌属也是由完整数据预先筛选出的，当前检查不构成独立验证。

## 已确认的解释限制

1. **PRJNA997934 年龄精度不支持半年分箱。** 原始 `SraRunTable.csv` 中，纳入低月龄组的 313 份样本均为 `AGE=0, age_cat=0Y`，占当前低月龄组 39.8%。当前清洗规则直接将 `age_numeric <=0.5` 归入 0–6 月龄。整岁 0Y 不能建立这些样本属于 0–6 月龄的证据，也不能据此确定哪些样本实际超过 6 个月。原论文按月龄分析婴儿，并以 11 月龄作为部分分析参照；当前本地 SRA 表没有月龄字段。整岁 3 岁记录也不能精确判断是否超过 36 个月。
2. **PRJEB35824 是纵向研究，632 份婴幼儿样本不是 632 个独立受试者。** 原研究采样时间为出生后 2 天、3 月、18 月、3 岁、5 岁。当前元数据中低月龄组为 2 天和 3 月，较大月龄组为 18 月和 3 岁，因此不能从本图推出“6 个月时突然转折”。项目整理文件说明无法获得跨访视的个体编号；分析元数据把该项目的 subject_id 填为 sample_id，不能据此恢复个体配对。
3. **模型已考虑项目，但未充分解决以上问题。** 两个年龄比较的实际 `run_config.tsv` 均包含性别、平台、扩增区域、BMI 组、大洲，以及 `(1|bioproject_id)`。项目随机截距不会纠正年龄误分，也不等于对同一个人多次采样的建模。三种方法使用同一批样本，多方法一致不等于多个独立队列验证。

项目内部方向一致，与婴幼儿口腔菌群成熟的解释相容。PRJEB35824 原论文也描述了 Actinomyces、Capnocytophaga、Corynebacterium、Leptotrichia、Neisseria 等从 18 月龄开始出现/增加；这是同一项目的既有结果，不能算额外独立验证。

## 论文使用建议

目前可描述为：这些菌属在较大月龄婴幼儿中呈较高相对丰度和检出率；主要项目内部也观察到同向变化。

正式以“0–6 vs 6–36 月龄”为结论前，应补充 PRJNA997934 的精确月龄；无法补齐时，将不能明确分组的样本排除出精细年龄分析或另作粗年龄分组敏感性分析。对可识别受试者的纵向项目使用个体相关性模型；无法识别的项目应做单独排除分析。之后重新拟合年龄比较，并检验 0–6 与 6–36 月龄的直接对比。此次核查未改变原始元数据、差异分析结果或原图。

## 文件与复现

- 复现：`/home/xuhuan/miniconda3/envs/R4/bin/Rscript script/age_differential_ova/04-audit-infant-project-effects.R`
- `sample_counts_by_age.tsv`、`sample_counts_by_project_age.tsv`、`infant_project_counts.tsv`：样本与项目分布。
- `within_project_age_comparison.tsv`：各共同项目中完整的均值、检出率、样本数。
- `descriptive_exclusion_sensitivity.tsv`：项目排除后的描述性结果。
- `infant_recorded_age_counts.tsv`、`PRJNA997934_deposited_age_counts.tsv`：实际年龄记录。
- `project_mean_abundance_audit.pdf`、`project_mean_abundance_audit.png`：项目内部与排除两大项目后的对照图；PNG 为 300 dpi。
- `figure_infant_age_project_boxplots.pdf`、`figure_infant_age_project_boxplots.png`：合并所有婴幼儿样本的单面板箱线图。每个菌属上方白色箱体为当前记录的 0–6 月龄，下方灰色箱体为 6–36 月龄；散点颜色区分 12 个 BioProject，图例在右侧。保留全部零值；箱体为原始相对丰度的中位数和四分位距，须为 1.5 倍 IQR 范围内最远的观测值，坐标采用可包含零值的 pseudo-log 转换。PNG 为 300 dpi，保存尺寸 200 × 165 mm。
- 箱线图复现：`/home/xuhuan/miniconda3/envs/R4/bin/Rscript script/age_differential_ova/05-plot-infant-project-boxplots.R`；对应的样本数据和箱线统计分别保存在 `infant_project_boxplot_samples.tsv`、`infant_project_boxplot_statistics.tsv`。图中共 1,329 份不同样本，每份样本在 10 个菌属中各出现一次。
- 重新计算的平均丰度已与原图保存的热图数据逐组核对一致。

原研究：

- PRJEB35824：[Scientific Reports (2020), DOI: 10.1038/s41598-020-66658-7](https://doi.org/10.1038/s41598-020-66658-7)；本地 `data/projects/PRJEB35824/parper.pdf`。
- PRJNA997934：[Cell (2024), DOI: 10.1016/j.cell.2024.07.008](https://doi.org/10.1016/j.cell.2024.07.008)；本地 `data/projects/PRJNA997934/paper.pdf`，婴儿月龄分析见 Figure S3C。
