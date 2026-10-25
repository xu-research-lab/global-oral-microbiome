# sex_differential_continent/pooled — pooled-continent sex DA (da_v2)

与父目录（按大洲分层）的差别：不再按大洲各自内部做 F_vs_M，而是把**所有大洲合并**
成一个 saliva F_vs_M 比较，`continent_group` 作为**协变量**进入模型矫正（无
sex×continent 交互项：Africa 无男性 saliva 样本，加交互会秩亏；`prepare_covars()`
会自动把 F-only 的大洲折叠为 "Other"，只贡献截距水平，不影响 sex 主效应）。

- 比较名（唯一）：`saliva_pooled_F_vs_M`（ref=F, target=M；log2FC>0 => M 富集）
- 三方法共识与父版一致：LinDA + Maaslin2 + ANCOMBC2；**robust = padj<0.05 在 ≥2 方法**
- scope：**只跑 hc**；结果根 `results/da_v2/sex_pooled_hc/`
  - `sex_{linda,maaslin2,ancombc2}_pairwise/saliva_pooled_F_vs_M/`
  - `sex_pooled_da_consensus/`（`consensus_robust.tsv` = n2 主输出）
- **不覆盖**父目录的脚本与结果；所有运行一律通过 sbatch 提交到计算节点
  （不要在登录节点直接跑驱动脚本）。
- pooled 专用常量直接内联在 01/02 里（比较名在 01 由 `paste0` 推导，
  与目录名构造路径一致），无额外 helpers 文件。

## 文件

| 文件 | 作用 |
|---|---|
| `01-run-three-method-da.R` | 三方法驱动（单一比较）；新增 `--smoke` 子抽样 |
| `02-merge-and-consensus.R` | 合并 + n2 共识（五阶段同父版） |
| `run_linda_maaslin2.sbatch` | LinDA+Maaslin2（120G / 2 天） |
| `run_ancombc2.sbatch` | ANCOMBC2 单作业（无 array，120G / 3 天——pooled HC ~10.5k 样本是长板） |
| `run_merge_consensus.sbatch` | 02 合并共识（60G / 6h；需 afterok 依赖） |
| `run_smoke.sbatch` | 端到端冒烟（子抽样，结果进 `results/da_v2_smoke_pooled`） |
| `submit_sex_pooled.sh` | 一键提交全链（afterok 依赖） |

## 推荐流程（全部在计算节点上）

1) **冒烟**（分钟级，先验证管线通了）：
   ```bash
   sbatch differential_abundance/sex_differential_continent/pooled/run_smoke.sbatch
   ```
   看日志确认三方法 + 共识全通过后清理：
   ```bash
   rm -rf results/da_v2_smoke_pooled
   ```
2) **全量**：
   ```bash
   bash differential_abundance/sex_differential_continent/pooled/submit_sex_pooled.sh
   ```
3) 监控：`squeue -u $USER`；日志 `results/logs/da_sex_pooled_*`

## 说明

- pooled 无 03 画图阶段（父版 `03-plot-consensus-bubble.R` 未适配 pooled）
- `--only` 仅为兼容保留；单比较无需使用
- Africa 男性 saliva 为 0（hc: F421/M0）→ 自动折叠为 "Other"
- 幂等：`analysis.done` sentinel 已存在的方法自动跳过