# 备选方案：OTU 聚类定 Support Floor（待讨论）

**提出日期**: 2026-07-31

## 问题

不同 BP 单独跑 DADA2，不同测序区域（V3-V4 / V4 / V6-V8 等）的同一菌 ASV 序列会有差异，单 ASV-level support floor 会误杀跨区域/跨 BP 变体。

## 数据证据

949K clean ASV 做 97% self-clustering：
- 336,484 clusters（2.8 ASVs/cluster）
- 31.6% clusters 含 ≥2 ASV → 确实有大量跨区域变体被合并

Cluster-level ≥2 BP, ≥50 reads → 62,658 clusters → 只需 SEPP 62K 代表序列（vs 34K 单 ASV 方案）

## 方法

1. 97% OTU self-clustering 所有 clean ASV
2. Cluster-level support floor（总 BP 数、总 reads、总 subjects）
3. 每个 cluster 选代表序列跑 SEPP
4. 成员继承 placement

## 潜在风险

- 97% identity ≠ 系统发育同源性，非同种菌可能落在同一 cluster
- 成员"继承"代表的 placement 不完全等同于独立 SEPP 结果

## 计算资源对比

| 方案 | SEPP 输入 | 预计耗时 |
|------|----------|----------|
| 全部 ASV SEPP | 608K | ~11 天 |
| 代表序列 SEPP | 62K | ~1 天 |
| 单 ASV support floor | 34K | ~0.5 天 |

**结论**: 等计算资源充足后，建议跑 62K 代表序列方案验证。
