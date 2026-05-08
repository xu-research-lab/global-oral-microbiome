# PRJEB24929 Pair ID 生成规则

## 配对类型
**时间点配对**

## 命名模板
`{subject_base}{timepoint_letter}`

## 数据来源
**SraRunTable.csv → Sample Name列 (提取前缀)**

## 示例
- `04A`, `04B`, `04C` → 受试者04，时间点A/B/C
- 提取逻辑: `04A_gut → 04A`

## 验证规则
- pair_id = subject编号 + 时间点字母 (无分隔符)

---
**配对状态**: ✅ 111正确配对
