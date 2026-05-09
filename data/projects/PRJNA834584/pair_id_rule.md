# PRJNA834584 Pair ID 生成规则

## 配对类型
**Subject ID配对**

## 命名模板
`{subject_id}_{marker}`

## 数据来源
**SraRunTable.csv → Sample Name列 + 16S标记**

## 示例
- `HM0245_16S` → 受试者HM0245，16S测序

## 验证规则
- pair_id = subject_id + "_16S"
- 大规模口腔-肠道传播研究

---
**配对状态**: ✅ 470正确配对
