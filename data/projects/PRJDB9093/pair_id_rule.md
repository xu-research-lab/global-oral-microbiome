# PRJDB9093 Pair ID 生成规则

## 配对类型
**时间点配对**

## 命名模板
`{subject_id}` 或 `{subject_id}_{timepoint}`

## 数据来源
**SraRunTable.csv → Sample Name列**

## 示例
- HC组: `HD01` (健康对照，单一采样)
- HIV组: `MB101_0w` (治疗前), `MB101_12w` (治疗12周), `MB101_24w` (治疗24周)

## 验证规则
- pair_id = 1 Gut + 1 Oral
- HC组无时间点后缀，HIV组有时间点

---
**配对状态**: ✅ 65正确配对
