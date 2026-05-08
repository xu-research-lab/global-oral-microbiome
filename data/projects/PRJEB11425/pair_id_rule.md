# PRJEB11425 Pair ID 生成规则

## 配对类型
**时间点配对**

## 命名模板
`{subject_id}_{timepoint}`

## 数据来源
**SraRunTable.csv → host_subject_id + 时间标识**

## 示例
- `1017-10_D10 or D/C` → 入院第10天
- `1017-11_ADMISSION` → 入院时

## 验证规则
- pair_id = subject_id + "_" + ADMISSION/D10

---
**配对状态**: ✅ 105正确配对
