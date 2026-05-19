# PRJEB27304 Pair ID 生成规则

## 配对类型
**时间点配对**

## 命名模板
`{subject_id}_{timepoint}`

## 数据来源
**SraRunTable.csv → Sample Name列**

## 示例
- `111_T1`, `111_T3`, `111_T4` → 受试者111，时间点T1/T3/T4

## 验证规则
- pair_id = subject_id + "_" + T编号
- T编号表示采样时间点

---
**配对状态**: ✅ 212正确配对
