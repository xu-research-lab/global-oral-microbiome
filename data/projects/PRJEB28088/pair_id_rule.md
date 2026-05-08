# PRJEB28088 Pair ID 生成规则

## 配对类型
**时间点配对**

## 命名模板
`{subject_id}_{days}_d`

## 数据来源
**SraRunTable.csv → host_subject_id + 天数信息**

## 示例
- `A7_1_d` → 受试者A7，第1天
- `A5_13_d` → 受试者A5，第13天

## 验证规则
- pair_id = subject + "_" + days + "_d"
- `d` = 天数标记

---
**配对状态**: ✅ 99正确配对
