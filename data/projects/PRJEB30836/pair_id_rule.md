# PRJEB30836 Pair ID 生成规则

## 配对类型
**Subject ID配对**

## 命名模板
`{prefix}{num}`

## 数据来源
**SraRunTable.csv → Sample Name列 (缩写转换)**

## 示例
- `Beb1`, `Beb2` → BEB01转换为Beb1

## 验证规则
- pair_id = sample_name缩写形式
- 无时间点后缀 (单一采样)

---
**配对状态**: ✅ 106正确配对
