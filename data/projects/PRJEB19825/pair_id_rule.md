# PRJEB19825 Pair ID 生成规则

## 配对类型
**年龄时间配对**

## 命名模板
`{subject_id}_{age}_{seq}`

## 数据来源
**SraRunTable.csv → host_subject_id + AGE + collection_timestamp**

## 示例
- `F4_33.13` → 受试者F4，年龄33.13岁
- `M3_33.14_1` → 受试者M3，年龄33.14岁，序号1

## 验证规则
- pair_id = subject + "_" + age + "_" + seq (可选)
- 同一pair_id必须来自同一采集日期

## 已修正错误
- `M3_33.14_3`: ERR1868669(Oral)+ERR1868675(Gut) ✅
- ERR1868672标记为 `M3_33.14_unpaired`

---
**配对状态**: ✅ 455正确配对 (修正后)
