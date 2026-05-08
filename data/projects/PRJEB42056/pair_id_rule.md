# PRJEB42056 Pair ID 生成规则

## 配对类型
**实验处理配对**

## 命名模板
`{condition}{subject_num}`

## 数据来源
**储存条件标识 + subject编号**

## 示例
- `A1`, `B1`, `C1` → 储存条件A/B/C，subject编号1
- 同一subject在不同储存条件下配对Oral+Gut

## 验证规则
- pair_id = 条件字母 + 编号 (无分隔符)
- 字母(A-G)表示储存方法

---
**配对状态**: ✅ 82正确配对
