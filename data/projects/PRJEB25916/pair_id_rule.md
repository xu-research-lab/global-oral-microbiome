# PRJEB25916 Pair ID 生成规则

## 配对类型
**年龄时间配对**

## 命名模板
`{subject_id}_{age}`

## 数据来源
**SraRunTable.csv → sample_name + Age**

## 示例
- `1-Jun_23` → 受试者1-Jun，年龄23岁 (young)
- `10-Jun_101` → 受试者10-Jun，年龄101岁 (centenarian)

## 验证规则
- pair_id = sample_name + "_" + Age
- site_code: O=Oral, G=Gut, F=Forehead, OW=Oral Wash

## ⚠️ 待修正错误
- 7个pair_id有1Gut+2Oral (重复测序Run)
- 原因: 同一生物样本被测序两次(BioSample ID相同)
- 建议: 删除重复Run或合并数据

---
**配对状态**: ⚠️ 53正确配对, 7待修正
