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

### 问题pair_id列表
| pair_id | 样本数 | 问题类型 |
|---------|--------|----------|
| 12-Jun_22 | 3 | 1 Gut + 2 Oral |
| 3-Jun_22 | 3 | 1 Gut + 2 Oral |
| 4-Jun_21 | 3 | 1 Gut + 2 Oral |
| 5-Jun_30 | 3 | 1 Gut + 2 Oral |
| 6-Jun_104 | 3 | 1 Gut + 2 Oral |
| 7-Jun_85 | 3 | 1 Gut + 2 Oral |
| 8-Jun_106 | 3 | 1 Gut + 2 Oral |

---
**配对状态**: ⚠️ 53正确配对, 7待修正 (共60个pair_id)
