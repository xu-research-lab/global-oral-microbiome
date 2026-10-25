# Random Forest OVO/OVA (Oral + Gut)

Project-level CV pipeline for region prediction.

## Project-level CV (v2-style, cross-study)

Same design as  
`/beegfs/huangxiaochang/global_oral/new_analysis/randomforest_v2`:

- Split unit: **BioProject / project** (whole study train or test)
- Folds: up to 5, `n_folds = min(5, n_pos_projects, n_neg_projects)`
- Repeats: 5 project shuffles
- Filter: site×region×project with ≥20 samples
- Groups: oral `bioproject_id`, gut `project`

```bash
python continent_rf/scripts/run_oral_gut_project_cv_repeat5.py
python continent_rf/scripts/plot_oral_gut_project_cv.py
```

Results: `continent_rf/result/oral/`

## Shared model / task design

| Item | Setting |
|------|---------|
| Target | continent / region (6 regions) |
| Oral sites | Saliva, Tongue, Tooth (separate models) |
| Gut | single site |
| OVA | each region vs Rest |
| OVO | each region pair |
| Model | RF 200 trees, `class_weight=balanced`, 8 jobs |
| Features | sample-normalized relative abundance |
