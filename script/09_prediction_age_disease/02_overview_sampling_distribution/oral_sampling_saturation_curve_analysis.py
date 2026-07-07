import pandas as pd
import numpy as np
import biom
import random
import os

# === 路径配置 ===
INPUT_DIR = "/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/input/oral_input"
OUTPUT_DIR = "/home/hexincheng/projects/hxc/oral_project/oral_plot_by_yl_20260602/oral_code_results_new/output/oral_output"
BIOM_FILE = os.path.join(INPUT_DIR, "biom_oral_qc.biom")
META_FILE = os.path.join(INPUT_DIR, "metadata_oral_qc.tsv")
OUT_FILE = os.path.join(OUTPUT_DIR, "rarefaction_pcofg.tsv")
TAX_FILE = os.path.join(INPUT_DIR, "taxonomy.tsv")

def get_tax_name(tax_list, prefix):
    """从分类学列表中提取指定层级的名称（例如 'p__' 代表门）"""
    for t in tax_list:
        if t.startswith(prefix):
            val = t[len(prefix):].strip()
            # 如果该层级未注释完整，则跳过不计入"已知Unique"中
            if val != "":
                return t 
    return None

def main():
    print("1. Loading metadata...")
    # 使用原生 Python 逐行解析，跳过引发报错的额外列
    valid_samples = set()
    with open(META_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        # 读取表头并定位 sample_id 所在的列索引
        header = f.readline().strip('\n').split('\t')
        if 'sample_id' not in header:
            raise ValueError("metadata 文件中未找到 sample_id 列")
        run_idx = header.index('sample_id')
        
        for line in f:
            parts = line.strip('\n').split('\t')
            # 只要该行长度足够包含 sample_id 列，就提取出来
            if len(parts) > run_idx:
                run_id = parts[run_idx].strip()
                if run_id:
                    valid_samples.add(run_id)

    print("2. Loading BIOM table (this may take a minute)...")
    table = biom.load_table(BIOM_FILE)

    # 筛选在 metadata 和 biom 表中同时存在的样本
    table_samples = set(table.ids(axis='sample'))
    keep_samples = list(valid_samples.intersection(table_samples))
    print(f"   Found {len(keep_samples)} valid matched samples.")

    print("3. Loading external Taxonomy file...")
    # 【新增】：读取独立的 taxonomy.tsv，存入字典以供快速查询
    tax_dict_map = {}
    with open(TAX_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        _ = f.readline() # 跳过第一行表头 (Feature ID  Taxon   Confidence 等)
        for line in f:
            parts = line.strip('\n').split('\t')
            if len(parts) >= 2:
                feat_id = parts[0].strip()
                taxon_str = parts[1].strip()
                tax_dict_map[feat_id] = taxon_str

    # 随机打乱样本顺序，模拟随机抽样累积的过程
    random.seed(42)
    random.shuffle(keep_samples)

    print("4. Parsing taxonomy for each feature...")
    feature_ids = table.ids(axis='observation')
    feature_to_tax = []

    # 预先解析好每个OTU/ASV对应的门纲目科属
    for fid in feature_ids:
        tax_dict = {'Phylum': None, 'Class': None, 'Order': None, 'Family': None, 'Genus': None}

        # 【修改】：直接从刚刚加载的字典里获取物种字符串
        if fid in tax_dict_map:
            tax_str = tax_dict_map[fid]
            tax = [t.strip() for t in tax_str.split(';')]

            tax_dict['Phylum'] = get_tax_name(tax, 'p__')
            tax_dict['Class'] = get_tax_name(tax, 'c__')
            tax_dict['Order'] = get_tax_name(tax, 'o__')
            tax_dict['Family'] = get_tax_name(tax, 'f__')
            tax_dict['Genus'] = get_tax_name(tax, 'g__')

        feature_to_tax.append(tax_dict)

    print("5. Calculating accumulation curve data...")
    # 转换为稀疏矩阵以大幅提升切片计算速度
    mat = table.matrix_data.tocsc()
    sample_id_to_idx = {sid: idx for idx, sid in enumerate(table.ids(axis='sample'))}
    keep_sample_indices = [sample_id_to_idx[sid] for sid in keep_samples]

    seen_phylum = set()
    seen_class = set()
    seen_order = set()
    seen_family = set()
    seen_genus = set()

    total_samples = len(keep_samples)
    
    # 设定需要记录结果的样本节点。
    steps = set(
        [1, 5, 10, 50, 100, 200, 500, 800] +
        list(range(1000, 10000, 1000)) +
        list(range(10000, total_samples, 5000)) +
        [total_samples]
    )

    results = []

    for i, s_idx in enumerate(keep_sample_indices):
        n_samples = i + 1

        # 提取当前样本中丰度大于0的Feature索引
        col = mat[:, s_idx]
        present_features = col.indices

        # 更新已发现的分类群集合 (Set自动去重)
        for f_idx in present_features:
            tax = feature_to_tax[f_idx]
            if tax['Phylum']: seen_phylum.add(tax['Phylum'])
            if tax['Class']: seen_class.add(tax['Class'])
            if tax['Order']: seen_order.add(tax['Order'])
            if tax['Family']: seen_family.add(tax['Family'])
            if tax['Genus']: seen_genus.add(tax['Genus'])

        # 在设定的节点记录当前的 Unique Taxa 数量
        if n_samples in steps:
            results.append({
                'Samples': n_samples,
                'Phylum': len(seen_phylum),
                'Class': len(seen_class),
                'Order': len(seen_order),
                'Family': len(seen_family),
                'Genus': len(seen_genus)
            })

    print("6. Saving results...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('Samples')
    df_results.to_csv(OUT_FILE, sep='\t', index=False)
    print(f"   Done! Data successfully saved to {OUT_FILE}")

if __name__ == "__main__":
    main()