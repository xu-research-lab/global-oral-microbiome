import os
import h5py
import pandas as pd
import numpy as np
import datetime

# 文件路径配置
# 输入文件
file_path = '/beegfs/dongbiao/oral/beta_oral_qc_dedup/unifrac/unifrac_unweight_pcoa.hdf5'
# 输出文件夹路径
out_dir = '/home/hexincheng/projects/hxc/oral_project/beta_diversity/beta_diversity_oral_qc_dedup/beta_pcoa_out_20260909/unifrac_unweight/'
os.makedirs(out_dir, exist_ok=True)

# 输出文件
out_raw_csv = f'{out_dir}unifrac-unweight-pcoa-raw-all-20260909.csv'
out_r_csv = f'{out_dir}unifrac-unweight-pcoa-R-input-20260909.csv'
log_file = f'{out_dir}data-conversion-check-20260909.log'

# 同时打印到屏幕并写入日志文件
def log_print(msg, file_handle):
    print(msg)
    file_handle.write(msg + '\n')

with open(log_file, 'w', encoding='utf-8') as f_log:
    log_print(f"=== HDF5 转 CSV 数据提取与校验日志 ===", f_log)
    log_print(f"时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f_log)
    log_print(f"输入文件: {file_path}", f_log)
    log_print(f"输出目录: {out_dir}", f_log)
    
    # ==========================================
    # 第一部分：提取与保存
    # ==========================================
    log_print("\n[1/3] 开始提取 HDF5 数据并保存...", f_log)
    with h5py.File(file_path, 'r') as f:
        # 提取核心数据
        sample_ids = [x.decode('utf-8') if isinstance(x, bytes) else x for x in f['order'][:]]
        coords = f['pcoa_samples'][:]
        eigvals = f['pcoa_eigvals'][:]
        prop_expl = f['pcoa_proportion_explained'][:]
        
    N = len(sample_ids)      # 样本数
    D = coords.shape[1]      # 主成分轴数
    max_rows = max(N, D)
    
    # ------------------ 保存 文件 1 (全数据) ------------------
    df_raw = pd.DataFrame(index=range(max_rows))
    df_raw.loc[0:N-1, 'SampleID'] = sample_ids
    for i in range(D):
        df_raw.loc[0:N-1, f'PC{i+1}'] = coords[:, i]
        
    df_raw.loc[0:D-1, 'Eigenvalues'] = eigvals
    df_raw.loc[0:D-1, 'Proportion_Explained'] = prop_expl
    df_raw.to_csv(out_raw_csv, index=False)
    log_print(f"✅ 文件 1 (全数据) 保存至: {out_raw_csv}", f_log)
    
    # ------------------ 保存 文件 2 (更新：补充特征值和解释度) ------------------
    pc_cols = [f'PC{i+1}' for i in range(D)]
    df_r = pd.DataFrame(coords, index=sample_ids, columns=pc_cols)
    df_r.index.name = 'SampleID'
    
    # === 新增逻辑：将 Eigenvalues 和 Proportion_Explained 补充到 df_r 中 ===
    # 因为 eigvals 只有 D 行（比如 10 行），而 df_r 有 N 行（比如 8 万行）
    # 我们先创建长度为 N 的带有空值 (np.nan) 的数组，再赋给 df_r
    padded_eigvals = np.pad(eigvals, (0, N - len(eigvals)), constant_values=np.nan) if N > len(eigvals) else eigvals[:N]
    padded_prop = np.pad(prop_expl, (0, N - len(prop_expl)), constant_values=np.nan) if N > len(prop_expl) else prop_expl[:N]
    
    df_r['Eigenvalues'] = padded_eigvals
    df_r['Proportion_Explained'] = padded_prop
    
    df_r.to_csv(out_r_csv)
    log_print(f"✅ 文件 2 (R作图格式，已整合 Eigenvalues/Proportion_Explained) 保存至: {out_r_csv}", f_log)
    
    # ==========================================
    # 第二部分：回读与一致性校验 (Sanity Check)
    # ==========================================
    log_print("\n[2/3] 开始重新读取 CSV 并进行一致性校验...", f_log)
    errors = 0
    
    # --- 校验 R_input.csv ---
    log_print("-> 检查 R_input.csv (样本坐标矩阵及附加列):", f_log)
    df_r_check = pd.read_csv(out_r_csv, index_col='SampleID')
    
    # 检查 SampleID 顺序
    if list(df_r_check.index) == sample_ids:
        log_print("   [PASS] 样本 ID 顺序与内容 100% 匹配", f_log)
    else:
        log_print("   [FAIL] 样本 ID 顺序或内容存在差异！", f_log)
        errors += 1
        
    # 检查坐标矩阵数值 (去掉后面新增的两列再做比较)
    coords_check = df_r_check[pc_cols].values
    if np.allclose(coords_check, coords, rtol=1e-05, atol=1e-08):
        log_print("   [PASS] 坐标数值一致 (通过浮点容差测试)", f_log)
    else:
        log_print("   [FAIL] 坐标数值存在明显差异！", f_log)
        errors += 1

    # 校验补充的解释度列
    prop_check_r = df_r_check['Proportion_Explained'].dropna().values
    if np.allclose(prop_check_r, prop_expl, rtol=1e-05, atol=1e-08):
        log_print("   [PASS] 附加列 (Proportion_Explained) 数值一致", f_log)
    else:
        log_print("   [FAIL] 附加列 (Proportion_Explained) 数值存在差异！", f_log)
        errors += 1

    # --- 校验 raw_all.csv ---
    log_print("-> 检查 raw_all.csv (特征值与解释度):", f_log)
    df_raw_check = pd.read_csv(out_raw_csv)
    eigvals_check = df_raw_check['Eigenvalues'].dropna().values
    prop_expl_check = df_raw_check['Proportion_Explained'].dropna().values
    
    if np.allclose(eigvals_check, eigvals, rtol=1e-05, atol=1e-08):
        log_print("   [PASS] 特征值 (Eigenvalues) 数值一致", f_log)
    else:
        log_print("   [FAIL] 特征值存在差异！", f_log)
        errors += 1
        
    if np.allclose(prop_expl_check, prop_expl, rtol=1e-05, atol=1e-08):
        log_print("   [PASS] 解释度 (Proportion_Explained) 数值一致", f_log)
    else:
        log_print("   [FAIL] 解释度存在差异！", f_log)
        errors += 1

    # ==========================================
    # 第三部分：总结
    # ==========================================
    log_print("\n[3/3] 校验结果总结:", f_log)
    if errors == 0:
        log_print("🎉 完美！所有 CSV 导出数据与原始 HDF5 数据完全一致。", f_log)
    else:
        log_print(f"⚠️ 警告！发现 {errors} 处数据不一致，请勿使用导出的文件！", f_log)
        
    # 附赠 R 标签
    r_label_tip = f'\n💡 供 R 绘图直接复制的轴标签参数:\nlabs(x = "PCoA1 ({prop_expl[0]*100:.2f}%)", y = "PCoA2 ({prop_expl[1]*100:.2f}%)")'
    log_print(r_label_tip, f_log)

print(f"\n📝 详细的校验过程已写入日志文件: \n   {log_file}")