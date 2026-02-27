from os import makedirs, path
import os
import numpy as np
import pandas as pd
from nibabel import save
from nilearn import image
from nimare import io, meta, correct
from nimare.diagnostics import FocusCounter
from nimare.workflows import PairwiseCBMAWorkflow
from nimare.reports.base import run_reports

def run_subtraction_with_diagnostics(
    text_file1,
    text_file2,
    voxel_thresh=0.001,
    n_iters=10000,
    output_dir="../Output/2_subtraction",
    random_seed=2025
):
    print(f'SUBTRACTION ANALYSIS: {path.basename(text_file1)} vs {path.basename(text_file2)}')
    
    # 确保输出目录存在
    makedirs(output_dir, exist_ok=True)
    if random_seed:
        np.random.seed(random_seed)

    # 1. 加载数据 (Sleuth text files)
    # 注意：确保你的环境能支持 ale_2mm 模板，否则 NiMARE 会自动下载
    dset1 = io.convert_sleuth_to_dataset(text_file=text_file1, target="ale_2mm")
    dset2 = io.convert_sleuth_to_dataset(text_file=text_file2, target="ale_2mm")

    # 2. 配置 Estimator (Subtraction)
    # 关键点：n_iters 和 voxel_thresh 必须在这里定义
    sub_estimator = meta.cbma.ALESubtraction(n_iters=n_iters, voxel_thresh=voxel_thresh)
    
    # 3. 配置 Corrector (FWE)
    # 关键修复：这里不要传 n_iters 或 voxel_thresh，它会自动从 Estimator 继承
    sub_corrector = correct.FWECorrector(method="montecarlo")

    # 4. 配置 Workflow 并集成 FocusCounter
    # display_second_group=True: 在表中显示两个组的贡献情况
    workflow = PairwiseCBMAWorkflow(
        estimator=sub_estimator,
        corrector=sub_corrector,
        diagnostics=FocusCounter(voxel_thresh=voxel_thresh, display_second_group=True)
    )

    # 5. 运行分析
    print(f"Fitting workflow with {n_iters} iterations (this may take a while)...")
    cres = workflow.fit(dset1, dset2)

    # 6. 保存 NIfTI 结果 (Z maps, etc.)
    name1 = path.basename(text_file1).replace(".txt", "")
    name2 = path.basename(text_file2).replace(".txt", "")
    prefix = f"{name1}_minus_{name2}"
    
    # save_maps 会自动保存未校正和校正后的图
    cres.save_maps(output_dir=output_dir, prefix=prefix)
    print(f"NIfTI maps saved with prefix: {prefix}")

    # 7. 提取并保存 FocusCounter CSV 结果 (双向导出)
    # 循环查找 'positive' (Group1 > Group2) 和 'negative' (Group2 > Group1) 的表
    found_tables = False
    for tail in ["positive", "negative"]:
        # 在结果字典中查找符合条件的键名
        diag_table_key = [k for k in cres.tables.keys() if "FocusCounter" in k and "counts" in k and tail in k]
        
        if diag_table_key:
            # 提取 DataFrame
            df_focus = cres.tables[diag_table_key[0]]
            
            # 构造文件名并保存
            csv_path = path.join(output_dir, f"{prefix}_focuscounter_{tail}.csv")
            df_focus.to_csv(csv_path, index=True)
            print(f"FocusCounter table ({tail}) saved to: {csv_path}")
            found_tables = True
    
    if not found_tables:
        print("No significant clusters found for FocusCounter tables (or n_iters was too low).")

    # 8. 生成 HTML 报告
    html_dir = path.join(output_dir, f"{prefix}_report")
    makedirs(html_dir, exist_ok=True)
    try:
        run_reports(cres, html_dir)
        print(f"HTML report generated at: {html_dir}")
    except Exception as e:
        print(f"Report generation failed: {e}")

    return cres

# --- 主程序运行入口 ---
if __name__ == "__main__":
    # 请根据实际情况修改文件路径
    file1 = "../1_Data/AnalysisData/InputData_ALE/patient.txt"
    file2 = "../1_Data/AnalysisData/InputData_ALE/Control_all.txt"
    
    # 检查文件是否存在
    if not path.exists(file1) or not path.exists(file2):
        print("Error: Input files not found. Please check the file paths.")
    else:
        run_subtraction_with_diagnostics(
            text_file1=file1,
            text_file2=file2,
            voxel_thresh=0.001,
            n_iters=100, # 正式跑建议 10000，测试可改小 (如 100)
            output_dir="../Output/2_subtraction"
        )