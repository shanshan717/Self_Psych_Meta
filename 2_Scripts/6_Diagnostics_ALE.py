# Import necessary modules
from os import makedirs, path
import os
import numpy as np
from nibabel import save
from nilearn import image
from nimare import correct, io, meta, utils
from nimare.diagnostics import FocusCounter, Jackknife
from nimare.reports.base import run_reports
from scipy.stats import norm

def run_ale_with_diagnostics(text_file, voxel_thresh=0.001, cluster_thresh=0.05,
                             random_seed=2025, n_iters=10000, output_dir="../3_Output",
                             run_focus=True, run_jackknife=True, html_report_dir=None):


    # Ensure output directory exists
    makedirs(output_dir, exist_ok=True)
    if random_seed:
        np.random.seed(random_seed)

    # Load data and run ALE
    dset = io.convert_sleuth_to_dataset(text_file=text_file, target="ale_2mm")
    ale = meta.cbma.ALE()
    res = ale.fit(dset)

    # FWE correction
    corr = correct.FWECorrector(method="montecarlo", voxel_thresh=voxel_thresh, n_iters=n_iters)
    cres = corr.transform(res)

    # Save unthresholded and corrected maps
    prefix = path.basename(text_file).replace(".txt", "")
    res.save_maps(output_dir=output_dir, prefix=prefix)
    cres.save_maps(output_dir=output_dir, prefix=prefix)

    # Threshold maps for visualization
    img_clust_mass = cres.get_map("z_desc-mass_level-cluster_corr-FWE_method-montecarlo")
    img_clust_size = cres.get_map("z_desc-size_level-cluster_corr-FWE_method-montecarlo")
    img_z = cres.get_map("z")
    img_ale = cres.get_map("stat")
    
    cluster_thresh_z = norm.ppf(1 - cluster_thresh / 2)
    img_clust_mass_thresh = image.threshold_img(img=img_clust_mass, threshold=cluster_thresh_z)
    img_clust_size_thresh = image.threshold_img(img=img_clust_size, threshold=cluster_thresh_z)
    
    img_mask = image.math_img("np.where(img > 0, 1, 0)", img=img_clust_size_thresh)
    img_z_thresh = image.math_img("img1 * img2", img1=img_mask, img2=img_z)
    img_ale_thresh = image.math_img("img1 * img2", img1=img_mask, img2=img_ale)

    # Save thresholded maps
    save(img=img_clust_mass_thresh, filename=f"{output_dir}/{prefix}_z_mass_level_thresh.nii.gz")
    save(img=img_clust_size_thresh, filename=f"{output_dir}/{prefix}_z_size_level_thresh.nii.gz")
    save(img=img_z_thresh, filename=f"{output_dir}/{prefix}_z_thresh.nii.gz")
    save(img=img_ale_thresh, filename=f"{output_dir}/{prefix}_stat_size_thresh.nii.gz")

    # Run diagnostics 
    diag_results = {}
    
    if run_focus:
        counter = FocusCounter(
            target_image="z_desc-size_level-cluster_corr-FWE_method-montecarlo",
            voxel_thresh=None
        )
        focus_res = counter.transform(cres)
        diag_results['focuscounter'] = focus_res
        # 输出表格示例
        # 保存所有 FocusCounter 表（最稳健）
        focus_dir = path.join(output_dir, "focuscounter")
        os.makedirs(focus_dir, exist_ok=True)

        for name, df in focus_res.tables.items():
            safe_name = name.replace(":", "_")
            csv_path = path.join(focus_dir,
                                f"{prefix}_focuscounter_{safe_name}.csv")
            df.to_csv(csv_path, index=False)
            print(f"Saved: {csv_path}")

    # Generate HTML report if path is provided
    if html_report_dir:
        makedirs(html_report_dir, exist_ok=True)
        
        
        try:
            run_reports(cres, html_report_dir)
            print(f"\nHTML report saved to {html_report_dir}")
        except Exception as e:
            import traceback
            print(f"\nReport generation failed: {e}")
            traceback.print_exc()

    return cres, diag_results

if __name__ == "__main__":
    text_file = "../1_Data/AnalysisData/InputData_ALE/hyperactivation.txt"
    file_prefix = path.basename(text_file).replace(".txt", "")
    output_dir = "../3_Output/7_other"
    html_report_path = path.join(output_dir, f"{file_prefix}_html_report")

    cres, diag_results = run_ale_with_diagnostics(
        text_file=text_file,
        voxel_thresh=0.001,
        cluster_thresh=0.05,
        random_seed=2025,
        n_iters=10000,
        output_dir=output_dir,
        run_focus=True,
        html_report_dir=html_report_path
    )