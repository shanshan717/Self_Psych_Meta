#!/usr/bin/env python3
"""Run robustness analyses (FSN + FocusCounter) for ALE Sleuth datasets.

Default targets follow project naming:
- healthy-patient  -> hyperactivation
- patient-healthy  -> hypoactivation
- healthy+patients -> hyper_and_hypo

Contrast analyses are intentionally excluded.
"""

from __future__ import annotations

import argparse
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
from glob import glob
from math import sqrt
from os import makedirs, path
from re import sub
from shutil import copy
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from nibabel import save
from nilearn import image, plotting, reporting
from nimare import correct, io, meta, utils
from nimare.diagnostics import FocusCounter
from scipy.stats import norm

TARGET_ALIASES: Dict[str, str] = {
    "healthy-patient": "hyperactivation",
    "healthy_patient": "hyperactivation",
    "control-patient": "hyperactivation",
    "control_patient": "hyperactivation",
    "patient-healthy": "hypoactivation",
    "patient_healthy": "hypoactivation",
    "patient-control": "hypoactivation",
    "patient_control": "hypoactivation",
    "healthy+patients": "hyper_and_hypo",
    "healthy_patients": "hyper_and_hypo",
    "control+patient": "hyper_and_hypo",
    "hyperactivation": "hyperactivation",
    "hypoactivation": "hypoactivation",
    "hyper_and_hypo": "hyper_and_hypo",
}

DEFAULT_TARGETS = ["healthy-patient", "patient-healthy", "healthy+patients"]


def resolve_targets(targets: Iterable[str]) -> List[str]:
    resolved = []
    for t in targets:
        key = t.strip().lower()
        if key not in TARGET_ALIASES:
            raise ValueError(f"Unsupported target '{t}'.")
        canon = TARGET_ALIASES[key]
        if canon not in resolved:
            resolved.append(canon)
    return resolved


def generate_null_dataset(
    text_file: str,
    space: str,
    k_null: int,
    random_seed: int,
    output_dir: str,
):
    temp = utils.get_template(space=space, mask="brain")
    x, y, z = np.where(temp.get_fdata() == 1.0)
    within_mni = image.coord_transform(x=x, y=y, z=z, affine=temp.affine)
    within_mni = np.array(within_mni).transpose()

    dset = io.convert_sleuth_to_dataset(text_file=text_file, target=space)

    random.seed(random_seed)

    nr_subjects_dset = [n[0] for n in dset.metadata["sample_sizes"]]
    nr_subjects_null = random.choices(nr_subjects_dset, k=k_null)

    nr_peaks_dset = dset.coordinates["study_id"].value_counts().tolist()
    nr_peaks_null = random.choices(nr_peaks_dset, k=k_null)

    idx_list = [
        random.sample(range(len(within_mni)), k=k_peaks) for k_peaks in nr_peaks_null
    ]
    peaks_null = [within_mni[idx] for idx in idx_list]

    makedirs(output_dir, exist_ok=True)
    text_file_basename = path.basename(text_file)
    null_file_basename = sub(
        pattern=".txt",
        repl="_plus_k" + str(k_null) + ".txt",
        string=text_file_basename,
    )
    null_file = path.join(output_dir, null_file_basename)
    copy(text_file, null_file)

    with open(null_file, mode="a", encoding="utf-8") as f:
        for i in range(k_null):
            f.write(
                "\n// nullstudy"
                + str(i + 1)
                + "\n// Subjects="
                + str(nr_subjects_null[i])
                + "\n"
            )
            np.savetxt(f, peaks_null[i], fmt="%.3f", delimiter="\t")

    dset_null = io.convert_sleuth_to_dataset(null_file, target=space)
    return dset_null


def compute_single_fsn(
    text_file: str,
    output_dir: str,
    random_null_seed: int,
    space: str,
    voxel_thresh: float,
    cluster_thresh: float,
    n_iters: int,
    k_max_factor: int,
    random_ale_seed: int,
) -> str:
    np.random.seed(random_ale_seed)

    ale = meta.cbma.ALE()
    corr = correct.FWECorrector(
        method="montecarlo", voxel_thresh=voxel_thresh, n_iters=n_iters
    )

    dset_orig = io.convert_sleuth_to_dataset(text_file=text_file, target=space)
    res_orig = ale.fit(dset_orig)
    cres_orig = corr.transform(res_orig)

    ids_orig = dset_orig.ids.tolist()
    k_max = len(ids_orig) * k_max_factor

    dset_null = generate_null_dataset(
        text_file=text_file,
        space=space,
        k_null=k_max,
        random_seed=random_null_seed,
        output_dir=output_dir,
    )

    img_fsn = cres_orig.get_map("z_desc-size_level-cluster_corr-FWE_method-montecarlo")
    cluster_thresh_z = norm.ppf(1 - cluster_thresh / 2)
    img_fsn = image.threshold_img(img_fsn, threshold=cluster_thresh_z)
    img_fsn = image.math_img("np.where(img > 0, 1, 0)", img=img_fsn)

    img_z = cres_orig.get_map("z")
    img_z = image.math_img("img1 * img2", img1=img_fsn, img2=img_z)

    for k in range(1, k_max):
        ids_null = [f"nullstudy{x}-" for x in range(1, k + 1)]
        ids = ids_orig + ids_null
        dset_k = dset_null.slice(ids)

        res_k = ale.fit(dset_k)
        cres_k = corr.transform(result=res_k)

        img_k = cres_k.get_map("z_desc-size_level-cluster_corr-FWE_method-montecarlo")
        img_k = image.threshold_img(img_k, threshold=cluster_thresh_z)
        img_k = image.math_img("np.where(img > 0, 1, 0)", img=img_k)

        count = str(k + 1)
        formula = "np.where(img_fsn + img_k == " + count + ", img_fsn + 1, img_fsn)"
        img_fsn = image.math_img(formula, img_fsn=img_fsn, img_k=img_k)

        if not np.any(img_k.get_fdata()):
            break

    prefix = path.basename(text_file).replace(".txt", "")
    filename_img = path.join(output_dir, f"{prefix}_fsn.nii.gz")
    save(img_fsn, filename=filename_img)

    tab_fsn = reporting.get_clusters_table(img_z, stat_threshold=0, min_distance=1000)
    inv_affine = np.linalg.inv(img_z.affine)
    x, y, z = [np.array(tab_fsn[col]) for col in ["X", "Y", "Z"]]
    x, y, z = image.coord_transform(x=x, y=y, z=z, affine=inv_affine)
    x, y, z = [arr.astype("int") for arr in [x, y, z]]
    tab_fsn["FSN"] = img_fsn.get_fdata()[x, y, z]

    filename_tab = path.join(output_dir, f"{prefix}_fsn.tsv")
    tab_fsn.to_csv(filename_tab, sep="\t", index=False)

    return output_dir


def aggregate_fsn(prefix: str, out_dir: str, ci_level: float = 0.05) -> None:
    fnames_maps = glob(path.join(out_dir, "filedrawer*", f"{prefix}_fsn.nii.gz"))
    if not fnames_maps:
        raise RuntimeError(f"No FSN maps found for {prefix} in {out_dir}")

    imgs_fsn = [image.load_img(fname) for fname in fnames_maps]
    img_mean = image.mean_img(imgs_fsn)
    fname_img_mean = path.join(out_dir, f"{prefix}_mean_fsn.nii.gz")
    save(img_mean, fname_img_mean)

    fnames_tabs = glob(path.join(out_dir, "filedrawer*", f"{prefix}_fsn.tsv"))
    tabs_fsn = [pd.read_csv(fname, delimiter="\t") for fname in fnames_tabs]
    tab_fsn = pd.concat(tabs_fsn)

    agg = tab_fsn.groupby("Cluster ID")["FSN"].agg(["mean", "count", "std"])
    z_crit = abs(norm.ppf(ci_level / 2))
    agg["se"] = [std / sqrt(count) if count > 0 else 0 for std, count in zip(agg["std"], agg["count"])]
    agg["ci_lower"] = agg["mean"] - z_crit * agg["se"]
    agg["ci_upper"] = agg["mean"] + z_crit * agg["se"]

    fname_agg = path.join(out_dir, f"{prefix}_mean_fsn.csv")
    agg.to_csv(fname_agg, float_format="%.3f")


def save_fsn_plot(prefix: str, out_dir: str, n_studies: int, ratio: float) -> None:
    img_fsn = image.load_img(path.join(out_dir, f"{prefix}_mean_fsn.nii.gz"))
    data = img_fsn.get_fdata()

    fsn_threshold = max(1, ratio * n_studies)

    img_below = image.new_img_like(
        img_fsn,
        ((data > 0) & (data < fsn_threshold)).astype(float),
    )
    img_above = image.new_img_like(
        img_fsn,
        (data >= fsn_threshold).astype(float),
    )

    disp = plotting.plot_glass_brain(None, display_mode="lyrz")
    disp.add_overlay(img_below, cmap="OrRd", vmin=0, vmax=1, colorbar=False)
    disp.add_overlay(img_above, cmap="RdYlGn", vmin=0, vmax=1, colorbar=False)
    disp.savefig(path.join(out_dir, f"{prefix}_fsn_thresholded_glass_brain.png"))
    plt.close("all")

    disp2 = plotting.plot_glass_brain(
        img_fsn,
        display_mode="lyrz",
        colorbar=True,
        cmap="RdYlGn",
        vmin=0,
        vmax=max(1, np.nanmax(data)),
    )
    disp2.savefig(path.join(out_dir, f"{prefix}_fsn_continuous_glass_brain.png"))
    plt.close("all")


def _pick_focus_table(focus_tables: Dict[str, pd.DataFrame]) -> Tuple[str, pd.DataFrame]:
    for key in focus_tables:
        if "tab-clust" in key or "clust" in key:
            return key, focus_tables[key]
    first = next(iter(focus_tables))
    return first, focus_tables[first]


def _pick_count_column(df: pd.DataFrame) -> str:
    for c in df.columns:
        if "count" in c.lower() and pd.api.types.is_numeric_dtype(df[c]):
            return c
    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if not numeric_cols:
        raise ValueError("No numeric columns in FocusCounter table for plotting.")
    return numeric_cols[0]


def run_focuscounter(
    text_file: str,
    output_dir: str,
    space: str,
    voxel_thresh: float,
    cluster_thresh: float,
    n_iters: int,
    random_seed: int,
) -> None:
    np.random.seed(random_seed)
    makedirs(output_dir, exist_ok=True)

    prefix = path.basename(text_file).replace(".txt", "")

    dset = io.convert_sleuth_to_dataset(text_file=text_file, target=space)
    ale = meta.cbma.ALE()
    res = ale.fit(dset)

    corr = correct.FWECorrector(method="montecarlo", voxel_thresh=voxel_thresh, n_iters=n_iters)
    cres = corr.transform(res)

    map_dir = path.join(output_dir, "maps")
    makedirs(map_dir, exist_ok=True)
    cres.save_maps(output_dir=map_dir, prefix=prefix)

    counter = FocusCounter(
        target_image="z_desc-size_level-cluster_corr-FWE_method-montecarlo",
        voxel_thresh=None,
    )
    focus_res = counter.transform(cres)

    table_dir = path.join(output_dir, "tables")
    makedirs(table_dir, exist_ok=True)
    for name, df in focus_res.tables.items():
        safe_name = name.replace(":", "_")
        csv_path = path.join(table_dir, f"{prefix}_focuscounter_{safe_name}.csv")
        df.to_csv(csv_path, index=False)

    # ALE thresholded glass brain (for direct robustness visual check)
    z_size = cres.get_map("z_desc-size_level-cluster_corr-FWE_method-montecarlo")
    cluster_thresh_z = norm.ppf(1 - cluster_thresh / 2)
    z_thr = image.threshold_img(z_size, threshold=cluster_thresh_z)
    disp = plotting.plot_glass_brain(z_thr, display_mode="lyrz", colorbar=True)
    disp.savefig(path.join(output_dir, f"{prefix}_focuscounter_ale_glass_brain.png"))
    plt.close("all")

    # FocusCounter bar plot
    key, df = _pick_focus_table(focus_res.tables)
    if not df.empty:
        count_col = _pick_count_column(df)
        label_col = "cluster_id" if "cluster_id" in df.columns else df.columns[0]

        vis_df = df[[label_col, count_col]].copy()
        vis_df = vis_df.sort_values(count_col, ascending=False).head(20)

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(vis_df[label_col].astype(str), vis_df[count_col])
        ax.invert_yaxis()
        ax.set_xlabel(count_col)
        ax.set_title(f"{prefix} FocusCounter ({key})")
        fig.tight_layout()
        fig.savefig(path.join(output_dir, f"{prefix}_focuscounter_top20.png"), dpi=200)
        plt.close(fig)


def run_target(
    target: str,
    args: argparse.Namespace,
) -> None:
    text_file = path.join(args.input_dir, f"{target}.txt")
    if not path.exists(text_file):
        raise FileNotFoundError(f"Input file not found: {text_file}")

    base_out = path.join(args.output_root, target)
    fsn_out = path.join(base_out, "fsn")
    focus_out = path.join(base_out, "focuscounter")
    makedirs(fsn_out, exist_ok=True)
    makedirs(focus_out, exist_ok=True)

    dset = io.convert_sleuth_to_dataset(text_file=text_file, target=args.space)
    n_studies = len(dset.ids)

    if args.run_fsn:
        random.seed(args.random_null_seed)
        seeds = random.sample(range(100000), k=args.nr_filedrawers)

        jobs = []
        for seed in seeds:
            filedrawer = f"filedrawer{seed}"
            out_dir = path.join(fsn_out, filedrawer)
            jobs.append((
                text_file,
                out_dir,
                seed,
                args.space,
                args.fsn_voxel_thresh,
                args.cluster_thresh,
                args.fsn_n_iters,
                args.k_max_factor,
                args.random_ale_seed,
            ))

        if args.fsn_workers > 1 and len(jobs) > 1:
            with ProcessPoolExecutor(max_workers=args.fsn_workers) as ex:
                futures = [ex.submit(compute_single_fsn, *job) for job in jobs]
                for fut in as_completed(futures):
                    fut.result()
        else:
            for job in jobs:
                compute_single_fsn(*job)

        aggregate_fsn(prefix=target, out_dir=fsn_out)
        save_fsn_plot(
            prefix=target,
            out_dir=fsn_out,
            n_studies=n_studies,
            ratio=args.fsn_threshold_ratio,
        )

    if args.run_focuscounter:
        run_focuscounter(
            text_file=text_file,
            output_dir=focus_out,
            space=args.space,
            voxel_thresh=args.focus_voxel_thresh,
            cluster_thresh=args.cluster_thresh,
            n_iters=args.focus_n_iters,
            random_seed=args.random_ale_seed,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Robustness analysis for ALE: FSN + FocusCounter"
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=DEFAULT_TARGETS,
        help="Targets to run: healthy-patient patient-healthy healthy+patients (or canonical names)",
    )
    parser.add_argument(
        "--input-dir",
        default="../1_Data/AnalysisData/InputData_ALE",
        help="Directory with Sleuth .txt input files",
    )
    parser.add_argument(
        "--output-root",
        default="../3_Output/9_Robustness",
        help="Root output directory",
    )
    parser.add_argument("--space", default="ale_2mm")

    parser.add_argument("--run-fsn", action="store_true")
    parser.add_argument("--run-focuscounter", action="store_true")

    parser.add_argument("--nr-filedrawers", type=int, default=5)
    parser.add_argument("--k-max-factor", type=int, default=5)
    parser.add_argument("--fsn-workers", type=int, default=1)

    parser.add_argument("--cluster-thresh", type=float, default=0.05)
    parser.add_argument("--fsn-voxel-thresh", type=float, default=0.001)
    parser.add_argument("--focus-voxel-thresh", type=float, default=0.001)

    parser.add_argument("--fsn-n-iters", type=int, default=100)
    parser.add_argument("--focus-n-iters", type=int, default=10000)

    parser.add_argument("--random-ale-seed", type=int, default=2025)
    parser.add_argument("--random-null-seed", type=int, default=2026)

    parser.add_argument(
        "--fsn-threshold-ratio",
        type=float,
        default=0.3,
        help="Threshold ratio for FSN visualization (ratio * N studies)",
    )

    args = parser.parse_args()

    if not args.run_fsn and not args.run_focuscounter:
        args.run_fsn = True
        args.run_focuscounter = True

    args.targets = resolve_targets(args.targets)
    return args


def main() -> None:
    args = parse_args()
    for target in args.targets:
        print(f"\n=== Running robustness analyses for {target} ===")
        run_target(target, args)
        print(f"=== Completed {target} ===")


if __name__ == "__main__":
    main()
