#!/usr/bin/env python
"""
Singe_cell_fret.py

Compute per-track FRET time traces (A/D) from:
- Cellpose mask TIFFs
- Donor channel TIFF stack
- Acceptor channel TIFF stack
- Track CSV from the sticky tracker (track_render_overlay.py)

Outputs:
- CSV with per-frame, per-track FRET ratios and intensities
- Optional Excel with one sheet per track
- Optional PNG plots of FRET vs time

Requires:
    pip install numpy pandas tifffile scipy scikit-image matplotlib
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter
from skimage.filters import threshold_sauvola
from skimage.morphology import dilation, erosion, disk
import matplotlib.pyplot as plt


# ==================== DEFAULTS ====================

DEFAULTS = {
    "min_frames": 45,
    "frame_interval_min": 10.0,       # minutes between frames
    "ring_dilate": 12,
    "ring_thick": 2,
    "use_sauvola": True,
    "sigma_smooth": 1.0,
    "sauv_window": 551,
    "sauv_k": 0.15,
    "erode_pixels": 1,
    "alpha_bleed": 0.0,
    "beta_bleed": 0.0,
    "write_excel": True,
    "save_plots": True,
    "plot_tag": "ensemble_like_hours",
}

__SHAPE_WARNED = False


# ==================== Helpers ====================

def _safe_excel_writer(path, write_excel=True):
    if not write_excel:
        return None
    try:
        return pd.ExcelWriter(path, engine="xlsxwriter")
    except Exception:
        try:
            return pd.ExcelWriter(path, engine="openpyxl")
        except Exception:
            return None


def _save_df_excel_or_csv(df, xlsx_path, csv_path, sheet_name="Sheet1", write_excel=True):
    w = _safe_excel_writer(xlsx_path, write_excel=write_excel)
    if w is None:
        df.to_csv(csv_path, index=False)
        return ("csv", csv_path)
    with w:
        df.to_excel(w, index=False, sheet_name=sheet_name)
    return ("xlsx", xlsx_path)


def _to_2d(a):
    a = np.asarray(a)
    if a.ndim == 2:
        return a
    if a.ndim > 2:
        a = np.squeeze(a)
        if a.ndim != 2:
            a = a.reshape((a.shape[0], a.shape[1]))
    return a


def load_mask(path):
    m = tifffile.imread(path)
    m = _to_2d(m).astype(np.int32, copy=False)
    return m


def load_frame(path):
    im = tifffile.imread(path)
    im = _to_2d(im).astype(np.float32, copy=False)
    return im


def align_mask_to_image(mask_bool, image_2d):
    global __SHAPE_WARNED
    H, W = image_2d.shape[:2]
    mh, mw = mask_bool.shape[:2]
    if (mh, mw) == (H, W):
        return mask_bool
    m = mask_bool[:H, :W]
    py, px = H - m.shape[0], W - m.shape[1]
    if py > 0 or px > 0:
        m = np.pad(
            m,
            ((0, max(0, py)), (0, max(0, px))),
            mode="constant",
            constant_values=False,
        )
    if not __SHAPE_WARNED:
        print(
            f"[WARN] Aligned mask {mh}x{mw} → image {H}x{W}. "
            f"Upstream frames likely cropped/resized differently."
        )
        __SHAPE_WARNED = True
    return m


def ring_mask(mask_bool, ring_dilate=8, ring_thick=4):
    outer = dilation(mask_bool, disk(ring_dilate))
    inner = dilation(mask_bool, disk(max(1, ring_dilate - ring_thick)))
    return outer ^ inner


# ==================== Core function ====================

def compute_fret_traces(
    mask_dir,
    donor_dir,
    acceptor_dir,
    tracks_csv=None,
    min_frames=DEFAULTS["min_frames"],
    frame_interval_min=DEFAULTS["frame_interval_min"],
    ring_dilate=DEFAULTS["ring_dilate"],
    ring_thick=DEFAULTS["ring_thick"],
    use_sauvola=DEFAULTS["use_sauvola"],
    sigma_smooth=DEFAULTS["sigma_smooth"],
    sauv_window=DEFAULTS["sauv_window"],
    sauv_k=DEFAULTS["sauv_k"],
    erode_pixels=DEFAULTS["erode_pixels"],
    alpha_bleed=DEFAULTS["alpha_bleed"],
    beta_bleed=DEFAULTS["beta_bleed"],
    write_excel=DEFAULTS["write_excel"],
    save_plots=DEFAULTS["save_plots"],
    plot_tag=DEFAULTS["plot_tag"],
):
    """
    Main FRET extraction pipeline.

    Parameters mirror the original CONFIG block but are now passed as arguments.
    """

    # ---------- 0) Resolve paths & load tracks ----------

    if tracks_csv is None:
        tracks_csv = os.path.join(mask_dir, "tracks_v3.7_sticky_new.csv")

    if not os.path.exists(tracks_csv):
        raise FileNotFoundError(f"Missing tracks CSV: {tracks_csv}")

    tracks_raw = pd.read_csv(tracks_csv)

    summary = (
        tracks_raw.groupby("track_id")
        .agg(
            frames_seen=("frame", "nunique"),
            first_frame=("frame", "min"),
            last_frame=("frame", "max"),
            n_rows=("frame", "count"),
            mean_area=("area", "mean"),
        )
        .reset_index()
        .sort_values("frames_seen", ascending=False)
    )

    keep_ids = set(
        summary.loc[summary["frames_seen"] >= min_frames, "track_id"].astype(int)
    )

    summ_xlsx = os.path.join(mask_dir, f"tracks_ge{min_frames}.xlsx")
    summ_csv = os.path.join(mask_dir, f"tracks_ge{min_frames}.csv")
    kind, path = _save_df_excel_or_csv(
        summary[summary["track_id"].isin(keep_ids)],
        summ_xlsx,
        summ_csv,
        sheet_name=f"tracks_ge{min_frames}",
        write_excel=write_excel,
    )
    print(
        f"[INFO] Persistent tracks (>= {min_frames} frames) → {path} ({kind}), "
        f"n={len(keep_ids)}"
    )
    if not keep_ids:
        print("[WARN] No tracks meet min_frames; exiting.")
        raise SystemExit(0)

    # ---------- 1) Gather files and sort by index ----------

    mask_files = sorted(glob.glob(os.path.join(mask_dir, "*cp_masks.tif")))
    if not mask_files:
        raise FileNotFoundError(f"No '*cp_masks.tif' under {mask_dir}")

    donor_files = sorted(glob.glob(os.path.join(donor_dir, "*.tif")))
    accept_files = sorted(glob.glob(os.path.join(acceptor_dir, "*.tif")))

    print(
        f"[INFO] Counts → masks:{len(mask_files)} "
        f"donor:{len(donor_files)} acceptor:{len(accept_files)}"
    )

    max_frame_in_csv = int(tracks_raw["frame"].max())
    min_required = max_frame_in_csv + 1
    if len(mask_files) < min_required:
        print(
            f"[WARN] You have frames up to {max_frame_in_csv} in CSV, "
            f"but only {len(mask_files)} mask files."
        )
    if len(donor_files) < min_required or len(accept_files) < min_required:
        print(
            "[WARN] Donor/Acceptor sequences shorter than CSV frames "
            f"({len(donor_files)}/{len(accept_files)} vs need ≥ {min_required})."
        )

    tracks = tracks_raw[tracks_raw["track_id"].isin(keep_ids)].copy()
    tracks = tracks.sort_values(["frame", "track_id"])

    rows = []
    missing_examples = []

    # ---------- 2) Loop frames (index-based pairing) ----------

    frames_to_process = sorted(tracks["frame"].unique())
    for f_idx in frames_to_process:
        if (
            f_idx >= len(mask_files)
            or f_idx >= len(donor_files)
            or f_idx >= len(accept_files)
        ):
            if len(missing_examples) < 8:
                missing_examples.append(f_idx)
            continue

        mpath = mask_files[f_idx]
        dpath = donor_files[f_idx]
        apath = accept_files[f_idx]

        mask = load_mask(mpath)
        donor = load_frame(dpath)
        accept = load_frame(apath)

        global __SHAPE_WARNED
        if donor.shape != accept.shape and not __SHAPE_WARNED:
            print(
                f"[WARN] Donor {donor.shape} and acceptor {accept.shape} differ; "
                f"using donor as reference shape."
            )
            __SHAPE_WARNED = True

        if accept.shape != donor.shape:
            H, W = donor.shape
            a = accept[:H, :W]
            py, px = H - a.shape[0], W - a.shape[1]
            if py > 0 or px > 0:
                a = np.pad(
                    a,
                    ((0, max(0, py)), (0, max(0, px))),
                    mode="edge",
                )
            accept_use = a
        else:
            accept_use = accept

        dfF = tracks[tracks["frame"] == f_idx]
        if dfF.empty:
            continue

        if use_sauvola:
            donor_s = gaussian_filter(donor, sigma_smooth)
            accept_s = gaussian_filter(accept_use, sigma_smooth)
            th = threshold_sauvola(accept_s, window_size=sauv_window, k=sauv_k)
            fg = accept_s > th
        else:
            donor_s, accept_s, fg = donor, accept_use, None

        for _, r in dfF.iterrows():
            lbl = int(r["label"])
            tid = int(r["track_id"])

            m = mask == lbl
            if m.sum() == 0:
                continue

            m = align_mask_to_image(m.astype(bool, copy=False), donor)

            ring = ring_mask(m, ring_dilate, ring_thick)

            if use_sauvola:
                m_core = erosion(m, disk(max(1, erode_pixels)))
                mc = m_core & fg
                if not mc.any():
                    mc = m_core
                useD, useA = donor_s, accept_s
            else:
                mc = m
                useD, useA = donor, accept_use

            if ring.any():
                D_bg = float(np.median(useD[ring]))
                A_bg = float(np.median(useA[ring]))
            else:
                D_bg = 0.0
                A_bg = 0.0

            if mc.any():
                D_fg = float(np.mean(useD[mc]))
                A_fg = float(np.mean(useA[mc]))
            else:
                D_fg = 0.0
                A_fg = 0.0

            D_corr = max(D_fg - D_bg, 1e-9)
            A_corr = max(A_fg - A_bg, 1e-9)

            A_bt = max(A_corr - alpha_bleed * D_corr, 1e-9)
            D_bt = max(D_corr - beta_bleed * A_corr, 1e-9)

            fret_ratio = A_bt / D_bt

            time_hr = (f_idx * frame_interval_min) / 60.0

            rows.append(
                {
                    "track_id": tid,
                    "frame": f_idx,
                    "time_hr": time_hr,
                    "donor_mean": D_fg,
                    "acceptor_mean": A_fg,
                    "donor_bg": D_bg,
                    "acceptor_bg": A_bg,
                    "donor_corr": D_corr,
                    "acceptor_corr": A_corr,
                    "fret_ratio": fret_ratio,
                }
            )

    if missing_examples:
        print(
            f"[WARN] Skipped {len(missing_examples)} frame(s) because one or more "
            f"files were missing by index."
        )
        print(f"       First few missing indices: {missing_examples[:8]}")

    # ---------- 3) Save numeric outputs ----------

    if not rows:
        raise RuntimeError(
            "No FRET rows were produced. "
            "Likely file pairing by index failed—check that donor/acceptor/masks "
            "are the same length and sorted so that index i is the same frame in each."
        )

    out = pd.DataFrame(rows).sort_values(["track_id", "frame"])

    tag = f"ge{min_frames}_{plot_tag}"
    fret_csv = os.path.join(mask_dir, f"fret_traces_{tag}.csv")
    out.to_csv(fret_csv, index=False)
    print(f"[INFO] Wrote time traces → {fret_csv}")

    if write_excel:
        fret_xlsx = os.path.join(mask_dir, f"fret_traces_{tag}.xlsx")
        w = _safe_excel_writer(fret_xlsx, write_excel=True)
        if w is None:
            print("[INFO] Excel engine not available; CSV only.")
        else:
            with w:
                for tid, g in out.groupby("track_id"):
                    g.to_excel(w, index=False, sheet_name=f"track_{tid}")
            print(f"[INFO] Also wrote Excel → {fret_xlsx}")

    # ---------- 4) Plotting ----------

    if save_plots:
        plot_dir = os.path.join(mask_dir, f"fret_trace_plots_{tag}")
        os.makedirs(plot_dir, exist_ok=True)
        for tid, g in out.groupby("track_id"):
            plt.figure()
            plt.plot(g["time_hr"], g["fret_ratio"], marker="o", linewidth=1)
            plt.xlabel("Time (hours)")
            plt.ylabel("FRET ratio (A/D)")
            plt.title(f"Track {tid} (n={len(g)} frames)")
            plt.tight_layout()
            plt.savefig(
                os.path.join(plot_dir, f"track_{tid}_{plot_tag}.png"),
                dpi=150,
            )
            plt.close()
        print(f"[INFO] Plots saved in {plot_dir}")

    print("[DONE]")
    return out, fret_csv


# ==================== CLI ====================

def build_argparser():
    p = argparse.ArgumentParser(
        description="Extract ensemble-like single-cell FRET traces from donor/acceptor images and sticky tracks."
    )
    p.add_argument(
        "--mask-dir",
        required=True,
        help="Directory containing Cellpose masks (*cp_masks.tif).",
    )
    p.add_argument(
        "--donor-dir",
        required=True,
        help="Directory containing donor channel TIFFs (one per frame).",
    )
    p.add_argument(
        "--acceptor-dir",
        required=True,
        help="Directory containing acceptor channel TIFFs (one per frame).",
    )
    p.add_argument(
        "--tracks-csv",
        default=None,
        help="Path to sticky tracks CSV (default: <mask-dir>/tracks_v3.7_sticky_new.csv).",
    )

    # core filters
    p.add_argument("--min-frames", type=int, help="Min frames to keep a track.")
    p.add_argument(
        "--frame-interval-min",
        type=float,
        help="Time between frames in minutes (for time axis).",
    )

    # ring background
    p.add_argument("--ring-dilate", type=int, help="Ring dilation radius.")
    p.add_argument("--ring-thick", type=int, help="Ring thickness in pixels.")

    # gating
    p.add_argument(
        "--no-sauvola",
        action="store_true",
        help="Disable Sauvola-based per-cell gating.",
    )
    p.add_argument("--sigma-smooth", type=float, help="Gaussian smoothing sigma.")
    p.add_argument("--sauv-window", type=int, help="Sauvola window size.")
    p.add_argument("--sauv-k", type=float, help="Sauvola k parameter.")
    p.add_argument("--erode-pixels", type=int, help="Pixels to erode core mask.")

    # bleed-through
    p.add_argument(
        "--alpha",
        type=float,
        help="Donor→acceptor bleed-through factor (ALPHA).",
    )
    p.add_argument(
        "--beta",
        type=float,
        help="Acceptor→donor bleed-through factor (BETA).",
    )

    # IO
    p.add_argument(
        "--no-excel",
        action="store_true",
        help="Disable writing Excel file (CSV only).",
    )
    p.add_argument(
        "--no-plots",
        action="store_true",
        help="Disable FRET time-course plot generation.",
    )
    p.add_argument(
        "--plot-tag",
        type=str,
        help="Tag string used in plot and output filenames.",
    )

    return p


def main():
    args = build_argparser().parse_args()

    compute_fret_traces(
        mask_dir=args.mask_dir,
        donor_dir=args.donor_dir,
        acceptor_dir=args.acceptor_dir,
        tracks_csv=args.tracks_csv,
        min_frames=args.min_frames
        if args.min_frames is not None
        else DEFAULTS["min_frames"],
        frame_interval_min=args.frame_interval_min
        if args.frame_interval_min is not None
        else DEFAULTS["frame_interval_min"],
        ring_dilate=args.ring_dilate
        if args.ring_dilate is not None
        else DEFAULTS["ring_dilate"],
        ring_thick=args.ring_thick
        if args.ring_thick is not None
        else DEFAULTS["ring_thick"],
        use_sauvola=(
            False if args.no_sauvola else DEFAULTS["use_sauvola"]
        ),
        sigma_smooth=args.sigma_smooth
        if args.sigma_smooth is not None
        else DEFAULTS["sigma_smooth"],
        sauv_window=args.sauv_window
        if args.sauv_window is not None
        else DEFAULTS["sauv_window"],
        sauv_k=args.sauv_k if args.sauv_k is not None else DEFAULTS["sauv_k"],
        erode_pixels=args.erode_pixels
        if args.erode_pixels is not None
        else DEFAULTS["erode_pixels"],
        alpha_bleed=args.alpha
        if args.alpha is not None
        else DEFAULTS["alpha_bleed"],
        beta_bleed=args.beta
        if args.beta is not None
        else DEFAULTS["beta_bleed"],
        write_excel=(
            False if args.no_excel else DEFAULTS["write_excel"]
        ),
        save_plots=(
            False if args.no_plots else DEFAULTS["save_plots"]
        ),
        plot_tag=args.plot_tag if args.plot_tag is not None else DEFAULTS["plot_tag"],
    )


if __name__ == "__main__":
    main()
