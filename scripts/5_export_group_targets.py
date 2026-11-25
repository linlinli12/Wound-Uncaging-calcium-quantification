"""Group-average fit targets from the arrival-time analysis ((I-min)/min), one row per group x ring.

Reads the outputs of notebooks 2 and 3 (results/wound_wave_analysis/<dataset>/wave_summary.csv, one per dataset,
results/uncaging_wave_analysis/wave_summary_all_stacks.csv and each dataset's ring_summary.csv) and, with --source,
results/data_timing_by_ring.csv from scripts/4_data_timing.py. Writes results/<name> (default group_targets.csv).

    python3 scripts/5_export_group_targets.py group_targets.csv --source --onset

Per group and ring: mean / SD / n of the per-dataset median arrival times, using only rings inside each
dataset's response reach; mean fraction of wedges responding. Per group: number of datasets and waves,
mean / SD response reach, mean / SD speed (outward waves only), representative geometry (mean circle /
wound radius), light duration for uncaging.
--source adds one row per group for the source: the lit circle (uncaging, offset -1) or the wound-edge ring
(wound, offset 0), from the half-max arrival of the responding sectors (t50_minnorm).
--onset (with --source) adds mean_onset_s / sd_onset_s / n_datasets_onset: the onset (first time at 10 % of the
rise, t10) at every ring, over the datasets that reach it.
Times are in seconds (7 s per frame); t = 0 at light-on (uncaging) or the wound frame (wound).
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = Path(os.environ.get("RESULTS_DIR", REPO_ROOT / "results"))
UNCAGE = RESULTS / "uncaging_wave_analysis"
WOUND = RESULTS / "wound_wave_analysis"
# output file name: first argument, default group_targets.csv
_ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
OUT = RESULTS / (_ARGS[0] if _ARGS else "group_targets.csv")
SOURCE = "--source" in sys.argv[1:]
ONSET = "--onset" in sys.argv[1:]
assert SOURCE or not ONSET, "--onset needs --source"
TIMING = RESULTS / "data_timing_by_ring.csv"

USE = {"ATP, uncage full": "fit_2B", "ATP, uncage 35 s": "fit_2B", "suramin, uncage full": "validate",
       "ctrl": "fit_2A", "heptanol": "fit_stage1", "POM1": "validate", "suramin": "validate",
       "CBX": "validate", "Gd3+": "validate"}
LIGHT_S = {"ATP, uncage full": 280.0, "ATP, uncage 35 s": 35.0, "suramin, uncage full": 280.0}


def group_rows(summary_csv, root, key, arr_col, radius_col, experiment):
    if (root / summary_csv).exists():
        s = pd.read_csv(root / summary_csv)
    else:   # wound: one wave_summary.csv per dataset (notebook 2, run once per dataset)
        s = pd.concat([pd.read_csv(f) for f in sorted(root.glob("*/wave_summary.csv"))], ignore_index=True)
    rows = []
    for g, sub in s.groupby(key, sort=False):
        arr, frac = [], []
        for _, d in sub.iterrows():
            ring = pd.read_csv(root / d.dataset / "ring_summary.csv")
            a = ring[arr_col].to_numpy(float)
            a = np.where(ring.reached.to_numpy(bool), a, np.nan)     # only rings inside this dataset's reach
            arr.append(a)
            frac.append(ring.fraction_wedges_responding.to_numpy(float))
        arr, frac = np.array(arr), np.array(frac)
        n = np.isfinite(arr).sum(0)
        mean = np.where(n > 0, np.nanmean(arr, axis=0), np.nan)
        sd = np.where(n > 1, np.nanstd(arr, axis=0, ddof=1), np.nan)
        waves = sub[sub.outward_wave]
        base = {"experiment": experiment, "group": g, "use": USE.get(g, "validate"),
                "datasets": ", ".join(sub.dataset), "n_datasets": len(sub), "n_outward_waves": int(len(waves)),
                "mean_reach_um": sub.response_reach_um.mean(),
                "sd_reach_um": sub.response_reach_um.std(ddof=1) if len(sub) > 1 else np.nan,
                "mean_speed_um_per_s": waves.speed_um_per_s.mean() if len(waves) else np.nan,
                "sd_speed_um_per_s": waves.speed_um_per_s.std(ddof=1) if len(waves) > 1 else np.nan,
                "mean_radius_um": sub[radius_col].mean(),
                "light_duration_s": LIGHT_S.get(g, np.nan),
                "time_zero": "light on (movie frame 0)" if experiment == "uncaging" else "wound frame"}
        for k, off in enumerate(range(25, 225, 25)):
            rows.append({**base, "offset_um": off, "mean_arrival_s": mean[k], "sd_arrival_s": sd[k],
                         "n_datasets_reaching_ring": int(n[k]), "mean_fraction_wedges_responding": frac[:, k].mean()})
    return rows


warnings.filterwarnings("ignore", category=RuntimeWarning)   # nanmean / nanstd of rings no dataset reaches
rows = group_rows("wave_summary_all_stacks.csv", UNCAGE, "group", "median_arrival_s", "circle_r_px", "uncaging")
rows += group_rows("wave_summary_all_datasets.csv", WOUND, "condition", "median_arrival_s_after_wound",
                   "hole_radius_at_wound_px", "wound")
df = pd.DataFrame(rows)
if SOURCE:
    T = pd.read_csv(TIMING)
    keep = ["experiment", "group", "use", "datasets", "n_datasets", "n_outward_waves", "mean_reach_um", "sd_reach_um",
            "mean_speed_um_per_s", "sd_speed_um_per_s", "mean_radius_um", "light_duration_s", "time_zero"]
    out = []
    for g, sub in df.groupby("group", sort=False):
        base = {k: sub.iloc[0][k] for k in keep}
        ring = -1 if base["experiment"] == "uncaging" else 0
        t = T[(T.group == g) & (T.ring_um == ring)]
        v = t[t.reached].t50_minnorm.to_numpy(float)
        out.append(pd.DataFrame([{**base, "offset_um": ring, "location": "lit circle" if ring < 0 else "edge ring",
                                  "mean_arrival_s": v.mean() if len(v) else np.nan,
                                  "sd_arrival_s": v.std(ddof=1) if len(v) > 1 else np.nan,
                                  "n_datasets_reaching_ring": len(v),
                                  "mean_fraction_wedges_responding": t.frac_resp.mean()}]))
        out.append(sub.assign(location="ring"))
    df = pd.concat(out, ignore_index=True)
    if ONSET:
        on = []
        for _, r in df.iterrows():
            t = T[(T.group == r.group) & (T.ring_um == r.offset_um) & T.reached].t10.to_numpy(float)
            t = t[np.isfinite(t)]
            on.append((t.mean() if len(t) else np.nan, t.std(ddof=1) if len(t) > 1 else np.nan, len(t)))
        df["mean_onset_s"], df["sd_onset_s"], df["n_datasets_onset"] = zip(*on)
df.to_csv(OUT, index=False)
print(f"wrote {OUT} ({len(df)} rows, {df.group.nunique()} groups)")
