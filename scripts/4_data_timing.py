"""Timing of the Ca2+ response in the image data (wound and uncaging): onset, half-max and peak per ring.

For every responding sector (11.25 deg wedge x 25-um ring, measured half) the (I-min)/min trace is rescaled from its
pre-event level (wound: mean of the 2 frames before the wound frame; uncaging: frame 0) to its peak, and
  t10 = onset (first time at 10 % of the rise), t50 = half-max, t90, tpeak, and the rise fraction at 7/14/21/28 s
are taken (frame-start clock, t = 0 at the wound frame / light-on, linear interpolation, clipped at 0). Each ring is
the median over its responding sectors. t50_minnorm = the arrival-time analysis half-max (half of
the sector's own (I-min)/min peak), median over the responding sectors; for rings >= 25 um it is the per-sector
version of analysis_halfmax, and it is the source (lit circle / edge ring) target of the model fits. ring_um = -1 = inside the uncaging circle (r < 0.8 R; 16 angles x 2 radial
bands, read from the AVG stacks). analysis_halfmax = the arrival time of the wave analysis (ring_summary.csv).

Reads the outputs of notebooks 2 (run once per wound dataset) and 3 (results/wound_wave_analysis/<dataset>,
results/uncaging_wave_analysis) and the
uncaging AVIs (data/uncaging); writes results/data_timing_by_ring.csv. One AVI / TIFF frame = 7 s.
Run after notebooks 2 and 3:  python3 scripts/4_data_timing.py   (paths: see README.md)
"""
import os
import numpy as np, pandas as pd, cv2
from pathlib import Path
REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = Path(os.environ.get("RESULTS_DIR", REPO_ROOT / "results"))
OUT = RESULTS / "data_timing_by_ring.csv"
B = Path(os.environ.get("UNCAGING_DATA_DIR", REPO_ROOT / "data" / "uncaging"))   # the uncaging AVIs
U = RESULTS / "uncaging_wave_analysis"
WD = RESULTS / "wound_wave_analysis"
STACK_FILES = {"ATP1_2": "101626 ATP uncaging (AVG_Stack)/uncage time (full)/AVG_Stack_ATP1_2.avi",
               "ATP3": "101626 ATP uncaging (AVG_Stack)/uncage time (35s or frame 50-400)/AVG_Stack_ATP3.avi",
               "ATP4": "101626 ATP uncaging (AVG_Stack)/uncage time (35s or frame 50-400)/AVG_Stack_ATP4.avi",
               "ATP5": "101626 ATP uncaging (AVG_Stack)/uncage time (35s or frame 50-400)/AVG_Stack_ATP5.avi",
               "Suramin1": "100126 uncaging suramin (AVG_Stack)/AVG_Stack_Suramin1.avi",
               "suramin 2 (wound right)": "100126 uncaging suramin (AVG_Stack)/AVG_Stack_ATP uncaging suramin 2.avi",
               "suramin 3 (wound right)": "100126 uncaging suramin (AVG_Stack)/AVG_Stack_ATP uncaging suramin 3.avi"}


def load_avi(path):
    cap = cv2.VideoCapture(str(path)); fr = []
    while True:
        ok, img = cap.read()
        if not ok:
            break
        fr.append(img[..., 0])
    return np.stack(fr).astype(float)


def half_max_arrival(trace):
    # peak and half-max arrival index along the last axis, interpolated between frames: the arrival-time
    # analysis definition (half of the trace's own (I-min)/min peak), as in observe.half_max_arrival
    filled = np.nan_to_num(trace, nan=-np.inf)
    peak = filled.max(axis=-1)
    thr = 0.5 * peak
    i = np.argmax(filled >= thr[..., None], axis=-1)
    v1 = np.take_along_axis(filled, i[..., None], -1)[..., 0]
    v0 = np.take_along_axis(filled, np.maximum(i - 1, 0)[..., None], -1)[..., 0]
    rising = (i > 0) & (v1 > v0)
    with np.errstate(invalid="ignore"):
        frac = np.divide(thr - v0, v1 - v0, out=np.zeros_like(thr), where=rising)
    return np.where(i > 0, i - 1 + frac, 0.0)


def sector_metrics(t, X, pre):
    # X[sector, time] normalized traces; pre = boolean mask of pre-event frames (none -> first frame)
    b = X[:, pre].mean(1) if pre.any() else X[:, 0]
    post = t >= 0
    Xp = np.where(post[None, :], X, -np.inf)
    k = Xp.argmax(1); P = X[np.arange(len(X)), k]
    F = (X - b[:, None]) / np.maximum(P - b, 1e-9)[:, None]
    def cross(level):
        out = np.full(len(X), np.nan)
        for s in range(len(X)):
            for i in np.nonzero(post)[0]:
                if F[s, i] >= level:
                    if i == 0 or F[s, i - 1] >= level or not post[i - 1] and False:
                        out[s] = t[i]
                    else:
                        out[s] = t[i - 1] + (level - F[s, i - 1]) / (F[s, i] - F[s, i - 1]) * (t[i] - t[i - 1])
                    break
        return np.maximum(out, 0.0)
    fr = lambda ts: F[:, np.isclose(t, ts)][:, 0] if np.isclose(t, ts).any() else np.full(len(X), np.nan)
    return dict(t10=cross(0.1), t50=cross(0.5), t90=cross(0.9), tpeak=t[k], f7=fr(7.0), f14=fr(14.0), f21=fr(21.0),
                f28=fr(28.0))


rows = []
def do(kind, name, group, root, tcol):
    S = pd.read_csv(root / name / "sector_intensity_over_time.csv"); S = S[S.in_measured_range]
    A = pd.read_csv(root / name / "sector_arrival_times.csv"); A = A[A.in_measured_range]
    acol = "arrival_s" if "arrival_s" in A.columns else "arrival_s_after_wound"
    R = pd.read_csv(root / name / "ring_summary.csv")
    arrcol = [c for c in R.columns if c.startswith("median_arrival") and c.endswith("_s") or c == "median_arrival_s_after_wound"][0]
    for off, g in S.groupby("offset_px"):
        ar = A[(A.offset_px == off) & A.responds]
        resp = set(ar.sub_wedge)
        if not resp:
            continue
        piv = g[g.sub_wedge.isin(resp)].pivot_table(index=tcol, columns="sub_wedge", values="normalized_intensity")
        t = piv.index.to_numpy(float); X = piv.to_numpy().T
        m = sector_metrics(t, X, t < 0)
        rr = R[R.offset_px == off]
        nw = g.sub_wedge.nunique()
        rows.append(dict(kind=kind, dataset=name, group=group, ring_um=off, n_resp=len(resp),
                         frac_resp=len(resp) / nw if rr.empty else float(rr.fraction_wedges_responding.iloc[0]),
                         reached=bool(rr.reached.iloc[0]) if not rr.empty else len(resp) / nw >= 0.5,
                         analysis_halfmax=float(rr[arrcol].iloc[0]) if not rr.empty else np.nan,
                         t50_minnorm=float(np.nanmedian(ar[acol])),
                         **{k: float(np.nanmedian(v)) for k, v in m.items()}))


us = pd.read_csv(U / "wave_summary_all_stacks.csv")
for _, s in us.iterrows():
    do("uncaging", s.dataset, s.group, U, "time_s")
    st = load_avi(B / STACK_FILES[s.dataset]); h, w = st.shape[1:]
    yy, xx = np.mgrid[:h, :w]; d = np.hypot(xx - s.circle_cx, yy - s.circle_cy)
    # inside the circle: 16 wedges x 2 radial bands (0-0.4R, 0.4-0.8R) as 'sectors', trace (I-min)/min each
    ang = (np.degrees(np.arctan2(yy - s.circle_cy, xx - s.circle_cx)) % 360 // 22.5).astype(int)
    X = []
    for a in range(16):
        for lo, hi in ((0, 0.4), (0.4, 0.8)):
            msk = (ang == a) & (d >= lo * s.circle_r_px) & (d < hi * s.circle_r_px)
            I = st[:, msk].mean(1); X.append((I - I.min()) / I.min())
    X = np.array(X); t = np.arange(X.shape[1]) * 7.0
    peak = X.max(1); resp = (peak >= 0.2) & (X[:, 0] < 0.5 * peak)
    m = sector_metrics(t, X[resp], np.zeros_like(t, bool)) if resp.any() else None
    rows.append(dict(kind="uncaging", dataset=s.dataset, group=s.group, ring_um=-1, n_resp=int(resp.sum()),
                     frac_resp=float(resp.mean()), reached=bool(resp.mean() >= 0.5), analysis_halfmax=np.nan,
                     t50_minnorm=float(np.median(half_max_arrival(X[resp]) * 7.0)) if resp.any() else np.nan,
                     **({k: float(np.nanmedian(v)) for k, v in m.items()} if m else {})))
# one wave_summary.csv per wound dataset (notebook 2, run once per dataset)
ws = pd.concat([pd.read_csv(f) for f in sorted(WD.glob("*/wave_summary.csv"))], ignore_index=True)
for _, s in ws.iterrows():
    do("wound", s.dataset, s.condition, WD, "time_after_wound_s")
D = pd.DataFrame(rows); D.to_csv(OUT, index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(D[["dataset", "ring_um", "n_resp", "reached", "t10", "t50", "t50_minnorm", "analysis_halfmax", "t90", "tpeak", "f7", "f14", "f21", "f28"]].round(1).to_string())
