"""Sheet-wave pixel selection, shared by the wound (notebook 2) and uncaging (notebook 3) analyses.

The Ca2+ movies show two kinds of signal:
  * the WAVE: one connected sheet of responding cells that grows out from the source (the wound edge or the
    uncaging circle), and
  * single cells that FLASH on their own, before, during and long after the wave, far from the front; they are
    also present without a wave (e.g. uncaging under suramin) and before the wound is made.
A sector (wedge x ring) average mixes the two, so a flash can set a sector's half-max "arrival". This module
finds, pixel by pixel, which pixels belong to the sheet and when the sheet reached them; the notebooks then
measure every sector on its sheet pixels only, with the same (I - min)/min half-max rule as before.

Steps (all on the frames the notebook analyses; t = frame index in that range):
  1. smooth_norm  : light Gaussian smoothing of every frame (SMOOTH_SIGMA_PX), then (I - min)/min per pixel
                    (min over the analysed frames) and each pixel's frame-to-frame noise (robust SD).
  2. tissue_mask  : pixels whose minimum over time is clearly above the dark background (outside the fin or the
                    illuminated field); ring sectors count only tissue pixels.
  3. active_map   : a pixel is ACTIVE in frame t when its (I - min)/min has RISEN by at least ACTIVE_RISE (and
                    ACTIVE_NOISE_K x its noise) above its own lowest value in the last RISE_WINDOW frames; a cell
                    that is bright and slowly dimming (e.g. bleaching) is not active. Specks smaller than a cell
                    are removed (binary opening).
  4. grow_sheet   : the sheet starts at the source (seed) and grows outward: in every frame, each active region
                    that touches the FRONT joins it. The front is the source (while it is on) and the pixels the
                    sheet reached in the last RECRUIT_WINDOW frames, so pixels reached long ago do not recruit:
                    cells flashing next to the sheet minutes after it passed cannot chain a false wave outward.
                    t_reach = the frame a pixel joined.
  5. classify     : a pixel that joined the sheet is TIMED unless it had already been active, without a break,
                    for more than LATE_TOL frames before it joined (its own earlier flash was still going on, so
                    its rise cannot be timed: EARLY). Active pixels the sheet never reached are FLASH_ONLY.
The notebooks then average every sector over its TIMED pixels only.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import gaussian

# method settings used in the paper (pixels / frames; 1 px = 1 um, 1 frame = 7 s)
SHEET_SETTINGS = dict(
    SMOOTH_SIGMA_PX=2.0,     # Gaussian smoothing of every frame before the pixel maps
    TISSUE_K=3.0,            # tissue = minimum over time > dark level + TISSUE_K x SD of the darkest 10 % of pixels
    ACTIVE_RISE=0.3,         # active = rise of (I - min)/min of at least this ...
    ACTIVE_NOISE_K=4.0,      # ... and of at least this many times the pixel's own frame-to-frame noise ...
    RISE_WINDOW=3,           # ... above its lowest value in the last RISE_WINDOW frames (3 frames = 21 s)
    OPEN_ITER=1,             # binary opening of the active map (removes specks of a few pixels)
    LATE_TOL=1,              # a pixel active for <= LATE_TOL frames before the sheet reached it is still timed
    RECRUIT_WINDOW=3,        # the front = source + pixels reached in the last RECRUIT_WINDOW frames (21 s)
    ONSET_ONLY=True,         # a pixel can join only in a frame in which it STARTS to be active (not while it stays active)
)


def smooth_norm(stack, sigma):
    """Smoothed frames, per-pixel minimum over time, (I - min)/min and the per-pixel frame-to-frame noise."""
    sm = np.stack([gaussian(f, sigma, preserve_range=True) for f in stack]) if sigma > 0 else stack.astype(float)
    fmin = sm.min(0)
    norm = (sm - fmin) / np.maximum(fmin, 1e-12)
    d = np.diff(norm, axis=0)
    noise = 1.4826 * np.median(np.abs(d - np.median(d, 0)), 0) / np.sqrt(2)
    return sm, fmin, norm, noise


def tissue_mask(fmin, k):
    """Pixels clearly brighter (minimum over time) than the dark background."""
    lo = np.percentile(fmin, 2)
    dark = fmin[fmin <= np.percentile(fmin, 10)]
    thr = lo + k * max(dark.std(), 0.01 * (np.percentile(fmin, 99) - lo))
    return ndi.binary_opening(fmin > thr, iterations=2), float(thr)


def active_map(norm, noise, mask, rise, noise_k, window, open_iter):
    """active[t]: the pixel's (I - min)/min rose by >= max(rise, noise_k x noise) above its minimum of the last
    `window` frames (frames t-window .. t)."""
    T = norm.shape[0]
    runmin = np.stack([norm[max(0, t - window):t + 1].min(0) for t in range(T)])
    act = ((norm - runmin) >= np.maximum(rise, noise_k * noise)) & mask
    if open_iter:
        act = np.stack([ndi.binary_opening(a, iterations=open_iter) for a in act])
    return act


def grow_sheet(act, mask, seed, t_start=0, recruit_window=None, onset_only=False):
    """The sheet grows from the seed. seed(t) -> bool map of the source in frame t, or None (source off).
    Each active region touching the front joins; the front = the source in frame t + the pixels reached in frames
    t - recruit_window .. t - 1 (all reached pixels if recruit_window is None).
    Returns t_reach (frame the pixel joined the sheet; -1 = never)."""
    T, H, W = act.shape
    reached = np.zeros((H, W), bool)
    t_reach = np.full((H, W), -1, int)
    for t in range(t_start, T):
        s = seed(t)
        src = np.zeros((H, W), bool) if s is None else (s & mask)
        new = src & ~reached
        reached |= new
        t_reach[new] = t
        front = src | (reached if recruit_window is None else (reached & (t_reach >= t - recruit_window)))
        if not front.any():
            continue
        cand = act[t] & ~act[t - 1] if (onset_only and t > 0) else act[t]
        lab, _ = ndi.label(cand | front)
        keep = np.unique(lab[front])
        grow = np.isin(lab, keep[keep > 0]) & ~reached & cand
        reached |= grow
        t_reach[grow] = t
    return t_reach


def classify(act, t_reach, late_tol, t_start=0):
    """TIMED / EARLY / FLASH_ONLY maps (see module docstring)."""
    T = act.shape[0]
    timed = np.zeros(t_reach.shape, bool)
    early = np.zeros(t_reach.shape, bool)
    for t in range(t_start, T):
        sel = t_reach == t
        if not sel.any():
            continue
        lo = t - late_tol - 1
        ongoing = np.all(act[lo:t], axis=0) if lo >= 0 else np.zeros(sel.shape, bool)
        early |= sel & ongoing
        timed |= sel & ~ongoing
    flash_only = act[t_start:].any(0) & (t_reach < 0)
    return timed, early, flash_only


def sheet_pixels(stack, seed, t_start=0, exclude=None, settings=None):
    """Run steps 1-5 on `stack` [frame, y, x]. `seed(t)` gives the source in frame t; `exclude` (bool map) removes
    pixels from the tissue (e.g. the wound hole). Returns a dict of the maps."""
    s = {**SHEET_SETTINGS, **(settings or {})}
    _, fmin, norm, noise = smooth_norm(stack, s["SMOOTH_SIGMA_PX"])
    mask, thr = tissue_mask(fmin, s["TISSUE_K"])
    if exclude is not None:
        mask &= ~exclude
    act = active_map(norm, noise, mask, s["ACTIVE_RISE"], s["ACTIVE_NOISE_K"], s["RISE_WINDOW"], s["OPEN_ITER"])
    t_reach = grow_sheet(act, mask, seed, t_start, s["RECRUIT_WINDOW"], s["ONSET_ONLY"])
    timed, early, flash_only = classify(act, t_reach, s["LATE_TOL"], t_start)
    return dict(tissue=mask, tissue_threshold=thr, active=act, t_reach=t_reach, timed=timed, early=early,
                flash_only=flash_only, settings=s)
