"""Engine Block 3 -- degenerate two-point donor/acceptor coverage signal for
length_tier=="awkward_mid" (100-1000bp) introns, where the multi-bin SHAPE
fitter (shape.py) can't be relied on -- most of its unfittability in this
range is a MASKING problem, not a raw-length one (73% of pooled-RETENTION
shape failures in `results/rescue_characterization_summary.md` Part 1 are
MASK_HEAVY), so a signal that only needs ~15bp unmasked at EACH end (not
several full bins) recovers most of the range.

`two_point_ratio` = donor-window / (donor+acceptor-window) mean coverage
share (SCHEMA.md SHAPE LAYER): true retention has both ends covered
comparably (ratio ~0.5); IPA/internal-end-use has the donor end covered and
the acceptor end dropped (ratio > 0.5, "donor-heavy"). As stated in
`rescue_characterization_summary.md` Part 3, this necessarily also matches
plain nascent 5'->3' read-through (same coverage pattern) -- it is a
retention-vs-NOT-retention rescue at these lengths, not an IPA-vs-nascent
discriminator. No IPA call is derived from this signal anywhere in the
engine (see block3.py).

WINDOW/MIN_UNMASKED are the exact parameters validated in that
characterization (AUC(donor_share, IPA vs RETENTION) = 0.677, computable for
432/458 = 94.3% of the awkward-middle RETENTION+IPA panel) -- unchanged here
so the engine reproduces that number rather than a re-tuned one.
"""
from __future__ import annotations

import numpy as np

from .context import GC_TRANSITION_MIN_DELTA, MAPPABILITY_LOW_THRESHOLD

WINDOW = 30
MIN_UNMASKED = 15


def donor_acceptor_windows(istart: int, iend: int, strand: str):
    """(donor_range, acceptor_range), genomic 1-based inclusive, strand-aware:
    donor = just inside the transcript-5' (donor splice site) boundary,
    acceptor = just inside the transcript-3' (acceptor splice site)
    boundary. Matches scripts/43_rescue_twopoint.py's window definition
    exactly."""
    if strand == "+":
        d = (istart, min(istart + WINDOW - 1, iend))
        a = (max(iend - WINDOW + 1, istart), iend)
    else:
        d = (max(iend - WINDOW + 1, istart), iend)
        a = (istart, min(istart + WINDOW - 1, iend))
    return d, a


def _window_slice(istart: int, window: tuple):
    ws, we = window
    return ws - istart, we - istart + 1


def window_mean(values: np.ndarray, free_bool: np.ndarray, istart: int, window: tuple):
    """(kept_bp, mean_cov) over the UNMASKED bases of one window. kept_bp is
    the actually-unmasked base count (vs. the window's raw span), the
    quantity `MIN_UNMASKED` gates on."""
    lo, hi = _window_slice(istart, window)
    seg_vals = values[lo:hi]
    seg_free = free_bool[lo:hi]
    kept_bp = int(seg_free.sum())
    seg_valid = seg_vals[seg_free]
    seg_valid = seg_valid[~np.isnan(seg_valid)]
    mean = float(np.mean(seg_valid)) if len(seg_valid) else float("nan")
    return kept_bp, mean


def compute_two_point(values: np.ndarray, free_bool: np.ndarray, istart: int, iend: int,
                       strand: str) -> dict:
    """Per-(intron, mask_policy) two-point result. Mirrors shape.compute_shape's
    return-dict style so block3.py can slot it into the same per-row assembly
    pattern block2.py uses for the full-shape tier."""
    d_range, a_range = donor_acceptor_windows(istart, iend, strand)
    donor_kept, donor_cov = window_mean(values, free_bool, istart, d_range)
    acceptor_kept, acceptor_cov = window_mean(values, free_bool, istart, a_range)
    fittable = (donor_kept >= MIN_UNMASKED and acceptor_kept >= MIN_UNMASKED
                and donor_cov == donor_cov and acceptor_cov == acceptor_cov
                and (donor_cov + acceptor_cov) > 0)
    ratio = (donor_cov / (donor_cov + acceptor_cov)) if fittable else float("nan")
    return dict(donor_range=d_range, acceptor_range=a_range,
                donor_kept_bp=donor_kept, acceptor_kept_bp=acceptor_kept,
                donor_cov=donor_cov, acceptor_cov=acceptor_cov,
                two_point_ratio=ratio, shape_fittable=bool(fittable))


def two_point_bias_flag(gc_arr, map_arr, free_bool: np.ndarray, istart: int,
                         d_range: tuple, a_range: tuple) -> tuple:
    """Adapts context.py's GC/mappability bias GATE (SCHEMA.md GC_MAPPABILITY)
    to the two-point method. The full-shape gate asks "does a GC/mappability
    transition coincide with the cliff bin" -- there is no bin-walk or cliff
    position here, so the natural equivalent is "does a GC/mappability
    difference exist BETWEEN the donor and acceptor windows" -- the exact two
    quantities `two_point_ratio` itself compares. Same thresholds as
    context.py (GC_TRANSITION_MIN_DELTA=0.15, MAPPABILITY_LOW_THRESHOLD=0.5),
    so "the same bias flag" carries over, just re-pointed at this method's
    two windows instead of adjacent coverage bins.

    Returns (shape_bias_flag, gc_transition_near_shape_feature), matching
    context.shape_bias_flag's return signature. Deliberately NOT consulted by
    retention_confidence anywhere (Block 2c principle carried forward: this
    flag is diagnostic, never coupled into the confidence score)."""

    def wmean(arr, window):
        if arr is None:
            return float("nan")
        lo, hi = _window_slice(istart, window)
        seg = arr[lo:hi][free_bool[lo:hi]]
        seg = seg[~np.isnan(seg)]
        return float(np.mean(seg)) if len(seg) else float("nan")

    gc_d, gc_a = wmean(gc_arr, d_range), wmean(gc_arr, a_range)
    map_d, map_a = wmean(map_arr, d_range), wmean(map_arr, a_range)

    gc_flag = bool(gc_d == gc_d and gc_a == gc_a and abs(gc_d - gc_a) >= GC_TRANSITION_MIN_DELTA)
    map_flag = bool((map_d == map_d and map_d < MAPPABILITY_LOW_THRESHOLD)
                     or (map_a == map_a and map_a < MAPPABILITY_LOW_THRESHOLD))
    return bool(gc_flag or map_flag), gc_flag
