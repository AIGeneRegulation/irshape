"""Engine Block 2 Part 1/3 -- GC content + mappability as first-class,
per-base reference tracks (see reference.py for how they're built/cached),
and the GC/mappability GATING that flags a shape feature as plausibly
technical rather than biological.

Per SCHEMA.md, the frozen contract has exactly one named boolean for this,
`gc_transition_near_shape_feature` -- kept strictly about GC. A mappability
drop at the cliff also feeds `shape_bias_flag` (the actual gate asked for in
Part 3, "GC transition OR mappability drop"), just without its own dedicated
column, since Block 1 didn't declare one; this is documented here rather than
silently expanding an already-named column's semantics.
"""
from __future__ import annotations

import numpy as np

from .binning import bin_ranges

TRANSITION_WINDOW = 500          # bp; matches this project's established
                                  # cliff-vs-annotation position tolerance
                                  # (results/shape_characterization_summary.md
                                  # Part 2, rescue_characterization Part 4)
GC_TRANSITION_MIN_DELTA = 0.15   # min |delta GC fraction| between adjacent
                                  # bins to call it a "transition" at all
MAPPABILITY_LOW_THRESHOLD = 0.5  # mean mappability below this, at the bin
                                  # nearest the cliff, is a low-mappability flag


def bin_track_full(values: np.ndarray, istart: int, iend: int, nbin: int):
    """Mean of `values` per genomic bin, over the FULL bin span (NOT masked-
    filtered -- GC/mappability are sequence properties that exist regardless
    of whether a position happens to be masked for coverage purposes)."""
    ranges = bin_ranges(istart, iend, nbin)
    means = []
    for (bs, be) in ranges:
        lo, hi = bs - istart, be - istart + 1
        seg = values[lo:hi]
        seg = seg[~np.isnan(seg)]
        means.append(float(np.mean(seg)) if len(seg) else float("nan"))
    return means, ranges


def masked_track_mean(values: np.ndarray, free_bool: np.ndarray) -> float:
    """Per-row (mask-policy-dependent) summary: mean over UNMASKED bases --
    used for `gc_content` and `mappability_score`."""
    v = values[free_bool]
    v = v[~np.isnan(v)]
    return float(np.mean(v)) if len(v) else float("nan")


def _largest_adjacent_transition(bin_means, ranges):
    """Largest |delta| between ADJACENT bins with real (non-NaN) means, and
    the genomic boundary position between them. (0.0, None) if <2 such bins."""
    best_delta, best_pos = 0.0, None
    prev_val, prev_end = None, None
    for (bs, be), val in zip(ranges, bin_means):
        if val == val:  # not NaN
            if prev_val is not None:
                d = abs(val - prev_val)
                if d > best_delta:
                    best_delta = d
                    best_pos = (prev_end + bs) // 2
            prev_val, prev_end = val, be
    return best_delta, best_pos


def _nearest_bin_value(bin_means, ranges, position):
    if position is None:
        return None, None
    best_i, best_d = None, None
    for i, (bs, be) in enumerate(ranges):
        mid = (bs + be) // 2
        d = abs(mid - position)
        if best_d is None or d < best_d:
            best_d, best_i = d, i
    if best_i is None or best_d > TRANSITION_WINDOW:
        return None, None
    return bin_means[best_i], best_d


def gc_transition_near_cliff(gc_values: np.ndarray, istart: int, iend: int, nbin: int,
                              cliff_position) -> bool:
    """SCHEMA.md `gc_transition_near_shape_feature`: does the sharpest
    adjacent-bin GC-content jump sit within TRANSITION_WINDOW bp of
    cliff_position, and is it large enough (>= GC_TRANSITION_MIN_DELTA) to
    call a "transition" at all (not just natural base-composition noise)?"""
    if cliff_position is None or nbin == 0:
        return False
    bin_means, ranges = bin_track_full(gc_values, istart, iend, nbin)
    delta, pos = _largest_adjacent_transition(bin_means, ranges)
    if pos is None or delta < GC_TRANSITION_MIN_DELTA:
        return False
    return abs(pos - cliff_position) <= TRANSITION_WINDOW


def mappability_drop_near_cliff(map_values: np.ndarray, istart: int, iend: int, nbin: int,
                                 cliff_position) -> bool:
    """A broadly low-mappability bin AT the cliff (not necessarily a sharp
    step -- a whole low-mappability stretch overlapping the cliff is itself
    the bias scenario, e.g. a repeat-masked region producing spurious
    coverage dropout that looks like an IPA cliff but isn't one)."""
    if cliff_position is None or nbin == 0:
        return False
    bin_means, ranges = bin_track_full(map_values, istart, iend, nbin)
    val, dist = _nearest_bin_value(bin_means, ranges, cliff_position)
    if val is None or val != val:
        return False
    return val < MAPPABILITY_LOW_THRESHOLD


def shape_bias_flag(gc_values, map_values, istart, iend, nbin, cliff_position) -> tuple:
    """Returns (shape_bias_flag, gc_transition_near_shape_feature). The
    overall gate ORs the GC-transition test with the mappability-drop test
    (Part 3's "GC transition OR a mappability drop"); only
    gc_transition_near_shape_feature is a named schema column."""
    gc_flag = gc_transition_near_cliff(gc_values, istart, iend, nbin, cliff_position) if gc_values is not None else False
    map_flag = mappability_drop_near_cliff(map_values, istart, iend, nbin, cliff_position) if map_values is not None else False
    return bool(gc_flag or map_flag), bool(gc_flag)
