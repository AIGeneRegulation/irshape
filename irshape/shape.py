"""Engine Block 2 Part 2 -- full coverage-SHAPE metrics (long-tier only,
shape_method="full_shape"): uniformity (entropy-based flatness) and
cliff_score/cliff_position (largest internal coverage step-down, the IPA/
alt-3' signature).

uniformity: Shannon entropy of the per-bin coverage distribution, normalized
by max entropy (Pielou evenness) -- the "robust entropy formulation" the
Block-2 task calls for, in place of Block-1-era scripts' `1 - CV`:
  - bounded in [0,1] without ad hoc clipping (CV is not, and needed
    max(0, ...) clamping in the old code)
  - degrades gracefully with a few zero-coverage bins (a bin with p_i=0
    contributes 0 to entropy) instead of blowing up when mean coverage is
    small, which is exactly the failure mode that produced many `Ai=0` rows
    for long, sparsely-covered introns in Block 1's classic output
  - well precedented for this exact "read-distribution uniformity" role in
    intron-retention tools (the iREAD-style evenness test this task's spec
    references), rather than an assumption of Gaussian-like spread.

cliff_score/cliff_position: same largest-single-adjacent-valid-bin-drop
definition validated in results/shape_evidence_final.md (AUC ~0.72-0.77 for
IPA-vs-background) -- unchanged from the historical scripts, since only
uniformity's formulation was asked to be revisited.
"""
from __future__ import annotations

import numpy as np

from .binning import MIN_BIN_FRAC, MIN_RAW_LEN, MIN_VALID_BINS, adaptive_nbin, bin_ranges


def bin_coverage(values: np.ndarray, free_bool: np.ndarray, istart: int, iend: int):
    """Bin per-base coverage (only over UNMASKED/free positions) into
    adaptive genomic bins. Returns (nbin, valid[bool per bin], cov[float per
    bin], mid[genomic midpoint per bin]) in genomic (not strand-oriented)
    order. A bin is "valid" if >= MIN_BIN_FRAC of its raw span survives
    masking and has real (non-NaN) coverage data."""
    length = iend - istart + 1
    if length < MIN_RAW_LEN:
        return 0, [], [], []
    nbin = adaptive_nbin(length)
    ranges = bin_ranges(istart, iend, nbin)
    valid, cov, mid = [], [], []
    for (bs, be) in ranges:
        lo, hi = bs - istart, be - istart + 1
        seg_vals = values[lo:hi]
        seg_free = free_bool[lo:hi]
        raw_len = hi - lo
        kept = int(seg_free.sum())
        frac = (kept / raw_len) if raw_len else 0.0
        seg_valid_vals = seg_vals[seg_free]
        seg_valid_vals = seg_valid_vals[~np.isnan(seg_valid_vals)]
        is_valid = (frac >= MIN_BIN_FRAC) and (len(seg_valid_vals) > 0)
        c = float(np.mean(seg_valid_vals)) if len(seg_valid_vals) else float("nan")
        valid.append(bool(is_valid))
        cov.append(c)
        mid.append((bs + be) // 2)
    return nbin, valid, cov, mid


def entropy_uniformity(cov_valid) -> float:
    """Pielou evenness (Shannon entropy / max entropy) of the valid-bin
    coverage distribution. 1.0 = perfectly flat (all bins equal, retention-
    like); 0.0 = all coverage concentrated in one bin (cliff-like). NaN if
    fewer than 2 valid bins or all-zero coverage (undefined)."""
    arr = np.clip(np.asarray(cov_valid, dtype=np.float64), 0, None)
    n = len(arr)
    total = arr.sum()
    if n < 2 or total <= 0:
        return float("nan")
    p = arr / total
    nz = p[p > 0]
    h = -np.sum(nz * np.log(nz))
    h_max = np.log(n)
    return float(np.clip(h / h_max, 0.0, 1.0))


def cliff(cov_valid, mid_valid, strand: str):
    """Largest single-step relative coverage drop between ADJACENT valid
    bins, walked 5'->3' (strand-oriented). Returns (cliff_score,
    cliff_position genomic 1-based, or None if <2 valid bins)."""
    cov = list(cov_valid)
    mid = list(mid_valid)
    if strand == "-":
        cov = cov[::-1]
        mid = mid[::-1]
    best_score, best_pos = 0.0, None
    for k in range(len(cov) - 1):
        before, after = cov[k], cov[k + 1]
        drop = (before - after) / max(before, 1.0)
        if drop > best_score:
            best_score = drop
            best_pos = (mid[k] + mid[k + 1]) // 2
    return float(best_score), best_pos


def compute_shape(values: np.ndarray, free_bool: np.ndarray, istart: int, iend: int, strand: str) -> dict:
    """Full per-(intron, mask_policy) shape result. `bin_valid`/`bin_mid` are
    returned too so the CONTEXT layer (GC/mappability) can bin on the exact
    same genomic ranges for the Part-3 bias gate."""
    nbin, valid, cov, mid = bin_coverage(values, free_bool, istart, iend)
    n_valid = sum(valid)
    fittable = (iend - istart + 1) >= MIN_RAW_LEN and n_valid >= MIN_VALID_BINS
    if not fittable:
        return dict(n_bins=nbin, n_valid_bins=n_valid, uniformity=float("nan"),
                    cliff_score=float("nan"), cliff_position=None, shape_fittable=False,
                    bin_valid=valid, bin_mid=mid)
    cov_valid = [c for v, c in zip(valid, cov) if v]
    mid_valid = [m for v, m in zip(valid, mid) if v]
    uniformity = entropy_uniformity(cov_valid)
    cliff_score, cliff_position = cliff(cov_valid, mid_valid, strand)
    return dict(n_bins=nbin, n_valid_bins=n_valid, uniformity=uniformity,
                cliff_score=cliff_score, cliff_position=cliff_position, shape_fittable=True,
                bin_valid=valid, bin_mid=mid)
