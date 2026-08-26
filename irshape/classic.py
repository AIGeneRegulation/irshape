"""CLASSIC (IRFinder-compatible) metric computation -- the only tier Engine
Block 1 implements. See ../SCHEMA.md for exact formulas/thresholds.

  Ai (intron_abundance_Ai) = median per-base coverage over this row's
      UNMASKED (free) intronic positions -- "robust intron depth" per
      CLAUDE.md ("numerator = median coverage over masked intronic positions").
  Ae = max(splice_left, splice_right)
  IRratio_classic = Ai / (Ai + Ae)

Everything here is computed independently for BOTH mask policies (classic,
best_perf) -- one row per (intron_id, mask_policy), per the frozen schema.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import mask as mask_mod
from . import schema as schema_mod

AI_MIN = 5.0
AE_MIN = 5.0
FLANK_MIN = 10.0            # CLAUDE.md mandatory high-confidence depth gate
MINOR_ISOFORM_FRAC = 0.3
NONUNIFORM_RATIO = 3.0


def _flank_values(adapter, flank):
    if flank is None:
        return np.array([], dtype=np.float64)
    chrom, s, e = flank
    return adapter.per_base_values(chrom, s, e)


def _nonuniform_warn(values: np.ndarray, free_bool: np.ndarray, strand: str) -> bool:
    """Coarse 2-half non-uniformity proxy (see SCHEMA.md warn_NonUniformCover
    -- the real multi-bin test is Block 2's SHAPE layer)."""
    L = len(values)
    if L < 4:
        return False
    mid = L // 2
    first_half = values[:mid][free_bool[:mid]]
    second_half = values[mid:][free_bool[mid:]]
    first_half = first_half[~np.isnan(first_half)]
    second_half = second_half[~np.isnan(second_half)]
    if len(first_half) == 0 or len(second_half) == 0:
        return False
    m1, m2 = np.median(first_half), np.median(second_half)
    # strand-orient so "5' half" vs "3' half" is consistent, though the ratio
    # test itself is symmetric and doesn't actually need orientation
    lo, hi = (m1, m2) if m1 <= m2 else (m2, m1)
    if hi <= 0:
        return False
    if lo <= 0:
        return True
    return (hi / lo) > NONUNIFORM_RATIO


def compute_row(intron_ref, other_mask: np.ndarray, same_mask: np.ndarray,
                 small_mask: np.ndarray, adapter, policy: str) -> dict:
    """One (intron_id, mask_policy) row, all Block-1 columns filled."""
    iid, chrom, start, end, strand = (intron_ref.intron_id, intron_ref.chrom,
                                       intron_ref.start, intron_ref.end, intron_ref.strand)
    L = end - start + 1

    values = adapter.per_base_values(chrom, start, end)
    masked_bool = mask_mod.masked_array(other_mask, same_mask, small_mask, policy)
    free_bool = ~masked_bool

    free_vals = values[free_bool]
    free_vals_valid = free_vals[~np.isnan(free_vals)]
    n_free = int(free_bool.sum())

    if n_free == 0:
        Ai = float("nan")
        coverage_fraction = float("nan")
    else:
        Ai = float(np.median(free_vals_valid)) if len(free_vals_valid) else float("nan")
        n_covered = int(np.sum(np.nan_to_num(free_vals, nan=0.0) > 0))
        coverage_fraction = n_covered / n_free

    jc = adapter.junction_counts(iid, chrom, start, end)
    Ae = max(jc.splice_left, jc.splice_right)
    if Ai != Ai:  # NaN
        irratio = float("nan")
    elif (Ai + Ae) > 0:
        irratio = Ai / (Ai + Ae)
    else:
        irratio = float("nan")

    left_vals = _flank_values(adapter, intron_ref.flank_left)
    right_vals = _flank_values(adapter, intron_ref.flank_right)
    flank_all = np.concatenate([left_vals, right_vals]) if (len(left_vals) or len(right_vals)) else np.array([])
    with np.errstate(invalid="ignore"):
        flank_cov = float(np.nanmean(flank_all)) if len(flank_all) else float("nan")

    ai_low = (Ai != Ai) or (Ai < AI_MIN)
    ae_low = Ae < AE_MIN
    flank_low = (flank_cov != flank_cov) or (flank_cov < FLANK_MIN)
    warn_low_cover = bool((ai_low and ae_low) or flank_low)
    warn_low_splicing = bool(ae_low)
    warn_minor_isoform = bool(jc.splice_exact < MINOR_ISOFORM_FRAC * max(jc.splice_left, jc.splice_right, 1.0))
    warn_nonuniform = _nonuniform_warn(values, free_bool, strand)

    mstats = mask_mod.mask_stats(other_mask, same_mask, small_mask, policy)
    n_warn = sum([warn_low_cover, warn_low_splicing, warn_minor_isoform, warn_nonuniform])
    confidence = 1.0 if n_warn == 0 else (0.5 if n_warn == 1 else 0.0)

    return dict(
        intron_id=iid, chrom=chrom, start=start, end=end, strand=strand,
        gene_id=intron_ref.gene_id, host_biotype=intron_ref.host_biotype, length=L,
        length_tier=intron_ref.length_tier,
        intron_abundance_Ai=Ai, splice_left=jc.splice_left, splice_right=jc.splice_right,
        splice_exact=jc.splice_exact, coverage_fraction=coverage_fraction,
        IRratio_classic=irratio,
        warn_LowCover=warn_low_cover, warn_LowSplicing=warn_low_splicing,
        warn_MinorIsoform=warn_minor_isoform, warn_NonUniformCover=warn_nonuniform,
        tag=intron_ref.tag,
        mask_policy=policy, **mstats,
        tier=intron_ref.length_tier, method="classic_ratio", confidence=confidence,
    )


def compute_table(bundle, adapter, intron_ids=None, progress=None) -> pd.DataFrame:
    """Full Block-1 table: 2 rows per intron (one per mask policy), schema-
    conformant (extra SHAPE/GC/IPA/CORRECTED/ANNOTATION columns present as null)."""
    ids = list(bundle.introns) if intron_ids is None else list(intron_ids)
    rows = []
    for k, iid in enumerate(ids):
        ref = bundle.introns[iid]
        om = bundle.other_mask[iid]
        sm = bundle.same_mask[iid]
        small = bundle.small_mask[iid]
        for policy in mask_mod.POLICIES:
            rows.append(compute_row(ref, om, sm, small, adapter, policy))
        if progress is not None and (k + 1) % 500 == 0:
            progress(f"classic metrics: {k + 1}/{len(ids)} introns")

    base = schema_mod.empty_table(len(rows))
    filled = pd.DataFrame(rows)
    for col in filled.columns:
        base[col] = filled[col].values
    for name, dtype in schema_mod.COLUMN_DTYPES.items():
        if name in filled.columns:
            base[name] = base[name].astype(dtype)
    return base
