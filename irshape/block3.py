"""Engine Block 3 orchestration: combine the two-point SHAPE layer
(twopoint.py) and the GC/mappability CONTEXT bias gate (twopoint.py's
window-adapted version of context.py's gate) into the Block-3 columns for
one (intron_id, mask_policy) row, length_tier=="awkward_mid" only
(100-1000bp -- see routing.py / SCHEMA.md).

No IPA call is made from this tier: `cliff_score`/`cliff_position`/
`ipa_cliff_score`/`ipa_position`/`ipa_annotation_support`/`ipa_candidate`
stay null. Two reasons, both already established by 2b/2c: (1) the
degenerate two-point signal has no bin-walk, so there is no cliff position
to even localize an IPA site at; (2) even the FULL shape method's cliff
signal was found to have no usable IPA precision on its own
(`results/ipa_calibration.md`) -- there is no basis to invent a weaker-tier
substitute. IPA stays ranking-only where it exists at all (long tier); this
tier does not attempt it.

`retention_confidence`: derived from how close `two_point_ratio` is to
perfect donor/acceptor BALANCE (0.5 = as retention-like as this signal can
read), then capped at this tier's own validated discriminative ceiling
(`TWO_POINT_AUC_CAP`, from `results/rescue_characterization_summary.md`
Part 3's AUC=0.677) -- see SCHEMA.md's Block 3 addendum for why a cap is
used here but wasn't needed for the long tier (whose uniformity AUC, ~0.94,
already exceeds any such ceiling). Independent of `shape_bias_flag`
everywhere (Block 2c principle carried forward: the bias gate is
diagnostic, never coupled into the confidence score).
"""
from __future__ import annotations

import numpy as np

from . import context as context_mod
from . import mask as mask_mod
from . import twopoint as twopoint_mod
from .block2 import apply_irratio_corrected  # noqa: F401  (re-exported for run-script convenience)

# results/rescue_characterization_summary.md Part 3: AUC(donor_share, IPA vs
# RETENTION) = 0.677 on the awkward-middle panel. Used as an explicit
# confidence CEILING: even a perfectly-balanced (ratio==0.5) two-point read
# is never asserted with more confidence than this tier's own measured
# discriminative power.
TWO_POINT_AUC_CAP = 0.68


def compute_twopoint_row(intron_ref, other_mask: np.ndarray, same_mask: np.ndarray,
                          small_mask: np.ndarray, gc_arr, map_arr, adapter,
                          mask_policy: str) -> dict:
    """All Block-3 columns for one (intron, mask_policy) row EXCEPT
    IRratio_corrected (needs IRratio_classic from the Block-1 classic table,
    merged in afterward by the caller -- see block2.apply_irratio_corrected,
    reused unchanged here since its formula is tier-agnostic)."""
    iid, chrom, istart, iend, strand = (intron_ref.intron_id, intron_ref.chrom,
                                         intron_ref.start, intron_ref.end, intron_ref.strand)
    values = adapter.per_base_values(chrom, istart, iend)
    masked_bool = mask_mod.masked_array(other_mask, same_mask, small_mask, mask_policy)
    free_bool = ~masked_bool

    tp = twopoint_mod.compute_two_point(values, free_bool, istart, iend, strand)

    gc_content = context_mod.masked_track_mean(gc_arr, free_bool) if gc_arr is not None else float("nan")
    mappability_score = (context_mod.masked_track_mean(map_arr, free_bool)
                          if map_arr is not None else float("nan"))

    bias_flag, gc_flag = False, False
    if tp["shape_fittable"]:
        bias_flag, gc_flag = twopoint_mod.two_point_bias_flag(
            gc_arr, map_arr, free_bool, istart, tp["donor_range"], tp["acceptor_range"]
        )

    ratio = tp["two_point_ratio"]
    retention_confidence = float("nan")
    if tp["shape_fittable"] and ratio == ratio:
        balance = float(np.clip(1.0 - 2.0 * abs(ratio - 0.5), 0.0, 1.0))
        retention_confidence = TWO_POINT_AUC_CAP * balance

    return dict(
        intron_id=iid, mask_policy=mask_policy,
        two_point_ratio=ratio, shape_fittable=tp["shape_fittable"], shape_method="two_point",
        gc_content=gc_content, mappability_score=mappability_score,
        gc_transition_near_shape_feature=gc_flag, shape_bias_flag=bias_flag,
        retention_confidence=retention_confidence,
    )
