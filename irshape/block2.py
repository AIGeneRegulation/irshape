"""Engine Block 2 orchestration: combine the SHAPE layer (shape.py), the GC/
mappability CONTEXT + bias gate (context.py), the IPA candidate/annotation
signal (annotation.py), and the corrected-output layer, into the Block-2
columns for one (intron_id, mask_policy) row.

Scope: length_tier=="long" only (>=1000bp) -- see reference.py's `min_length`
param, which the Block-2 scripts set to 1000 when building the reference
bundle these functions consume.

Engine Block 2b (correctness fix, see results/ipa_calibration.md and
results/engine_block2b_validation.md): the IPA cliff signal has NO usable
precision as a classifier -- the best rule found (cliff>=0.80 AND
polyA-near-cliff) only reaches ~9.3% precision, and the prior binary
`ipa_flag` (cliff>=0.80 alone) was 96.7% false positives among confident-
long-read-adjudicable introns. Consequences of that finding, implemented
here:
  - `IRratio_corrected` no longer applies any IPA-driven down-weight -- it
    is simply the best_perf-mask row's classic-formula retention estimate
    (median Ai / (Ai+Ae) under the validated best_perf mask), i.e. equal to
    that row's own `IRratio_classic`. See `apply_irratio_corrected`.
  - `ipa_flag` and `ipa_confidence` are RETIRED as determinations: this
    module no longer computes real values for them (Block 2b's engine
    scripts leave them null in the output table, per the frozen schema's
    "declared, not computed" convention) -- they stay in the schema for
    backward compatibility, not as fields callers should read.
  - `ipa_candidate` replaces them as an explicitly HYPOTHESIS-GENERATING
    signal: cliff>=0.80 AND NOT shape_bias_flag AND polyA-near-cliff (the
    single best rule found in calibration), with its measured
    precision/recall (~9.3% / ~32.1%) documented in SCHEMA.md, not implied
    to be a validated call.
  - `cliff_score`/`ipa_cliff_score`/`cliff_position`/`ipa_annotation_support`
    are unchanged in VALUE but are now documented (SCHEMA.md) as
    ranking/candidate signals, not calls.

Engine Block 2c (correctness fix, see results/engine_block2b_validation.md
Part 4b and results/engine_block2c_validation.md): Block 2b's validation
found `retention_confidence` hard-zeroed whenever `shape_bias_flag` was
True, which inverted rank order and crashed its own AUC (0.548) far below
the underlying `uniformity` signal it is derived from (0.943) -- the
GC/mappability bias gate was validated ONLY as an IPA-cliff discriminator
(results/engine_block2_validation.md Part 3), never as a retention-estimate
gate. `retention_confidence` is now a clean function of `uniformity` alone
(soft-floor ramp at the validated 0.65 operating point) -- `shape_bias_flag`
no longer touches it in any way. `shape_bias_flag` itself is unchanged and
remains available as its own independent column for callers who want to
filter on it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import annotation as annotation_mod
from . import context as context_mod
from . import mask as mask_mod
from . import shape as shape_mod

IPA_CANDIDATE_CLIFF_THRESHOLD = 0.80   # results/ipa_calibration.md Part 2, rule b
RETENTION_UNIFORMITY_FLOOR = 0.65       # results/shape_evidence_final.md operating point


def compute_shape_row(intron_ref, other_mask: np.ndarray, same_mask: np.ndarray,
                       small_mask: np.ndarray, gc_arr, map_arr, adapter, mask_policy: str,
                       polya_hits: dict, alt3_hits: dict) -> dict:
    """All Block-2 columns EXCEPT IRratio_corrected (needs IRratio_classic
    from the Block-1 classic table, merged in afterward by the caller -- see
    `apply_irratio_corrected`) and minor_intron_U12/tier/method (computed
    once per intron / already set by Block 1, not per mask-policy row)."""
    iid, chrom, istart, iend, strand = (intron_ref.intron_id, intron_ref.chrom,
                                         intron_ref.start, intron_ref.end, intron_ref.strand)
    values = adapter.per_base_values(chrom, istart, iend)
    masked_bool = mask_mod.masked_array(other_mask, same_mask, small_mask, mask_policy)
    free_bool = ~masked_bool

    shp = shape_mod.compute_shape(values, free_bool, istart, iend, strand)

    gc_content = context_mod.masked_track_mean(gc_arr, free_bool) if gc_arr is not None else float("nan")
    mappability_score = (context_mod.masked_track_mean(map_arr, free_bool)
                          if map_arr is not None else float("nan"))

    bias_flag, gc_flag = False, False
    if shp["shape_fittable"] and shp["cliff_position"] is not None:
        bias_flag, gc_flag = context_mod.shape_bias_flag(
            gc_arr, map_arr, istart, iend, shp["n_bins"], shp["cliff_position"]
        )

    cliff_score = shp["cliff_score"]
    has_cliff_signal = shp["shape_fittable"] and cliff_score == cliff_score

    ipa_position = shp["cliff_position"]
    polya_positions = list(polya_hits.get(iid, []))
    ann_positions = polya_positions + list(alt3_hits.get(iid, []))
    ann_support = annotation_mod.nearest_within(ann_positions, ipa_position) if ipa_position is not None else False
    polya_near_cliff = annotation_mod.nearest_within(polya_positions, ipa_position) if ipa_position is not None else False

    # ipa_candidate: the single best rule found in results/ipa_calibration.md
    # Part 2 (rule b) -- precision ~9.3%, recall ~32.1% on the confident-
    # long-read-adjudicable long-tier population. HYPOTHESIS-GENERATING
    # candidate nomination, NOT a validated call -- see SCHEMA.md.
    ipa_candidate = bool(has_cliff_signal and cliff_score >= IPA_CANDIDATE_CLIFF_THRESHOLD
                          and not bias_flag and polya_near_cliff)

    retention_confidence = float("nan")
    if shp["shape_fittable"] and shp["uniformity"] == shp["uniformity"]:
        retention_confidence = float(np.clip(
            (shp["uniformity"] - RETENTION_UNIFORMITY_FLOOR) / (1.0 - RETENTION_UNIFORMITY_FLOOR),
            0.0, 1.0,
        ))

    return dict(
        intron_id=iid, mask_policy=mask_policy,
        uniformity=shp["uniformity"], cliff_score=cliff_score, cliff_position=shp["cliff_position"],
        shape_fittable=shp["shape_fittable"], shape_method="full_shape",
        gc_content=gc_content, mappability_score=mappability_score,
        gc_transition_near_shape_feature=gc_flag, shape_bias_flag=bias_flag,
        ipa_cliff_score=cliff_score, ipa_position=ipa_position,
        ipa_annotation_support=ann_support, ipa_candidate=ipa_candidate,
        retention_confidence=retention_confidence,
        polyA_site_present=iid in polya_hits, alt3ss_present=iid in alt3_hits,
    )


def apply_irratio_corrected(irratio_classic, mask_policy: str):
    """SCHEMA.md CORRECTED OUTPUT, Engine Block 2b: IRratio_corrected is
    anchored to the best_perf mask ONLY (null on classic rows, by
    construction) and is now EXACTLY that row's IRratio_classic -- the
    classic-principle ratio (median Ai / (Ai+Ae)) under the validated
    best_perf mask. NO cliff/IPA adjustment is applied (removed in Block 2b
    -- see results/ipa_calibration.md: the best IPA rule found reaches only
    ~9.3% precision, so the prior down-weight corrupted far more retention
    estimates than it fixed). IRratio_classic itself is never modified (kept
    alongside for audit); the two columns are expected to be numerically
    identical on best_perf rows until a validated correction exists."""
    if mask_policy != "best_perf":
        return float("nan")
    if irratio_classic is None or pd.isna(irratio_classic):
        return float("nan")
    return float(irratio_classic)


def compute_u12_candidates(introns: dict, fasta_path: str) -> dict:
    """{iid: bool} minor-intron (U12-type) candidate flag -- a motif-
    similarity proxy (see annotation.py module docstring for why this is not
    a curated-database lookup), computed once per intron (policy-independent)."""
    import pysam

    fa = pysam.FastaFile(fasta_path)
    try:
        out = {}
        for iid, ref in introns.items():
            if ref.strand == "+":
                seq_5ss = fa.fetch(ref.chrom, ref.start - 1, ref.start - 1 + 8)
                seq_3ss_last2 = fa.fetch(ref.chrom, ref.end - 2, ref.end)
            else:
                # transcript direction is genomic 3'->5'; reverse-complement
                raw_5 = fa.fetch(ref.chrom, ref.end - 8, ref.end)
                seq_5ss = _revcomp(raw_5)
                raw_3 = fa.fetch(ref.chrom, ref.start - 1, ref.start + 1)
                seq_3ss_last2 = _revcomp(raw_3)
            out[iid] = annotation_mod.is_u12_candidate(seq_5ss, seq_3ss_last2)
        return out
    finally:
        fa.close()


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]
