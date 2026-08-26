"""Per-base intron masking: classify every intronic base by why it would be
excluded from a coverage measurement, then derive both frozen mask policies
from that one classification. See ../SCHEMA.md "Mask policies".

Ported/generalized from scripts/irlib.py:mask_categories_sorted (same
two-pointer sweep), but returns three independent boolean arrays
(other_mask, same_mask, small_mask) instead of a single priority-ordered
category code, because deriving BOTH `classic` and `best_perf` masks (and the
MASK PROVENANCE fraction columns) needs same_mask and other_mask separately
-- irlib's category code collapses "same-gene-only overlap" into the same
bucket as "no overlap at all", which is correct for the single mask irlib.py
computes but loses information the two-policy contract here needs.
"""
from __future__ import annotations

import numpy as np

POLICIES = ("classic", "best_perf")


def classify_bases_sorted(sorted_introns, exons_sorted):
    """sorted_introns: [(iid, istart, iend, own_gene_id), ...] sorted by istart.
    exons_sorted: irshape.gtf.exon_index(gtf)[chrom] (ALL genes, any biotype),
    sorted by start.

    Returns {iid: (other_mask, same_mask) } where each is a numpy bool array
    of length (iend-istart+1), 1-based-inclusive-intron-order (index 0 = istart).
    `small_mask` (snoRNA/miRNA-host) is layered on separately by
    `apply_smallrna` since it needs a different index (gene_id -> is-smallRNA),
    not another sweep over exons_sorted.
    """
    n = len(exons_sorted)
    lo = 0
    out = {}
    for (iid, istart, iend, own_gid) in sorted_introns:
        while lo < n and exons_sorted[lo][1] < istart:
            lo += 1
        L = iend - istart + 1
        other_mask = np.zeros(L, dtype=bool)
        same_mask = np.zeros(L, dtype=bool)
        j = lo
        while j < n and exons_sorted[j][0] <= iend:
            s, e, gid = exons_sorted[j]
            if e >= istart:
                a, b = max(s, istart), min(e, iend)
                if b >= a:
                    if gid == own_gid:
                        same_mask[a - istart: b - istart + 1] = True
                    else:
                        other_mask[a - istart: b - istart + 1] = True
            j += 1
        out[iid] = (other_mask, same_mask)
    return out


def apply_smallrna(sorted_introns, smallrna_spans_sorted):
    """{iid: small_mask} -- bases where the OTHER-gene exon covering them
    belongs to a snoRNA/miRNA gene. Independent sweep keyed by gene span
    (not exon-level), since a snoRNA/miRNA gene record is what matters, not
    its individual exons."""
    n = len(smallrna_spans_sorted)
    lo = 0
    out = {}
    for (iid, istart, iend, _own_gid) in sorted_introns:
        while lo < n and smallrna_spans_sorted[lo][1] < istart:
            lo += 1
        L = iend - istart + 1
        small_mask = np.zeros(L, dtype=bool)
        j = lo
        while j < n and smallrna_spans_sorted[j][0] <= iend:
            s, e, _gid = smallrna_spans_sorted[j]
            if e >= istart:
                a, b = max(s, istart), min(e, iend)
                if b >= a:
                    small_mask[a - istart: b - istart + 1] = True
            j += 1
        out[iid] = small_mask
    return out


def masked_array(other_mask: np.ndarray, same_mask: np.ndarray, small_mask: np.ndarray,
                  policy: str) -> np.ndarray:
    """Boolean array, True = excluded from measurement, under `policy`."""
    if policy == "classic":
        return other_mask | same_mask
    if policy == "best_perf":
        return other_mask & ~(same_mask & ~small_mask)
    raise ValueError(f"unknown mask policy {policy!r}, expected one of {POLICIES}")


def mask_stats(other_mask: np.ndarray, same_mask: np.ndarray, small_mask: np.ndarray,
               policy: str) -> dict:
    """MASK PROVENANCE fractions for one intron under one policy. See
    SCHEMA.md for exact semantics (components are not required to sum to
    masked_fraction)."""
    L = len(other_mask)
    if L == 0:
        return dict(masked_fraction=float("nan"), mask_other_gene_frac=float("nan"),
                    mask_same_gene_frac=float("nan"))
    masked = masked_array(other_mask, same_mask, small_mask, policy)
    return dict(
        masked_fraction=float(masked.sum()) / L,
        mask_other_gene_frac=float((masked & other_mask).sum()) / L,
        mask_same_gene_frac=float(same_mask.sum()) / L,
    )


def other_gene_strand_tag(sorted_introns, exon_strand_index_sorted) -> dict:
    """tag {clean/anti-near/other-overlap}: whether OTHER-gene exons overlapping
    this intron are same-strand (other-overlap) or exclusively opposite-strand
    (anti-near), or absent entirely (clean).

    exon_strand_index_sorted: [(start,end,gene_id,strand), ...] sorted by start,
    ALL genes. sorted_introns: [(iid, istart, iend, own_gene_id, own_strand), ...]
    sorted by istart.
    """
    n = len(exon_strand_index_sorted)
    lo = 0
    out = {}
    for (iid, istart, iend, own_gid, own_strand) in sorted_introns:
        while lo < n and exon_strand_index_sorted[lo][1] < istart:
            lo += 1
        same_strand_hit = False
        opp_strand_hit = False
        j = lo
        while j < n and exon_strand_index_sorted[j][0] <= iend:
            s, e, gid, strand = exon_strand_index_sorted[j]
            if gid != own_gid and e >= istart:
                if strand == own_strand:
                    same_strand_hit = True
                else:
                    opp_strand_hit = True
            j += 1
        if not same_strand_hit and not opp_strand_hit:
            out[iid] = "clean"
        elif same_strand_hit:
            out[iid] = "other-overlap"
        else:
            out[iid] = "anti-near"
    return out
