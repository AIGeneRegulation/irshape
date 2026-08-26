"""Adaptive intron binning shared by the SHAPE layer (coverage) and the
GC/mappability CONTEXT layer, so coverage bins and GC/mappability bins align
1:1 -- required for the Part 3 "does a GC/mappability transition coincide
with this cliff" gate.

LENIENT operating point from results/shape_evidence_final.md /
scripts/46_mask_recovery_sweep.py: NBIN=clip(L//15,4,30), MIN_RAW_LEN=20,
MIN_BIN_FRAC=0.25, MIN_VALID_BINS=3 -- the setting whose AUC numbers Engine
Block 2 validates against (Part 6).
"""
NBIN_DIV, NBIN_LO, NBIN_HI = 15, 4, 30
MIN_RAW_LEN = 20
MIN_BIN_FRAC = 0.25
MIN_VALID_BINS = 3


def adaptive_nbin(length: int) -> int:
    return max(NBIN_LO, min(NBIN_HI, length // NBIN_DIV))


def bin_ranges(istart: int, iend: int, nbin: int):
    """[(bs, be), ...] genomic 1-based inclusive per-bin ranges spanning
    [istart, iend], nbin bins, in genomic (not strand-oriented) order."""
    length = iend - istart + 1
    ranges = []
    for b in range(nbin):
        bs = istart + (b * length) // nbin
        be = istart + ((b + 1) * length) // nbin - 1
        if b == nbin - 1:
            be = iend
        ranges.append((bs, be))
    return ranges
