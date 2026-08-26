"""Shared adapter interface: irshape.classic (and later blocks' shape/GC
layers) consume coverage and junction data through this interface only, so
the same downstream code runs unchanged whether the source is a BAM (per-read,
Block-2+ path) or precomputed recount3 bigWig+junction-flatfile products (the
scale path, functional from Block 1).

Both quantities are returned RAW (unmasked, un-aggregated-by-policy) --
irshape.mask / irshape.classic apply the two mask policies afterward. This
keeps the adapter's job small and identical across input types: "give me the
per-base coverage over this genomic range" and "give me the three junction
counts for this intron."
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np


@dataclass
class JunctionCounts:
    splice_left: float    # sum of junction depth over ALL junctions sharing this intron's donor coord
    splice_right: float   # sum of junction depth over ALL junctions sharing this intron's acceptor coord
    splice_exact: float   # junction depth at the exact (start,end) pair


class CoverageAdapter(ABC):
    """Per-intron coverage + junction access, pooled across whatever
    replicates/samples this adapter was constructed with (mean across
    replicates at the per-base level, matching project convention -- see
    scripts/35_shape_shortread_metrics.py's per-SRR-then-mean pattern)."""

    @abstractmethod
    def per_base_values(self, chrom: str, start: int, end: int) -> np.ndarray:
        """1-based inclusive [start, end] -> float array of length
        (end - start + 1), raw per-base coverage (pooled mean across
        replicates), NaN for positions with no data from any replicate."""
        raise NotImplementedError

    @abstractmethod
    def junction_counts(self, intron_id: str, chrom: str, start: int, end: int) -> JunctionCounts:
        """Pooled (mean across replicate rails/samples) junction counts for
        this intron's donor/acceptor/exact-match junctions."""
        raise NotImplementedError
