"""BAM input adapter -- STUBBED in Engine Block 1.

Per PART 2 of the Block 1 task: both adapters share `CoverageAdapter`, but
only `recount3_adapter` is made functional this block (it's what Part 3's
validation runs against). This stub fixes the intended interface and
constructor signature so Block 2+ can implement it without changing any
downstream caller.

Intended implementation (not built here):
  - pysam.AlignmentFile, one pass per BAM, restricted to the region of
    interest (fetch per gene, like scripts/03_longread_pir.py's stage-2
    gene-batching, not per-intron -- reduces redundant fetches for
    multi-intron genes).
  - Paired-end aware: a fragment's coverage/junction evidence should be
    counted once per fragment, not once per mate; mate pairs need to be
    merged (by query name) before per-base accumulation.
  - Spliced CIGAR (`N` ops) define observed junctions directly -- no
    external junction source needed for a BAM input, unlike the recount3
    path. `splice_left`/`splice_right`/`splice_exact` are all derivable from
    the same one-pass junction tally (group by (donor,acceptor), then filter
    by shared donor / shared acceptor / exact match).
  - Strand-aware via the library's read1/read2 + alignment-strand convention
    (this project's stranded protocols; see irlib.py `to_bam_chrom` for the
    matching contig-naming convention already in use for SG-NEx BAMs).
  - UNIQUE reads only by default: exclude `is_secondary`, `is_supplementary`,
    and (new vs. the current project scripts) filter on a mapping-quality or
    `NH:i:1` tag threshold, configurable.
  - Multi-mapping handling: NOT implemented here. A future block should add
    KMA-style (or similar) fractional multi-mapper redistribution rather than
    dropping multi-mappers outright, since intron-proximal reads are
    disproportionately multi-mapping (retained introns often share sequence
    with unspliced pre-mRNA from paralogs/pseudogenes). Flagged as future
    work, not a Block-1 gap that blocks anything -- the recount3 path (which
    the Block-1 validation actually uses) has already-pooled bigWig coverage
    that is not re-derivable at the read level anyway.
"""
from __future__ import annotations

import numpy as np

from .base import CoverageAdapter, JunctionCounts


class BamAdapter(CoverageAdapter):
    def __init__(self, bam_paths, unique_only: bool = True, min_mapq: int = 1):
        self.bam_paths = list(bam_paths)
        self.unique_only = unique_only
        self.min_mapq = min_mapq

    def per_base_values(self, chrom: str, start: int, end: int) -> np.ndarray:
        raise NotImplementedError(
            "BamAdapter is stubbed in Engine Block 1 -- see module docstring "
            "for the intended pysam-based implementation. Use recount3_adapter "
            "for the functional Block-1 path."
        )

    def junction_counts(self, intron_id: str, chrom: str, start: int, end: int) -> JunctionCounts:
        raise NotImplementedError(
            "BamAdapter is stubbed in Engine Block 1 -- see module docstring."
        )
