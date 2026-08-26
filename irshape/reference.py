"""Reference-build step: annotated introns -> per-mask-policy clean intron
regions (+ flank exons, + GC/mappability base tracks), as one cached artifact
keyed on (genome_build, annotation_version, gene/length scope[, gc/mappability
readlen]).

This is deliberately separate from the adapters: building it means one GTF
scan + one genome-wide exon-index sweep (+ a genome FASTA/mappability-bigwig
read when GC/mappability context is requested), which is expensive enough
that every engine run should reuse a cached bundle instead of rebuilding it.
Adapters only ever *read* a ReferenceBundle.

Block 2 adds GC content + mappability as first-class per-base reference
tracks (Engine Block 2 Part 1) -- both are static properties of the genome
build (+ a k-mer/read length for mappability), not of any RNA-seq sample, so
they belong in this cache alongside the masks, not recomputed per engine run.
Mappability source: UCSC's mirror of the Hoffman-lab Umap k100 multi-read
mappability track for GRCh38
(https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw,
verified reachable via HEAD before fetching, one-shot login-node wget) --
k=100 is the closest standard Umap k-mer length to this project's real A549
read length (101bp, per results/rescue_characterization_summary.md).
"""
from __future__ import annotations

import hashlib
import os
import pickle
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from . import gtf as gtf_mod
from . import mask as mask_mod

SMALLRNA_TYPES = {"snoRNA", "miRNA"}


@dataclass
class IntronRef:
    intron_id: str
    gene: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    host_biotype: str
    length: int
    length_tier: str
    tag: str
    flank_left: tuple    # (chrom, start, end) 1-based inclusive, or None
    flank_right: tuple   # (chrom, start, end) 1-based inclusive, or None


@dataclass
class ReferenceBundle:
    genome_build: str
    annotation_version: str
    gtf_path: str
    scope_key: str
    introns: dict = field(default_factory=dict)          # iid -> IntronRef
    other_mask: dict = field(default_factory=dict)        # iid -> np.bool_[L]
    same_mask: dict = field(default_factory=dict)         # iid -> np.bool_[L]
    small_mask: dict = field(default_factory=dict)        # iid -> np.bool_[L]
    gc_base: dict = field(default_factory=dict)           # iid -> np.float32[L], 1.0=G/C, 0.0=A/T, NaN=N
    mappability_base: dict = field(default_factory=dict)  # iid -> np.float32[L], 0-1
    mappability_readlen: int = None                        # k-mer length of the mappability track used
    mappability_source: str = None                          # provenance string

    def cache_key(self) -> str:
        return f"{self.genome_build}__{self.annotation_version}__{self.scope_key}"

    def has_context(self) -> bool:
        return bool(self.gc_base) or bool(self.mappability_base)


def length_tier(length: int) -> str:
    """sub_read <100bp / awkward_mid 100-1000bp / long >=1000bp.
    Thresholds per SCHEMA.md / results/rescue_characterization_summary.md."""
    if length < 100:
        return "sub_read"
    if length < 1000:
        return "awkward_mid"
    return "long"


def _scope_key(want_genes, min_length, max_length, with_context, mappability_readlen) -> str:
    if want_genes is None:
        base = "genome_wide"
    else:
        h = hashlib.sha1("\n".join(sorted(want_genes)).encode()).hexdigest()[:12]
        base = f"genes_{len(want_genes)}_{h}"
    if min_length:
        base += f"__minlen{min_length}"
    if max_length:
        base += f"__maxlen{max_length}"
    if with_context:
        base += f"__ctx_k{mappability_readlen}"
    return base


def cache_path(cache_dir: str, genome_build: str, annotation_version: str, scope_key: str) -> str:
    return os.path.join(cache_dir, f"irshape_ref__{genome_build}__{annotation_version}__{scope_key}.pkl")


def build_reference(gtf: str, genome_build: str, annotation_version: str,
                     cache_dir: str, want_genes=None, protein_coding_only: bool = True,
                     force: bool = False, min_length: int = None, max_length: int = None,
                     fasta_path: str = None, mappability_bigwig: str = None,
                     mappability_readlen: int = 100,
                     mappability_source: str = None) -> ReferenceBundle:
    """Build (or load from cache) a ReferenceBundle.

    want_genes: restrict the INTRON UNIVERSE to these gene names (None = every
    protein-coding gene, genome-wide). The masking sweep always considers
    exons from ALL genes/biotypes regardless of `want_genes`, since an
    unrelated gene outside the scope can still mask an in-scope intron.

    min_length: drop introns shorter than this (applied AFTER intron
    enumeration, before the mask sweep, so mask/context arrays are never built
    for introns outside scope) -- e.g. Engine Block 2 passes 1000 to scope the
    bundle to length_tier=="long" only, matching Block 1's routing.

    max_length: drop introns `>=` this (same filtering point as min_length).
    Added in Engine Block 3 so the awkward_mid tier (100-999bp) can be scoped
    to exactly its own range -- e.g. min_length=100, max_length=1000 -- without
    also building mask/GC/mappability arrays for the ~116k long-tier introns
    that would otherwise pass a min_length=100-only filter. None (default) =
    no upper bound, matching pre-Block-3 behavior exactly.

    fasta_path / mappability_bigwig: if given, also compute + cache per-base
    GC and mappability tracks (Block 2 Part 1). Omit either to skip that
    track (stays an empty dict on the bundle, matching Block 1's "not
    implemented" placeholder behavior).
    """
    with_context = bool(fasta_path or mappability_bigwig)
    scope_key = _scope_key(want_genes, min_length, max_length, with_context, mappability_readlen)
    path = cache_path(cache_dir, genome_build, annotation_version, scope_key)
    if not force and os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)

    if protein_coding_only:
        all_genes = gtf_mod.gene_spans(gtf, want_types={"protein_coding"})
        scope_genes = set(all_genes) if want_genes is None else (set(want_genes) & set(all_genes))
    else:
        all_genes = gtf_mod.gene_spans(gtf)
        scope_genes = set(all_genes) if want_genes is None else set(want_genes)

    introns_by_gene = gtf_mod.transcript_introns(gtf, want_genes=scope_genes)
    tx_exons = gtf_mod.transcript_exons(gtf, want_genes=scope_genes)
    exidx = gtf_mod.exon_index(gtf)                                   # ALL genes
    smallrna_spans = gtf_mod.gene_type_spans(gtf, SMALLRNA_TYPES)     # ALL genes

    # exon index carrying strand + gene_id, for the clean/anti-near/other-overlap tag
    exon_strand_idx = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            exon_strand_idx[f[0]].append((int(f[3]), int(f[4]), gtf_mod._attr(f[8], "gene_id"), f[6]))
    for c in exon_strand_idx:
        exon_strand_idx[c].sort()

    by_chrom = defaultdict(list)         # chrom -> [(iid,istart,iend,gene_id)]
    by_chrom_strand = defaultdict(list)  # chrom -> [(iid,istart,iend,gene_id,strand)]
    meta = {}
    flanks = {}

    for gene, ints in introns_by_gene.items():
        gdat = tx_exons.get(gene)
        if not gdat:
            continue
        gene_type = all_genes.get(gene, (None,) * 6)[5]
        _, strand, exons, gid = gdat
        for (chrom, istart, iend, strand2, idx) in ints:
            if min_length is not None and (iend - istart + 1) < min_length:
                continue
            if max_length is not None and (iend - istart + 1) >= max_length:
                continue
            iid = f"{gene}:i{idx}:{istart}-{iend}"
            by_chrom[chrom].append((iid, istart, iend, gid))
            by_chrom_strand[chrom].append((iid, istart, iend, gid, strand2))
            meta[iid] = dict(gene=gene, gene_id=gid, chrom=chrom, start=istart, end=iend,
                              strand=strand2, host_biotype=gene_type)
            fl, fr = None, None
            for k in range(len(exons) - 1):
                tx_istart = exons[k][1] + 1
                tx_iend = exons[k + 1][0] - 1
                if tx_istart == istart and tx_iend == iend:
                    fl = (chrom, exons[k][0], exons[k][1])
                    fr = (chrom, exons[k + 1][0], exons[k + 1][1])
                    break
            # 5'/3' orientation: for '-' strand genes the upstream-in-transcript
            # exon is the genomically-later one; flank_left/flank_right below
            # are GENOMIC left/right (matches the flank BED convention already
            # used project-wide in universe_flank.bed), not transcript 5'/3'.
            flanks[iid] = (fl, fr)

    other_mask, same_mask, small_mask, tags = {}, {}, {}, {}
    for chrom, items in by_chrom.items():
        items_sorted = sorted(items, key=lambda t: t[1])
        cls = mask_mod.classify_bases_sorted(items_sorted, exidx.get(chrom, []))
        for iid, (om, sm) in cls.items():
            other_mask[iid] = om
            same_mask[iid] = sm
        sm_small = mask_mod.apply_smallrna(items_sorted, smallrna_spans.get(chrom, []))
        small_mask.update(sm_small)

    for chrom, items in by_chrom_strand.items():
        items_sorted = sorted(items, key=lambda t: t[1])
        tags.update(mask_mod.other_gene_strand_tag(items_sorted, exon_strand_idx.get(chrom, [])))

    introns = {}
    for iid, m in meta.items():
        L = m["end"] - m["start"] + 1
        fl, fr = flanks[iid]
        introns[iid] = IntronRef(
            intron_id=iid, gene=m["gene"], gene_id=m["gene_id"], chrom=m["chrom"],
            start=m["start"], end=m["end"], strand=m["strand"], host_biotype=m["host_biotype"],
            length=L, length_tier=length_tier(L), tag=tags.get(iid, "clean"),
            flank_left=fl, flank_right=fr,
        )

    gc_base, mappability_base = {}, {}
    if fasta_path:
        gc_base = _build_gc_base(fasta_path, introns)
    if mappability_bigwig:
        mappability_base = _build_mappability_base(mappability_bigwig, introns)

    bundle = ReferenceBundle(
        genome_build=genome_build, annotation_version=annotation_version, gtf_path=gtf,
        scope_key=scope_key, introns=introns, other_mask=other_mask, same_mask=same_mask,
        small_mask=small_mask, gc_base=gc_base, mappability_base=mappability_base,
        mappability_readlen=mappability_readlen if mappability_bigwig else None,
        mappability_source=mappability_source if mappability_bigwig else None,
    )
    os.makedirs(cache_dir, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(bundle, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return bundle


def _build_gc_base(fasta_path: str, introns: dict) -> dict:
    """{iid: np.float32[L]} -- 1.0 at G/C bases, 0.0 at A/T, NaN at N/other."""
    import pysam

    fa = pysam.FastaFile(fasta_path)
    try:
        out = {}
        for iid, ref in introns.items():
            seq = fa.fetch(ref.chrom, ref.start - 1, ref.end).upper()
            arr = np.full(len(seq), np.nan, dtype=np.float32)
            b = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
            arr[(b == ord("G")) | (b == ord("C"))] = 1.0
            arr[(b == ord("A")) | (b == ord("T"))] = 0.0
            out[iid] = arr
        return out
    finally:
        fa.close()


def _build_mappability_base(bigwig_path: str, introns: dict) -> dict:
    """{iid: np.float32[L]} -- per-base mappability score (0-1) from a bigWig
    track (e.g. Umap k100)."""
    import pyBigWig

    bw = pyBigWig.open(bigwig_path)
    try:
        chroms = bw.chroms()
        out = {}
        for iid, ref in introns.items():
            if ref.chrom not in chroms:
                out[iid] = np.full(ref.length, np.nan, dtype=np.float32)
                continue
            vals = bw.values(ref.chrom, ref.start - 1, ref.end, numpy=True)
            out[iid] = vals.astype(np.float32)
        return out
    finally:
        bw.close()
