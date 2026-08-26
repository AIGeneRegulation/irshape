"""GENCODE GTF parsing.

Ported from the project's scripts/irlib.py (same logic: best-transcript selection
MANE_Select > Ensembl_canonical > longest, 1-based inclusive intron coordinates
= upstream_exon_end+1 .. downstream_exon_start-1), generalized to take an
explicit `gtf` path instead of a hardcoded project path, and stripped of
project-specific network/A549 helpers (those live in irshape.adapters, not here).
"""
from __future__ import annotations

import re
from collections import defaultdict


def _attr(field: str, key: str):
    m = re.search(key + r' "([^"]+)"', field)
    return m.group(1) if m else None


def gene_spans(gtf: str, want_genes=None, want_types=None):
    """{gene_name: (chrom, start, end, strand, gene_id, gene_type)}.

    want_types: restrict to a set of GENCODE gene_type values (None = all)."""
    out = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            gn = _attr(f[8], "gene_name")
            if want_genes is not None and gn not in want_genes:
                continue
            gt = _attr(f[8], "gene_type")
            if want_types is not None and gt not in want_types:
                continue
            out[gn] = (f[0], int(f[3]), int(f[4]), f[6], _attr(f[8], "gene_id"), gt)
    return out


def transcript_introns(gtf: str, want_genes=None):
    """{gene_name: [(chrom, istart, iend, strand, intron_idx_tx_order), ...]}
    for each gene's best transcript (MANE_Select > Ensembl_canonical > longest)."""
    tx_exons = defaultdict(list)
    tx_info = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] not in ("exon", "transcript"):
                continue
            gn = _attr(f[8], "gene_name")
            if want_genes is not None and gn not in want_genes:
                continue
            tid = _attr(f[8], "transcript_id")
            if f[2] == "transcript":
                tag = f[8]
                score = 2 if "MANE_Select" in tag else (1 if "Ensembl_canonical" in tag else 0)
                tx_info[tid] = (gn, f[0], f[6], score, int(f[4]) - int(f[3]))
            else:
                tx_exons[tid].append((int(f[3]), int(f[4])))
    best = {}
    for tid, (gn, chrom, strand, score, length) in tx_info.items():
        key = (score, length)
        if gn not in best or key > best[gn][0]:
            best[gn] = (key, tid)
    introns = {}
    for gn, (_, tid) in best.items():
        ex = sorted(tx_exons.get(tid, []))
        if len(ex) < 2:
            continue
        chrom = tx_info[tid][1]
        strand = tx_info[tid][2]
        ints = []
        for i in range(len(ex) - 1):
            istart = ex[i][1] + 1
            iend = ex[i + 1][0] - 1
            if iend < istart:
                continue
            idx = i + 1 if strand == "+" else (len(ex) - 1 - i)
            ints.append((chrom, istart, iend, strand, idx))
        introns[gn] = ints
    return introns


def transcript_exons(gtf: str, want_genes=None):
    """{gene: (chrom, strand, [(start,end),...] sorted, gene_id)} for each
    gene's best transcript (same selection as transcript_introns)."""
    tx_exons = defaultdict(list)
    tx_info = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] not in ("exon", "transcript"):
                continue
            gn = _attr(f[8], "gene_name")
            if want_genes is not None and gn not in want_genes:
                continue
            tid = _attr(f[8], "transcript_id")
            if f[2] == "transcript":
                score = 2 if "MANE_Select" in f[8] else (1 if "Ensembl_canonical" in f[8] else 0)
                tx_info[tid] = (gn, f[0], f[6], score, int(f[4]) - int(f[3]), _attr(f[8], "gene_id"))
            else:
                tx_exons[tid].append((int(f[3]), int(f[4])))
    best = {}
    for tid, (gn, chrom, strand, score, length, gid) in tx_info.items():
        key = (score, length)
        if gn not in best or key > best[gn][0]:
            best[gn] = (key, tid)
    out = {}
    for gn, (_, tid) in best.items():
        chrom, strand, gid = tx_info[tid][1], tx_info[tid][2], tx_info[tid][5]
        out[gn] = (chrom, strand, sorted(tx_exons.get(tid, [])), gid)
    return out


def exon_index(gtf: str):
    """{chrom: [(start,end,gene_id), ...] sorted by start} across ALL genes
    (any biotype) -- the mask needs to see exons from genes outside `want_genes`
    too, since an intron can be masked by an unrelated overlapping gene."""
    by = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            by[f[0]].append((int(f[3]), int(f[4]), _attr(f[8], "gene_id")))
    for c in by:
        by[c].sort()
    return by


def gene_type_spans(gtf: str, want_types):
    """{chrom: [(start,end,gene_id), ...] sorted} restricted to `want_types`
    (a set of gene_type strings, e.g. {"snoRNA","miRNA"})."""
    by = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "gene":
                continue
            gt = _attr(f[8], "gene_type")
            if gt not in want_types:
                continue
            by[f[0]].append((int(f[3]), int(f[4]), _attr(f[8], "gene_id")))
    for c in by:
        by[c].sort()
    return by


def to_bam_chrom(chrom: str) -> str:
    """GENCODE 'chrN' -> Ensembl-style BAM contig naming ('N', 'MT')."""
    if chrom == "chrM":
        return "MT"
    return chrom[3:] if chrom.startswith("chr") else chrom
