"""Engine Block 2 Part 4 -- annotation corroboration for the IPA call, and a
minor-intron (U12-type) candidate flag.

PolyASite 2.0 + GENCODE alt-3': same fetch/logic as
scripts/44_rescue_annotation_concordance.py (PolyASite already local,
GENCODE alt-3' derived from the exon index already built for masking -- no
new fetch). `polyA_site_present`/`alt3ss_present` reuse that script's
two-pointer "any site anywhere inside the intron" sweep; `ipa_annotation_support`
additionally requires a site within `context.TRANSITION_WINDOW` bp of the
candidate cliff position, since "anywhere inside a 10kb intron" is too loose
to corroborate a specific cliff call.

minor_intron_U12: the curated U12DB (genome.crg.es/datasets/u12) turned out
to be a legacy interactive CGI search form, not a bulk-downloadable flat
file -- querying it per-intron would be an API loop (CLAUDE.md violation),
and no dump URL was found. Instead this uses a motif-similarity proxy against
the well-established U12-type splice-site consensus (Levine & Durbin 2001;
Sheth et al. 2006): explicitly a computational classifier, NOT a curated
database lookup -- documented as approximate in
results/engine_block2_validation.md, not oversold.
"""
from __future__ import annotations

import gzip
from bisect import bisect_left
from collections import defaultdict

from . import gtf as gtf_mod
from .context import TRANSITION_WINDOW

U12_5SS_CONSENSUS = "GTATCCTT"   # GT-AG subtype, ~98% of U12-type introns
U12_5SS_MAX_MISMATCH = 2         # allow <=2 mismatches (>=75% identity) at positions 3-8


def load_polya_sites(polya_bed_gz: str) -> dict:
    """{chrom: [(pos, strand), ...] sorted by pos} from PolyASite 2.0's
    atlas.clusters bed.gz (col0=chrom digit/X/Y, col3="chrom:pos:strand")."""
    by_chrom = defaultdict(list)
    with gzip.open(polya_bed_gz, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c = f[0]
            if not (c.isdigit() or c in ("X", "Y")):
                continue
            chrom = "chr" + c
            pos = int(f[3].split(":")[1])
            strand = f[5]
            by_chrom[chrom].append((pos, strand))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def sweep_polya(sorted_introns, polya_by_chrom, strand_of) -> dict:
    """{iid: [pos, ...]} of same-strand PolyASite positions strictly inside
    (istart, iend). sorted_introns: [(iid,istart,iend,gene_id)] sorted by istart."""
    out = {}
    for chrom, items in sorted_introns.items():
        sites = polya_by_chrom.get(chrom, [])
        n = len(sites)
        lo = 0
        for (iid, istart, iend, _gid) in items:
            while lo < n and sites[lo][0] < istart:
                lo += 1
            hits = []
            j = lo
            while j < n and sites[j][0] <= iend:
                pos, sstrand = sites[j]
                if sstrand == strand_of[iid]:
                    hits.append(pos)
                j += 1
            if hits:
                out[iid] = hits
    return out


def sweep_alt3(sorted_introns, exon_index_by_chrom) -> dict:
    """{iid: [pos, ...]} of same-gene, other-transcript exon boundaries
    (GENCODE alt-3'/alt-end candidates) strictly inside (istart, iend).
    exon_index_by_chrom: irshape.gtf.exon_index(gtf) output (ALL transcripts,
    gene_id-tagged -- already built for masking, reused here as-is)."""
    out = {}
    for chrom, items in sorted_introns.items():
        exons = exon_index_by_chrom.get(chrom, [])
        n = len(exons)
        lo = 0
        for (iid, istart, iend, own_gid) in items:
            while lo < n and exons[lo][1] < istart:
                lo += 1
            hits = []
            j = lo
            while j < n and exons[j][0] <= iend:
                s, e, gid = exons[j]
                if gid == own_gid:
                    if istart < s < iend:
                        hits.append(s)
                    if istart < e < iend:
                        hits.append(e)
                j += 1
            if hits:
                out[iid] = sorted(set(hits))
    return out


def nearest_within(positions, target, window=TRANSITION_WINDOW) -> bool:
    if not positions or target is None:
        return False
    positions = sorted(positions)
    i = bisect_left(positions, target)
    for j in (i - 1, i):
        if 0 <= j < len(positions) and abs(positions[j] - target) <= window:
            return True
    return False


def is_u12_candidate(seq_5ss: str, seq_3ss_last2: str) -> bool:
    """seq_5ss: first 8 intron bases (5' donor side, transcript direction).
    seq_3ss_last2: last 2 intron bases (3' acceptor side, transcript direction).

    AT-AC boundary (seq_5ss[:2]=="AT" and seq_3ss_last2=="AC") is treated as a
    strong U12 candidate on its own -- major-class (U2-type) AT-AC introns are
    exceedingly rare. Otherwise, a GT-AG intron is a (weaker) candidate if its
    first 8bp is within U12_5SS_MAX_MISMATCH of the U12_5SS_CONSENSUS."""
    if len(seq_5ss) < 8 or len(seq_3ss_last2) < 2:
        return False
    seq_5ss = seq_5ss.upper()
    seq_3ss_last2 = seq_3ss_last2.upper()
    if seq_5ss[:2] == "AT" and seq_3ss_last2 == "AC":
        return True
    if seq_5ss[:2] == "GT":
        mismatches = sum(1 for a, b in zip(seq_5ss, U12_5SS_CONSENSUS) if a != b)
        if mismatches <= U12_5SS_MAX_MISMATCH:
            return True
    return False
