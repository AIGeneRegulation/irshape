"""recount3 input adapter -- FUNCTIONAL in Engine Block 1 (the scale path).

Two local, pre-fetched inputs only (never network, never the snaptron REST
API in a loop -- see CLAUDE.md "Network: file-fetch vs API-loop"):
  1. one or more recount3 base-sum bigWigs (already downloaded on the login
     node) -- per-base coverage, read directly via pyBigWig (NOT megadepth:
     megadepth has no median op, and `intron_abundance_Ai` is specified as a
     *median* over masked positions, so this adapter needs per-base values,
     not per-interval means).
  2. a pre-parsed junction table, produced ONCE by
     `parse_recount3_junction_flatfiles` from a study's own
     `sra.junctions.<SRP>.ALL.{ID,RR,MM}.gz` flat files (never the REST API).
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass

import numpy as np

from .base import CoverageAdapter, JunctionCounts

try:
    import pyBigWig
except ImportError:  # pragma: no cover - surfaced at call time, not import time
    pyBigWig = None


class Recount3Adapter(CoverageAdapter):
    """bigwig_paths: local paths to one or more replicate bigWigs, pooled by
    elementwise mean (NaN-safe: a position is NaN only if every replicate
    lacks data there). junction_counts_by_intron: {intron_id: JunctionCounts},
    typically the output of `parse_recount3_junction_flatfiles`."""

    def __init__(self, bigwig_paths, junction_counts_by_intron: dict):
        if pyBigWig is None:
            raise ImportError("pyBigWig is required for Recount3Adapter (pip install pyBigWig)")
        if not bigwig_paths:
            raise ValueError("Recount3Adapter needs at least one bigwig path")
        self.bigwig_paths = list(bigwig_paths)
        self._handles = [pyBigWig.open(p) for p in self.bigwig_paths]
        self.junction_counts_by_intron = junction_counts_by_intron

    def close(self):
        for h in self._handles:
            h.close()
        self._handles = []

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def per_base_values(self, chrom: str, start: int, end: int) -> np.ndarray:
        L = end - start + 1
        per_replicate = np.full((len(self._handles), L), np.nan, dtype=np.float64)
        for i, h in enumerate(self._handles):
            if chrom not in h.chroms():
                continue
            vals = h.values(chrom, start - 1, end, numpy=True)
            per_replicate[i, :] = vals
        with np.errstate(invalid="ignore"):
            pooled = np.nanmean(per_replicate, axis=0)
        return pooled

    def junction_counts(self, intron_id: str, chrom: str, start: int, end: int) -> JunctionCounts:
        return self.junction_counts_by_intron.get(
            intron_id, JunctionCounts(splice_left=0.0, splice_right=0.0, splice_exact=0.0)
        )


def pool_junction_counts(per_rail: dict, iid: str, rail_subset: list) -> JunctionCounts:
    """Aggregate one intron's per-rail splice sums (from
    `parse_recount3_junction_flatfiles(..., return_per_rail=True)`) into a
    pooled JunctionCounts, MEAN over `rail_subset` -- matching the project's
    existing junction_depth() rail-averaging convention (irlib.py). Lets one
    flat-file parse serve multiple pools (e.g. whole-cell/cytosol/nucleus)
    without re-scanning the (large) MM.gz once per pool."""
    roles = per_rail.get(iid, {})

    def mean_over(role):
        vals = [roles.get(role, {}).get(r, 0) for r in rail_subset]
        return (sum(vals) / len(vals)) if vals else 0.0

    return JunctionCounts(splice_left=mean_over("left"), splice_right=mean_over("right"),
                           splice_exact=mean_over("exact"))


def parse_recount3_junction_flatfiles(id_gz: str, rr_gz: str, mm_gz: str,
                                       introns: dict, rails: list,
                                       chunk_size: int = 20_000_000,
                                       progress=None, return_per_rail: bool = False):
    """Parse a recount3 study's own junction flat files (never the snaptron
    REST API -- see CLAUDE.md) into per-intron splice_left/right/exact.

    introns: {intron_id: (chrom, start, end)} (1-based inclusive, snaptron/
        recount3-junction convention -- matches irlib.py's intron coords).
    rails: recount3 rail_id strings (as found in `id_gz`) for ALL samples
        that might be needed -- pass the full superset (e.g. all 4 A549
        fraction rails) even if a given pool only uses a subset; MM.gz is
        scanned once regardless of how many pools will be derived from it.

    Returns {intron_id: JunctionCounts} (MEAN over ALL of `rails` -- the
    single-pool convenience path) by default. If `return_per_rail=True`,
    instead returns {intron_id: {"left": {rail: count}, "right": {...},
    "exact": {...}}} so a caller can derive several pools (different rail
    subsets) from one parse via `pool_junction_counts`.

    An intron absent from the flat files entirely (no junction ever observed
    at its donor, acceptor, or exact coords in this study) gets all-zero
    counts.

    A single junction row can satisfy more than one role for the same intron
    (the exact intron junction IS also "a junction sharing this donor" and
    "a junction sharing this acceptor") -- SpliceExact is a subset counted
    within SpliceLeft/SpliceRight, not exclusive of them, matching IRFinder
    semantics.
    """
    import pandas as pd

    def log(msg):
        if progress is not None:
            progress(msg)

    rail_to_col = {}
    with gzip.open(id_gz, "rt") as fh:
        fh.readline()
        for i, line in enumerate(fh, start=1):
            rail_to_col[line.strip()] = i
    missing = [r for r in rails if r not in rail_to_col]
    if missing:
        raise ValueError(f"rail(s) not found in {id_gz}: {missing}")
    target_cols = {rail_to_col[r] for r in rails}
    log(f"resolved {len(target_cols)} rails to MM columns")

    exact_lookup, donor_lookup, acceptor_lookup = {}, {}, {}
    for iid, (chrom, start, end) in introns.items():
        exact_lookup.setdefault((chrom, start, end), []).append(iid)
        donor_lookup.setdefault((chrom, start), []).append(iid)
        acceptor_lookup.setdefault((chrom, end), []).append(iid)

    row_roles = {}
    with gzip.open(rr_gz, "rt") as fh:
        fh.readline()
        row = 0
        for line in fh:
            row += 1
            f = line.rstrip("\n").split("\t")
            chrom, jstart, jend = f[0], int(f[1]), int(f[2])
            roles = []
            for iid in exact_lookup.get((chrom, jstart, jend), ()):
                roles.append((iid, "exact"))
            for iid in donor_lookup.get((chrom, jstart), ()):
                roles.append((iid, "left"))
            for iid in acceptor_lookup.get((chrom, jend), ()):
                roles.append((iid, "right"))
            if roles:
                row_roles[row] = roles
            if row % 5_000_000 == 0:
                log(f"RR.gz: {row} rows scanned, {len(row_roles)} matched so far")
    log(f"RR.gz done: {row} junctions total, {len(row_roles)} matched to the intron set")
    target_rows = set(row_roles)

    sums = {}   # (intron_id, role, col) -> total count
    n_scanned = 0
    with gzip.open(mm_gz, "rt") as fh:
        fh.readline()
        fh.readline()
        fh.readline()   # dims line, not needed here
        reader = pd.read_csv(fh, sep="\t", header=None, names=["row", "col", "val"],
                              dtype={"row": "int32", "col": "int32", "val": "int64"},
                              chunksize=chunk_size)
        for chunk in reader:
            n_scanned += len(chunk)
            sub = chunk[chunk["col"].isin(target_cols)]
            if not sub.empty:
                sub = sub[sub["row"].isin(target_rows)]
                for r, c, v in sub.itertuples(index=False):
                    for iid, role in row_roles[r]:
                        key = (iid, role, c)
                        sums[key] = sums.get(key, 0) + int(v)
            log(f"MM.gz: {n_scanned} entries scanned, {len(sums)} accumulator entries so far")
    log(f"MM.gz done: {n_scanned} entries scanned")

    col_to_rail = {rail_to_col[r]: r for r in rails}
    if return_per_rail:
        per_rail = {}
        for iid in introns:
            per_rail[iid] = {
                role: {col_to_rail[c]: sums.get((iid, role, c), 0) for c in target_cols}
                for role in ("exact", "left", "right")
            }
        return per_rail

    out = {}
    for iid in introns:
        per_role = {}
        for role in ("exact", "left", "right"):
            vals = [sums.get((iid, role, c), 0) for c in target_cols]
            per_role[role] = (sum(vals) / len(vals)) if vals else 0.0
        out[iid] = JunctionCounts(splice_left=per_role["left"], splice_right=per_role["right"],
                                   splice_exact=per_role["exact"])
    return out
