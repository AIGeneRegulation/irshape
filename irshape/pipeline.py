"""Engine Packaging: full multi-tier pipeline orchestration.

This wires together modules that already exist and are independently
validated -- classic.py (Block 1), block2.py (Block 2, length_tier=="long"),
block3.py (Block 3, length_tier=="awkward_mid"), annotation.py (PolyASite/
alt-3'/U12) -- into the single per-length-tier sequence the project's numbered
scripts (50/58/64) already run by hand, so `irshape run` can produce one
complete schema-conformant table without the caller knowing about tiers.
No formula from those modules is reproduced or modified here -- this file
only calls them, in the same order and with the same column-merge pattern
those scripts use.

Length-tier cutoffs (100bp, 1000bp) are the frozen SCHEMA.md/routing.py
boundaries, not a tunable -- they are not exposed as CLI options.
"""
from __future__ import annotations

from collections import defaultdict

import pandas as pd

from . import annotation as annotation_mod
from . import block2 as block2_mod
from . import block3 as block3_mod
from . import classic as classic_mod
from . import gtf as gtf_mod
from . import reference as ref_mod
from . import schema as schema_mod

TIER_BOUNDS = {
    "long":        dict(min_length=1000, max_length=None),
    "awkward_mid": dict(min_length=100,  max_length=1000),
    "sub_read":    dict(min_length=None, max_length=100),
}
# Only the long/awkward_mid tiers use the GC/mappability context tracks (the
# multi-bin shape fitter and the two-point donor/acceptor windows both read
# them); sub_read's classic-ratio-only path never touches gc_base/mappability_base
# (see SCHEMA.md's tier routing table), so building it with context would
# just do a genome-wide FASTA/mappability sweep for 9k introns nothing reads.
TIERS_WITH_CONTEXT = {"long", "awkward_mid"}


def _log(progress, msg):
    if progress is not None:
        progress(msg)


def build_tier_bundle(tier: str, gtf: str, genome_build: str, annotation_version: str,
                       cache_dir: str, want_genes, fasta_path: str = None,
                       mappability_bigwig: str = None, mappability_readlen: int = 100,
                       mappability_source: str = None, force: bool = False):
    bounds = TIER_BOUNDS[tier]
    with_context = tier in TIERS_WITH_CONTEXT
    return ref_mod.build_reference(
        gtf=gtf, genome_build=genome_build, annotation_version=annotation_version,
        cache_dir=cache_dir, want_genes=want_genes, force=force,
        min_length=bounds["min_length"], max_length=bounds["max_length"],
        fasta_path=fasta_path if with_context else None,
        mappability_bigwig=mappability_bigwig if with_context else None,
        mappability_readlen=mappability_readlen,
        mappability_source=mappability_source,
    )


def _annotation_context(gtf: str, fasta_path: str, polya_bed: str, bundle, progress=None):
    """PolyASite / GENCODE-alt-3' / U12 lookups for the long tier (Block 2
    Part 4) -- identical sweep to scripts/58_engine_block2b_corrected_run.py.
    Any of gtf/fasta_path/polya_bed may be absent (a minimal bundle); the
    corresponding annotation columns then stay at their "no hit" default
    rather than erroring, since these are corroboration signals, not
    required inputs to the classic/shape computation."""
    exidx = gtf_mod.exon_index(gtf) if gtf else {}
    polya_by_chrom = annotation_mod.load_polya_sites(polya_bed) if polya_bed else {}

    by_chrom = defaultdict(list)
    strand_of = {}
    for iid, ref in bundle.introns.items():
        by_chrom[ref.chrom].append((iid, ref.start, ref.end, ref.gene_id))
        strand_of[iid] = ref.strand
    for c in by_chrom:
        by_chrom[c].sort(key=lambda t: t[1])

    polya_hits = annotation_mod.sweep_polya(by_chrom, polya_by_chrom, strand_of) if polya_bed else {}
    alt3_hits = annotation_mod.sweep_alt3(by_chrom, exidx) if gtf else {}
    _log(progress, f"introns with PolyASite site inside: {len(polya_hits)}")
    _log(progress, f"introns with GENCODE alt-3' site inside: {len(alt3_hits)}")

    u12_flags = {}
    if fasta_path:
        u12_flags = block2_mod.compute_u12_candidates(bundle.introns, fasta_path)
        _log(progress, f"U12 candidates flagged: {sum(u12_flags.values())}/{len(u12_flags)}")
    return polya_hits, alt3_hits, u12_flags


_BLOCK2_COLS = ["uniformity", "cliff_score", "cliff_position", "shape_fittable", "shape_method",
                "gc_content", "mappability_score", "gc_transition_near_shape_feature", "shape_bias_flag",
                "ipa_cliff_score", "ipa_position", "ipa_annotation_support", "ipa_candidate",
                "retention_confidence", "polyA_site_present", "alt3ss_present", "minor_intron_U12"]

_BLOCK3_COLS = ["two_point_ratio", "shape_fittable", "shape_method",
                "gc_content", "mappability_score", "gc_transition_near_shape_feature",
                "shape_bias_flag", "retention_confidence"]


def _long_tier_table(bundle, adapter, gtf, fasta_path, polya_bed, progress=None) -> pd.DataFrame:
    classic_table = classic_mod.compute_table(bundle, adapter, progress=progress)
    polya_hits, alt3_hits, u12_flags = _annotation_context(gtf, fasta_path, polya_bed, bundle, progress)

    rows = []
    n = 0
    for iid, ref in bundle.introns.items():
        om, sm, small = bundle.other_mask[iid], bundle.same_mask[iid], bundle.small_mask[iid]
        gc_arr = bundle.gc_base.get(iid)
        map_arr = bundle.mappability_base.get(iid)
        for policy in ("classic", "best_perf"):
            row = block2_mod.compute_shape_row(ref, om, sm, small, gc_arr, map_arr, adapter, policy,
                                                polya_hits, alt3_hits)
            row["minor_intron_U12"] = u12_flags.get(iid, False)
            rows.append(row)
        n += 1
        if progress is not None and n % 10000 == 0:
            progress(f"  [long] shape/context/ipa: {n}/{len(bundle.introns)} introns")

    b2 = pd.DataFrame(rows).set_index(["intron_id", "mask_policy"])
    table = classic_table.set_index(["intron_id", "mask_policy"])
    for col in _BLOCK2_COLS:
        table[col] = b2[col]
    table["IRratio_corrected"] = [
        block2_mod.apply_irratio_corrected(irc, pol)
        for irc, pol in zip(table["IRratio_classic"], table.index.get_level_values("mask_policy"))
    ]
    table["method"] = "classic_ratio+full_shape"
    return table.reset_index()


def _awkward_mid_tier_table(bundle, adapter, progress=None) -> pd.DataFrame:
    classic_table = classic_mod.compute_table(bundle, adapter, progress=progress)
    rows = []
    n = 0
    for iid, ref in bundle.introns.items():
        om, sm, small = bundle.other_mask[iid], bundle.same_mask[iid], bundle.small_mask[iid]
        gc_arr = bundle.gc_base.get(iid)
        map_arr = bundle.mappability_base.get(iid)
        for policy in ("classic", "best_perf"):
            rows.append(block3_mod.compute_twopoint_row(ref, om, sm, small, gc_arr, map_arr, adapter, policy))
        n += 1
        if progress is not None and n % 10000 == 0:
            progress(f"  [awkward_mid] two-point/context: {n}/{len(bundle.introns)} introns")

    b3 = pd.DataFrame(rows).set_index(["intron_id", "mask_policy"])
    table = classic_table.set_index(["intron_id", "mask_policy"])
    for col in _BLOCK3_COLS:
        table[col] = b3[col]
    table["IRratio_corrected"] = [
        block2_mod.apply_irratio_corrected(irc, pol)
        for irc, pol in zip(table["IRratio_classic"], table.index.get_level_values("mask_policy"))
    ]
    table["method"] = "classic_ratio+two_point"
    return table.reset_index()


def run_full_engine(gtf: str, genome_build: str, annotation_version: str, cache_dir: str,
                     adapter, want_genes=None, fasta_path: str = None,
                     mappability_bigwig: str = None, mappability_readlen: int = 100,
                     mappability_source: str = None, polya_bed: str = None,
                     progress=None) -> pd.DataFrame:
    """Run Block 1 + Block 2 (long) + Block 3 (awkward_mid) across all three
    length tiers and concatenate into one schema-conformant table -- the
    packaged equivalent of running scripts 50 + 58 + 64 and pasting their
    outputs together. `adapter` is reused as-is for all three tier bundles
    (junction/coverage lookups are keyed by intron_id and don't depend on
    which bundle scope built the reference)."""
    tables = []
    for tier in ("long", "awkward_mid", "sub_read"):
        _log(progress, f"[{tier}] building reference bundle...")
        bundle = build_tier_bundle(
            tier, gtf, genome_build, annotation_version, cache_dir, want_genes,
            fasta_path, mappability_bigwig, mappability_readlen, mappability_source,
        )
        _log(progress, f"[{tier}] {len(bundle.introns)} introns")
        if not bundle.introns:
            continue
        if tier == "long":
            t = _long_tier_table(bundle, adapter, gtf, fasta_path, polya_bed, progress)
        elif tier == "awkward_mid":
            t = _awkward_mid_tier_table(bundle, adapter, progress)
        else:
            t = classic_mod.compute_table(bundle, adapter, progress=progress)
        # cast to the frozen nullable dtypes per-tier, before concatenation --
        # otherwise pandas has to guess a common dtype across tiers for
        # columns that are all-NA in one tier's table (e.g. IPA columns are
        # null for the whole awkward_mid/sub_read tables), which triggers a
        # FutureWarning and risks a silent downcast to plain object/float.
        for name, dtype in schema_mod.COLUMN_DTYPES.items():
            if name in t.columns:
                t[name] = t[name].astype(dtype)
        tables.append(t)

    table = pd.concat(tables, ignore_index=True)
    for name, dtype in schema_mod.COLUMN_DTYPES.items():
        table[name] = table[name].astype(dtype)
    table = table[schema_mod.COLUMN_NAMES]
    return table
