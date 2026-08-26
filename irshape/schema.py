"""Frozen output column contract for the irshape engine.

See ../SCHEMA.md for full prose documentation (semantics, formulas, thresholds,
the row-per-mask-policy design rationale). This module is the machine-checkable
counterpart: a typed, ordered column list and an `empty_table()` builder so every
block emits a schema-conformant table even for columns it does not yet compute.

Row grain: one row per (intron_id, mask_policy). See SCHEMA.md "Row grain".
"""
from __future__ import annotations

import pandas as pd

# Each entry: (name, dtype, group, computed_by).
# dtype uses pandas nullable extension types throughout so "not yet computed"
# columns are representable (pd.NA) without silently downcasting to object/float.
#
# computed_by: which block first populates real (non-null) values for this
# column. "1" = Engine Block 1. Sentinel 999 = RETIRED: the column stays in
# the frozen contract (name/dtype preserved) but is deliberately never
# populated going forward -- see SCHEMA.md's Block 2b note (ipa_flag /
# ipa_confidence were found to be miscalibrated determinations and were
# retired in favor of the ranking-only ipa_candidate column).
COLUMNS = [
    # ---- IDENTITY (policy-independent, repeated per mask_policy row) ----
    ("intron_id",      "string",  "IDENTITY", 1),
    ("chrom",          "string",  "IDENTITY", 1),
    ("start",          "Int64",   "IDENTITY", 1),
    ("end",            "Int64",   "IDENTITY", 1),
    ("strand",         "string",  "IDENTITY", 1),
    ("gene_id",        "string",  "IDENTITY", 1),
    ("host_biotype",   "string",  "IDENTITY", 1),
    ("length",         "Int64",   "IDENTITY", 1),
    ("length_tier",    "string",  "IDENTITY", 1),   # sub_read / awkward_mid / long

    # ---- CLASSIC (IRFinder-compatible; policy-dependent) ----
    ("intron_abundance_Ai",  "Float64", "CLASSIC", 1),
    ("splice_left",          "Float64", "CLASSIC", 1),
    ("splice_right",         "Float64", "CLASSIC", 1),
    ("splice_exact",         "Float64", "CLASSIC", 1),
    ("coverage_fraction",    "Float64", "CLASSIC", 1),
    ("IRratio_classic",      "Float64", "CLASSIC", 1),
    ("warn_LowCover",        "boolean", "CLASSIC", 1),
    ("warn_LowSplicing",     "boolean", "CLASSIC", 1),
    ("warn_MinorIsoform",    "boolean", "CLASSIC", 1),
    ("warn_NonUniformCover", "boolean", "CLASSIC", 1),
    ("tag",                  "string",  "CLASSIC", 1),   # clean / anti-near / other-overlap

    # ---- MASK PROVENANCE (policy-dependent) ----
    ("mask_policy",           "string",  "MASK_PROVENANCE", 1),  # classic / best_perf
    ("masked_fraction",       "Float64", "MASK_PROVENANCE", 1),
    ("mask_other_gene_frac",  "Float64", "MASK_PROVENANCE", 1),
    ("mask_same_gene_frac",   "Float64", "MASK_PROVENANCE", 1),

    # ---- SHAPE LAYER (declared; Block 2) ----
    ("uniformity",      "Float64", "SHAPE", 2),
    ("gradient_score",  "Float64", "SHAPE", 2),   # UNVALIDATED, see SCHEMA.md
    ("cliff_score",     "Float64", "SHAPE", 2),
    ("cliff_position",  "Int64",   "SHAPE", 2),
    ("two_point_ratio", "Float64", "SHAPE", 2),
    ("shape_method",    "string",  "SHAPE", 2),   # full_shape / two_point / junction_only
    ("shape_fittable",  "boolean", "SHAPE", 2),

    # ---- GC / MAPPABILITY BIAS (computed; Block 2) ----
    ("gc_content",                       "Float64", "GC_MAPPABILITY", 2),
    ("mappability_score",                "Float64", "GC_MAPPABILITY", 2),
    ("gc_transition_near_shape_feature",  "boolean", "GC_MAPPABILITY", 2),
    ("shape_bias_flag",                   "boolean", "GC_MAPPABILITY", 2),

    # ---- IPA (ranking/candidate signals ONLY -- see SCHEMA.md; ipa_flag/
    # ipa_confidence RETIRED in Block 2b, kept for contract compatibility) ----
    ("ipa_cliff_score",         "Float64", "IPA", 2),
    ("ipa_position",            "Int64",   "IPA", 2),
    ("ipa_annotation_support",  "boolean", "IPA", 2),
    ("ipa_candidate",           "boolean", "IPA", 2),   # added Block 2b; see SCHEMA.md
    ("ipa_flag",                "boolean", "IPA", 999),  # RETIRED Block 2b
    ("ipa_confidence",          "Float64", "IPA", 999),  # RETIRED Block 2b

    # ---- CORRECTED OUTPUT (IRratio_corrected/retention_confidence computed
    # Block 2 (long tier); IRratio_shrunk declared, Block 3+) ----
    ("IRratio_corrected",    "Float64", "CORRECTED", 2),
    ("IRratio_shrunk",       "Float64", "CORRECTED", 3),
    ("retention_confidence", "Float64", "CORRECTED", 2),

    # ---- ANNOTATION (declared; Block 2+) ----
    ("minor_intron_U12",   "boolean", "ANNOTATION", 2),
    ("polyA_site_present", "boolean", "ANNOTATION", 2),
    ("alt3ss_present",     "boolean", "ANNOTATION", 2),

    # ---- PROVENANCE (this block; minimal) ----
    ("tier",       "string",  "PROVENANCE", 1),
    ("method",     "string",  "PROVENANCE", 1),
    ("confidence", "Float64", "PROVENANCE", 1),
]

COLUMN_NAMES = [c[0] for c in COLUMNS]
COLUMN_DTYPES = {c[0]: c[1] for c in COLUMNS}
COLUMN_GROUPS = {c[0]: c[2] for c in COLUMNS}
COLUMN_COMPUTED_BY = {c[0]: c[3] for c in COLUMNS}

MASK_POLICIES = ("classic", "best_perf")
LENGTH_TIERS = ("sub_read", "awkward_mid", "long")
SHAPE_METHODS = ("full_shape", "two_point", "junction_only")
TAGS = ("clean", "anti-near", "other-overlap")


def empty_table(n_rows: int = 0) -> pd.DataFrame:
    """A schema-conformant DataFrame: every frozen column present with the
    correct nullable dtype, filled with n_rows of nulls. Block-1 code fills in
    IDENTITY/CLASSIC/MASK_PROVENANCE/PROVENANCE values; everything else stays
    null, exactly as SCHEMA.md specifies for this block."""
    return pd.DataFrame(
        {name: pd.array([pd.NA] * n_rows, dtype=dtype) for name, dtype in COLUMN_DTYPES.items()},
        columns=COLUMN_NAMES,
    )


def columns_for_block(block: int) -> list[str]:
    """Columns whose computed_by <= block (i.e. real values are expected by then)."""
    return [name for name, _, _, cb in COLUMNS if cb <= block]


def validate(df: pd.DataFrame, block: int | None = None) -> list[str]:
    """Return a list of contract violations (empty list = conformant).

    Checks: all frozen columns present, no unexpected extra columns, dtypes
    match, and (if `block` given) columns computed by an earlier-or-equal
    block contain no nulls."""
    problems = []
    missing = [c for c in COLUMN_NAMES if c not in df.columns]
    if missing:
        problems.append(f"missing columns: {missing}")
    extra = [c for c in df.columns if c not in COLUMN_DTYPES]
    if extra:
        problems.append(f"unexpected columns not in the frozen contract: {extra}")
    for name in COLUMN_NAMES:
        if name not in df.columns:
            continue
        want = COLUMN_DTYPES[name]
        got = str(df[name].dtype)
        if got != want:
            problems.append(f"column {name}: expected dtype {want}, got {got}")
    if block is not None:
        for name in columns_for_block(block):
            if name in df.columns and df[name].isna().all() and len(df) > 0:
                problems.append(f"column {name} is computed by block {COLUMN_COMPUTED_BY[name]} "
                                 f"(<= requested block {block}) but is entirely null")
    return problems
