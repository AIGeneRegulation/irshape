# irshape column contract (frozen at Engine Block 1)

This document is the frozen output schema for the `irshape` engine. It is the contract
every later block (shape, GC/mappability, IPA, shrinkage) builds against — column
*names*, *types*, and *semantics* declared here do not change without a new schema
version. Columns can be added in later blocks; existing ones are not renamed or
repurposed.

Corresponding typed definition: `schema.py` (`COLUMNS`, `empty_table()`).

## Row grain: one row per (intron_id, mask_policy)

The engine emits **two rows per intron** — one for `mask_policy="classic"`, one for
`mask_policy="best_perf"` — rather than wide columns per policy. This was chosen because:

- Every coverage-derived quantity (`intron_abundance_Ai`, `IRratio_classic`,
  `coverage_fraction`, the `warn_*` flags, `masked_fraction` and its breakdown) is
  policy-dependent — computed from a *different set of masked/unmasked bases* depending
  on which mask is active. Wide-format twin columns (`IRratio_classic_classic`,
  `IRratio_classic_bestperf`, ...) would double every coverage column and get worse, not
  better, as more mask-dependent tiers (shape, GC bias) are added in later blocks.
  Long format keeps one column name per *quantity* and lets `mask_policy` be an ordinary
  filter/group-by key.
- `IDENTITY` columns are policy-independent and are simply repeated across the two rows
  for that `intron_id` (small, deliberate redundancy — trivial to `drop_duplicates` on
  `intron_id` if a caller only wants identity fields).
- Downstream default consumption is: `df[df.mask_policy == "best_perf"]` for the
  corrected/production path, `df[df.mask_policy == "classic"]` for audit/comparison.

`intron_id` + `mask_policy` together are the primary key of the table.

## Mask policies

Both policies mask *other-gene* exon overlap. They differ only in whether a base that is
covered by BOTH an other-gene exon AND a same-gene (alternative-isoform) exon gets masked.
Own-gene-only overlap (no other-gene exon at that base) is **never** masked under either
policy — an alternative-isoform exon of the intron's own gene does not erase intronic
signal for this measurement.

| policy | masks other-gene-only overlap | masks same-gene-only overlap | masks other+same-gene coincident overlap | masks snoRNA/miRNA-host overlap |
|---|---|---|---|---|
| `classic` (IRFinder-original, audit) | yes | **yes** | yes | yes |
| `best_perf` (default, corrected) | yes | no | **no** | yes |

`best_perf` is the mask validated in `results/mask_characterization_summary.md`
(unmasking same-gene-overlap positions recovers 49–53/311 previously-unmeasurable
RETENTION introns, AUC(uniformity, retention-vs-IPA) = 0.9813 on the recovered set, zero
false-inflation against A549 long-read truth on the 54 checkable cases). `classic` is the
literal traditional IRFinder mask (mask any exonic overlap regardless of gene of origin) —
known from Part A of this project to destroy the CLK1 intron-4 retention signal, kept here
**only** as an audit/comparison baseline, not as a candidate default.

Per-base classification used to build both masks (`irshape/mask.py`, ported from
`scripts/irlib.py:mask_categories_sorted`):
- `other_mask`: an exon from a gene ≠ this intron's own gene overlaps this base.
- `same_mask`: an exon from this intron's own gene (i.e. an alt-isoform exon) overlaps
  this base.
- `small_mask`: `other_mask` is true and the covering other-gene exon's `gene_type` is
  `snoRNA` or `miRNA`.

```
classic_masked   = other_mask OR same_mask
best_perf_masked = other_mask AND NOT (same_mask AND NOT small_mask)
```

## Column groups

For each column: **type**, **unit/range**, **computed by** (which engine block first
populates it; `null` = declared here, not yet computed), and semantics.

### IDENTITY (this block; policy-independent, repeated per mask_policy row)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `intron_id` | str | `GENE:i<idx>:<start>-<end>` | Block 1 | matches `results/intron_universe.tsv` convention |
| `chrom` | str | `chrN` | Block 1 | |
| `start` | int | 1-based inclusive genomic coord | Block 1 | = upstream exon end + 1 |
| `end` | int | 1-based inclusive genomic coord | Block 1 | = downstream exon start − 1 |
| `strand` | str | `+`/`-` | Block 1 | |
| `gene_id` | str | Ensembl gene ID | Block 1 | |
| `host_biotype` | str | GENCODE `gene_type` of the intron's own gene | Block 1 | e.g. `protein_coding` |
| `length` | int | bp, ≥ 1 | Block 1 | `end - start + 1` |
| `length_tier` | enum | `sub_read` / `awkward_mid` / `long` | Block 1 | `sub_read`: length < 100; `awkward_mid`: 100 ≤ length < 1000; `long`: length ≥ 1000. Thresholds from `results/rescue_characterization_summary.md` (empirical A549 read length 101bp; IPA essentially absent < 150bp, dominant > 1000bp). |

### CLASSIC — IRFinder-compatible (this block; policy-dependent)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `intron_abundance_Ai` | float | mean depth (reads/bp-equivalent), ≥ 0 or NaN | Block 1 | **median** per-base coverage over this row's masked-free (i.e. retained/unmasked) intronic positions — "robust intron depth" per CLAUDE.md; NaN if zero free positions |
| `splice_left` | float | mean junction count over the pool's rails, ≥ 0 | Block 1 | sum of ALL junction depths sharing this intron's donor coordinate (`start`), any acceptor — "any splicing leaving this 5' site" |
| `splice_right` | float | mean junction count over the pool's rails, ≥ 0 | Block 1 | sum of ALL junction depths sharing this intron's acceptor coordinate (`end`), any donor |
| `splice_exact` | float | mean junction count over the pool's rails, ≥ 0 | Block 1 | junction depth at the exact `(start,end)` pair — the canonical intron-excision junction |
| `coverage_fraction` | float | 0–1 or NaN | Block 1 | fraction of this row's unmasked intronic bases with per-base coverage > 0 (breadth of coverage); NaN if zero free positions |
| `IRratio_classic` | float | 0–1 or NaN | Block 1 | `Ai / (Ai + Ae)` where `Ae = max(splice_left, splice_right)`; NaN if `Ai` is NaN and `Ae` = 0 |
| `warn_LowCover` | bool | | Block 1 | `True` if `Ai < 5` **and** `Ae < 5` (both intron depth and total splicing evidence are low — insufficient data to trust the ratio either way), **or** if flanking-exon coverage < 10 (CLAUDE.md project-wide mandatory depth gate, "informative reads ≥ 20 AND flanking-exon coverage ≥ 10"; flank coverage is computed internally for this flag but is not itself a persisted column in Block 1) |
| `warn_LowSplicing` | bool | | Block 1 | `True` if `Ae < 5` — too little total splice evidence at either site to trust the denominator |
| `warn_MinorIsoform` | bool | | Block 1 | `True` if `splice_exact < 0.3 * max(splice_left, splice_right, 1)` — most splicing observed at this donor/acceptor does **not** use this exact intron, i.e. the annotated intron itself is a minor isoform choice, so `IRratio_classic` conflates retention with alternative splicing |
| `warn_NonUniformCover` | bool | | Block 1 | `True` if the ratio of median coverage between the 5' half and 3' half of the unmasked intron differs by > 3× — a coarse 2-region proxy (the multi-bin SHAPE-layer non-uniformity test is not implemented until Block 2; this flag exists so Block-1 consumers have *some* non-uniformity signal) |
| `tag` | enum | `clean` / `anti-near` / `other-overlap` | Block 1 | `clean`: `mask_other_gene_frac == 0` (no other-gene exon touches this intron at all); `other-overlap`: overlapping other-gene exon(s) share this intron's strand; `anti-near`: all overlapping other-gene exon(s) are on the opposite strand (antisense) |

### MASK PROVENANCE (this block; policy-dependent)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `mask_policy` | enum | `classic` / `best_perf` | Block 1 | row-grain key, see above |
| `masked_fraction` | float | 0–1 | Block 1 | fraction of raw intron length (`length`) excluded from measurement under this row's policy |
| `mask_other_gene_frac` | float | 0–1 | Block 1 | fraction of raw intron length masked where an other-gene exon overlaps (component; components can co-occur at the same base, so `mask_other_gene_frac + mask_same_gene_frac` is not required to equal `masked_fraction`) |
| `mask_same_gene_frac` | float | 0–1 | Block 1 | fraction of raw intron length where a same-gene (alt-isoform) exon overlaps, regardless of whether that base is actually masked under this row's policy (provenance, not a masked/free indicator by itself) |

### SHAPE LAYER (declared here; Block 2)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `uniformity` | float | 0–1, higher = more uniform | Block 2 | `1 - CV` of binned coverage over the unmasked intron |
| `gradient_score` | float | unitless, signed | Block 2 | **UNVALIDATED** — slope of bin coverage vs. position / mean coverage; carry the flag through to any consumer until Block 2 validates it |
| `cliff_score` | float | 0–1 | Block 2 | largest single-step relative coverage drop between adjacent valid bins |
| `cliff_position` | int | genomic 1-based coord | Block 2 | midpoint of the cliff-defining bin boundary |
| `two_point_ratio` | float | 0–1 | Block 2 | donor-window / (donor+acceptor-window) coverage share, for `awkward_mid` introns (see `rescue_characterization_summary.md` Part 3) |
| `shape_method` | enum | `full_shape` / `two_point` / `junction_only` | Block 2 | routing is **determined** in Block 1 from `length_tier` (see below) but the column itself is null until Block 2 actually runs the corresponding fitter |
| `shape_fittable` | bool | | Block 2 | whether the chosen `shape_method` produced a usable result for this intron (enough valid bins / unmasked windows) |

Block-1 routing rule (computed for reporting in `engine_block1_validation.md`, **not**
persisted as the `shape_method` column value until Block 2):
`sub_read → junction_only`, `awkward_mid → two_point`, `long → full_shape`.

### GC / MAPPABILITY BIAS (declared here; Block 2/3)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `gc_content` | float | 0–1 | Block 2/3 | GC fraction of the unmasked intron sequence |
| `mappability_score` | float | 0–1 | Block 2/3 | no local mappability/repeat track exists yet (confirmed absent in `mask_characterization_summary.md` Part 2, "P4... not separately computable"); requires acquiring one |
| `gc_transition_near_shape_feature` | bool | | Block 2/3 | whether a GC/mappability transition coincides with a detected shape feature (cliff, etc.) |
| `shape_bias_flag` | bool | | Block 2/3 | `True` when `gc_transition_near_shape_feature` — flags a shape feature as possibly technical, not biological |

### IPA (Block 2; **RANKING/CANDIDATE SIGNALS ONLY as of Block 2b — see addendum below,
none of these is a validated classification**)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `ipa_cliff_score` | float | 0–1 | Block 2 | RANKING signal, identical to `cliff_score`. NOT a calibrated score — no threshold on it reaches usable precision (see Block 2b addendum) |
| `ipa_position` | int | genomic 1-based coord | Block 2 | candidate internal poly(A)/alt-3'-end position (wherever the ranked cliff sits) |
| `ipa_annotation_support` | bool | | Block 2 | PolyASite 2.0 / GENCODE alt-3' site within 500bp of `ipa_position` — a corroboration SIGNAL, not by itself a call |
| `ipa_candidate` | bool | | Block 2b | **HYPOTHESIS-GENERATING candidate nomination, NOT a validated IPA call.** `cliff_score>=0.80 AND NOT shape_bias_flag AND` a PolyASite site within 500bp of `ipa_position`. Measured on the confident-long-read-adjudicable long-tier population (`results/ipa_calibration.md`): **precision ≈0.093, recall ≈0.321** — i.e. roughly 9 in 10 candidates are NOT actually IPA by long-read. Use only to rank/prioritize candidates for follow-up, never to assert an intron IS IPA |
| `ipa_flag` | bool | | **RETIRED (Block 2b)** | always null. Was `cliff_score>=0.80 AND NOT shape_bias_flag`; found to be 96.7% false positive against confident long-read labels (`results/ipa_calibration.md` Part 4). Column kept only for schema-contract compatibility — do not read it |
| `ipa_confidence` | float | 0–1 | **RETIRED (Block 2b)** | always null. Was a 3-level (0.9/0.6/0.3) encoding whose EMPIRICAL precision was 6.8%/1.7%/1.4% respectively — not remotely calibrated to its nominal values. Column kept only for schema-contract compatibility — do not read it |

### CORRECTED OUTPUT (Block 2, long tier; Block 3+ elsewhere)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `IRratio_corrected` | float | 0–1 or NaN | Block 2 (long tier) | **default production metric.** As of Block 2b: EQUALS the `best_perf`-mask row's `IRratio_classic` (median `Ai`/(Ai+Ae) under the validated best_perf mask) — **no IPA/cliff adjustment is applied** (removed; see Block 2b addendum for why). `mask_policy != "best_perf"` rows are null by construction |
| `IRratio_shrunk` | float | 0–1 or NaN | Block 3+ | `IRratio_corrected` after empirical-Bayes shrinkage toward a length/coverage-matched prior |
| `retention_confidence` | float | 0–1 | Block 2 (long tier), fixed Block 2c | soft-floor ramp on `uniformity` ONLY (see Block 2c addendum) — no longer zeroed by `shape_bias_flag`. Does NOT depend on any IPA column or on `shape_bias_flag` |

### ANNOTATION (declared here; Block 2+)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `minor_intron_U12` | bool | | Block 2+ | U12-type (minor spliceosome) intron, from splice-site motif/annotation |
| `polyA_site_present` | bool | | Block 2+ | PolyASite 2.0 / GENCODE site anywhere inside the intron (locus-level prior, not per-sample; see `rescue_characterization_summary.md` Part 4 caveat) |
| `alt3ss_present` | bool | | Block 2+ | an alternative 3' splice site from another annotated transcript falls inside this intron |

### PROVENANCE (this block; minimal, expanded later)

| column | type | unit/range | computed by | notes |
|---|---|---|---|---|
| `tier` | enum | mirrors `length_tier` for now (`sub_read`/`awkward_mid`/`long`) | Block 1 | which measurement regime this intron falls into; kept distinct from `length_tier` because later blocks may re-tier on evidence quality, not just length |
| `method` | str | `"classic_ratio"` (Block 1); more values added in Block 2+ | Block 1 | which computation actually produced this row's headline ratio, so downstream confidence is legible without re-deriving it |
| `confidence` | float | 0–1 | Block 1 (coarse) | Block 1 sets this from `warn_*` flags only (`1.0` if no warnings, `0.5` if exactly one, `0.0` if ≥ 2) as a placeholder; Block 3's `retention_confidence` supersedes this once corrections exist |

## Adapter input contract

Both input adapters (`bam_adapter`, `recount3_adapter`) return the same intermediate
shape before classic-metric computation: per intron, per mask_policy — masked-free
per-base coverage values (for `Ai`/`coverage_fraction`/`warn_NonUniformCover`), flanking
exon coverage (for the `warn_LowCover` mandatory depth gate), and junction counts
(`splice_left`, `splice_right`, `splice_exact`). See `irshape/adapters/base.py`.

## What Block 1 does NOT compute

Everything in SHAPE LAYER, GC / MAPPABILITY BIAS, IPA, CORRECTED OUTPUT, and ANNOTATION
is declared above with types/ranges/semantics but is emitted as `null`/`NaN`/`pd.NA` in
`results/engine_block1_baseline_A549.tsv`. `schema.py:empty_table()` fills these with the
correct dtype's null value so downstream code can rely on the column existing and typed
correctly, even before it is populated.

## Block 2 addendum: concrete formulas (long tier only, `length_tier=="long"`, >=1000bp)

Block 2 populates the previously-null SHAPE, GC_MAPPABILITY, IPA, CORRECTED, and
ANNOTATION columns, but **only** for `length_tier=="long"` rows (`shape_method` stays
null for `sub_read`/`awkward_mid` rows until Block 3). No column is renamed or
repurposed; this section adds the semantics the Block-1 declaration deferred.

- **`uniformity`** — Shannon entropy of the per-bin coverage distribution over
  `best_perf`/`classic`-masked-free positions, normalized by max entropy (Pielou
  evenness): `H = -sum(p_i * log(p_i))`, `p_i = cov_i / sum(cov)` over valid bins,
  `uniformity = H / log(n_valid_bins)`. Chosen over the earlier `1 - CV` formulation
  (Block-1-era scripts) because it's bounded in [0,1] without clipping and degrades
  gracefully with zero-coverage bins. Binning: adaptive `NBIN = clip(L//15, 4, 30)`, a
  bin is valid if ≥25% of its raw span survives masking (the LENIENT operating point
  from `results/shape_evidence_final.md`). NaN if <2 valid bins.
- **`cliff_score` / `cliff_position`** — largest single-step relative coverage drop
  between adjacent valid bins, walked 5'→3' (unchanged from the validated historical
  definition). `ipa_cliff_score`/`ipa_position` are currently identical to
  `cliff_score`/`cliff_position` (Block 2 does not yet add IPA-specific shape logic
  beyond reusing the general cliff signal).
- **`shape_fittable`** — `length >= 20` and `n_valid_bins >= 3`.
- **`gc_content` / `mappability_score`** — mean GC fraction / mean Umap k100
  mappability over this row's unmasked bases (mask-policy-dependent, like `Ai`).
- **`gc_transition_near_shape_feature`** — the largest adjacent-bin GC-fraction jump
  (≥0.15 absolute, over the SAME bins as coverage) falls within 500bp of `cliff_position`.
- **`shape_bias_flag`** — `gc_transition_near_shape_feature` **OR** a mappability drop
  (mean mappability <0.5 in the bin nearest `cliff_position`, within 500bp). Only the
  GC half of this OR has its own named column, per the frozen contract; the
  mappability half is folded into `shape_bias_flag` directly.
- **`ipa_flag` / `ipa_confidence`** — **RETIRED as of Block 2b; see the Block 2b addendum
  below.** (Originally: `ipa_flag = shape_fittable AND cliff_score >= 0.80 AND NOT
  shape_bias_flag`; `ipa_confidence` a 3-level 0.9/0.6/0.3 encoding. Both were found
  badly miscalibrated and are no longer populated.)
- **`ipa_annotation_support`** — a PolyASite 2.0 or GENCODE alt-3' site within 500bp of
  `cliff_position` (same tolerance as the project's established cliff-vs-annotation
  concordance checks). A corroboration signal, not a call.
- **`retention_confidence`** — **see Block 2c addendum below; the `shape_bias_flag`
  hard-zero documented here for Block 2/2b was found to invert rank order and has been
  removed.** As of Block 2c: `clip((uniformity - 0.65) / 0.35, 0, 1)` (0.65 is the
  validated soft floor from `results/shape_evidence_final.md`), unconditionally. NaN if
  not `shape_fittable`. Independent of `shape_bias_flag` and of any IPA column.
- **`IRratio_corrected`** — **see Block 2b addendum below; the formula documented in
  Block 2 (an IPA-driven down-weight) was found to over-correct and has been removed.**
- **`polyA_site_present` / `alt3ss_present`** — any PolyASite 2.0 / GENCODE alt-3' site
  anywhere inside the intron (the "reach" definition, not proximity-gated — see
  `ipa_annotation_support` for the proximity-gated version).
- **`minor_intron_U12`** — a **motif-similarity proxy**, not a curated-database lookup:
  the curated U12DB (genome.crg.es/datasets/u12) turned out to be a legacy interactive
  search form with no bulk-downloadable flat file, so querying it would require an API
  loop (disallowed by CLAUDE.md). Flags AT-AC-boundary introns, or GT-AG introns whose
  first 8bp is within 2 mismatches of the canonical U12 5'SS consensus `GTATCCTT`
  (Levine & Durbin 2001; Sheth et al. 2006). Approximate; not validated against a
  gold-standard U12 set in this block.
- **`method`** — `"classic_ratio+full_shape"` for long-tier rows (vs. Block 1's plain
  `"classic_ratio"`).

See `results/engine_block2_validation.md` for the AUC validation of `uniformity` and
`cliff_score` against pooled long-read labels, and the GC/mappability gate's
false-positive-enrichment check.

## Block 2b addendum: IPA over-correction fix (correctness fix, not a new feature)

Follow-up calibration (`results/ipa_calibration.md`) found that no cliff-score threshold,
up to and including `cliff_score>=0.999`, reaches usable precision against confident
long-read labels (peak ≈7% around 0.94–0.98; annotation corroboration lifts the best case
to ≈9.3% precision at ≈32.1% recall). The Block 2 `ipa_flag` rule (`cliff_score>=0.80 AND
NOT shape_bias_flag`) was **96.7% false positive** on the confident-long-read-adjudicable
long-tier population, and `IRratio_corrected`'s IPA-driven down-weight was firing (and
numerically shrinking the retention estimate) on the great majority of those false
positives. This addendum documents the fix (Engine Block 2b); no column is renamed, and
`IRratio_classic` was never touched.

- **`IRratio_corrected`** — the IPA-driven down-weight is **removed entirely**.
  `IRratio_corrected = IRratio_classic` on `best_perf` rows (i.e. exactly the median-`Ai`,
  best_perf-mask, classic-formula ratio — the two validated Block-1 corrections, mask
  policy and the median estimator, and nothing else); null on `classic` rows, unchanged.
  The two columns are expected to be numerically **identical** on `best_perf` rows until a
  validated correction (e.g. Block 3's shrinkage) exists — this is intentional, not a bug:
  there is currently no validated basis for `IRratio_corrected` to differ from
  `IRratio_classic` on the long tier.
- **`ipa_flag` / `ipa_confidence`** — retired (schema dtype `999`, see `schema.py`).
  Always null starting Block 2b. Kept in the frozen column list only so old code reading
  the schema doesn't break; new code must not read them as determinations.
- **`ipa_candidate`** (new, boolean, IPA group) — the single best rule found in
  calibration: `cliff_score>=0.80 AND NOT shape_bias_flag AND` a PolyASite 2.0 site within
  500bp of `ipa_position`. Measured precision ≈0.093, recall ≈0.321 on the confident-
  long-read-adjudicable long-tier population — explicitly a **hypothesis-generating
  candidate nomination**, not a validated call. Ship `cliff_score`/`ipa_cliff_score` as
  the primary continuous ranking signal; use `ipa_candidate` only to shortlist, never to
  assert.

See `results/engine_block2b_validation.md` for the re-validation (corrected now
equal-or-better than classic in long-read agreement; the previously-over-corrected
introns' estimates restored; `retention_confidence` unaffected; CLK1 intron 4 note).

## Block 2c addendum: decouple `retention_confidence` from `shape_bias_flag` (correctness fix)

Block 2b's re-validation (`results/engine_block2b_validation.md` Part 4b) found
`retention_confidence`'s own AUC (retention-vs-IPA, strict `pooled_cls=="RETENTION"`
population) was **0.548** — far below the `uniformity` signal it is derived from
(**0.943**, bit-identical to Block 2's own reported figure). Root cause: `block2.py` set
`retention_confidence = 0.0` whenever `shape_bias_flag` was True, a hard floor override
(not a monotonic clip) that inverted rank order for genuinely high-uniformity introns
that happened to be GC/mappability bias-flagged. The bias gate was validated ONLY as an
IPA-cliff discriminator (`results/engine_block2_validation.md` Part 3: enrichment of
long-read false positives among bias-flagged high-`cliff_score` introns) — it was never
validated as a retention-estimate gate, and coupling it to `retention_confidence` was
unjustified by any evidence.

- **`retention_confidence`** — the `shape_bias_flag` hard-zero is **removed entirely**.
  Now purely `clip((uniformity - 0.65) / 0.35, 0, 1)` if `shape_fittable` else NaN — a
  clean function of `uniformity` alone. `uniformity` itself is unchanged (not
  recomputed); only the derived `retention_confidence` formula changed.
- **`shape_bias_flag`** — unchanged, independent column. Still available for callers who
  want to filter out GC/mappability-confounded rows themselves; no longer implicitly
  applied inside `retention_confidence`.
- No other column changes. `IRratio_corrected`, `ipa_flag`/`ipa_confidence` (retired),
  `ipa_candidate`, and all Block 2b corrections are untouched by this fix.

See `results/engine_block2c_validation.md` for the re-validation (`retention_confidence`
AUC restored to ~uniformity's own ~0.94/~0.90; everything else confirmed unmoved).

## Block 3 addendum: `length_tier=="awkward_mid"` (100-1000bp), `shape_method="two_point"`

Block 3 populates the previously-null SHAPE/GC_MAPPABILITY/CORRECTED columns for the
`awkward_mid` tier (64,315 introns genome-wide, `results/intron_universe.tsv`) -- the
length range too short/mask-heavy for the multi-bin full-shape fitter (Block 2) but too
long for the classic junction ratio alone to need no shape signal (`sub_read`, which needs
no Block 2/3 shape layer at all -- AUC 0.955 from the classic ratio directly, per
`results/rescue_characterization_summary.md` Part 2). This is the LAST length tier; after
this block every row in `intron_universe.tsv` has a routed, computed `shape_method`.

- **`two_point_ratio`** -- `donor_cov / (donor_cov + acceptor_cov)`, mean `best_perf`/
  `classic`-masked-free coverage in a 30bp window just inside the donor (transcript-5')
  boundary vs. a 30bp window just inside the acceptor (transcript-3') boundary,
  strand-aware. True retention reads as balanced (~0.5); IPA/internal-end-use reads as
  donor-heavy (>0.5) -- see `results/rescue_characterization_summary.md` Part 3 (median
  0.483 RETENTION vs. 0.792 IPA). As stated there, this also matches plain nascent
  5'->3' read-through (same coverage pattern) -- it is a retention-vs-NOT-retention
  signal at this length range, not an IPA-vs-nascent discriminator, and **no IPA call is
  derived from it** (see below).
- **`shape_fittable`** -- both the donor and acceptor windows have >=15bp of actually
  unmasked coverage (`irshape/twopoint.py`'s `MIN_UNMASKED`; window target 30bp,
  `WINDOW`). Unchanged from the validated characterization parameters, so the engine
  reproduces that AUC/coverage rather than a re-tuned one.
- **`gc_content` / `mappability_score`** -- same semantics/computation as Block 2 (mean
  over this row's unmasked bases across the WHOLE intron, not just the two windows --
  these are sequence properties, independent of which shape method applies).
- **`gc_transition_near_shape_feature` / `shape_bias_flag`** -- adapted from Block 2's
  bin-walk-vs-cliff-position gate to this method's two discrete windows: flags when the
  donor and acceptor windows themselves differ by >=0.15 GC fraction, or either window's
  mean mappability is <0.5 (same thresholds as Block 2, `irshape/twopoint.py`'s
  `two_point_bias_flag`) -- i.e. "is there a technical GC/mappability difference between
  the two things `two_point_ratio` is comparing." Per the Block 2c principle, this flag
  is diagnostic only and is **never** consulted by `retention_confidence`.
- **`ipa_cliff_score` / `ipa_position` / `ipa_annotation_support` / `ipa_candidate`** --
  **left null for this tier.** The two-point method has no bin-walk, so there is no cliff
  position to localize an IPA call from, and Block 2b already found the full-shape
  method's OWN cliff signal has no usable IPA precision by itself
  (`results/ipa_calibration.md`) -- there is no evidentiary basis to invent a weaker-tier
  substitute. IPA stays ranking-only where it exists at all (long tier only); this tier
  does not attempt it, consistent with 2b/2c.
- **`retention_confidence`** -- `TWO_POINT_AUC_CAP * clip(1 - 2*|two_point_ratio - 0.5|,
  0, 1)` if `shape_fittable` else NaN, where `TWO_POINT_AUC_CAP = 0.68`
  (`irshape/block3.py`), rounded from this tier's own characterization AUC (0.677,
  `results/rescue_characterization_summary.md` Part 3: AUC(donor_share, IPA vs
  RETENTION) on the awkward-middle panel). The cap is an explicit, deliberate
  ceiling -- even a perfectly-balanced (`two_point_ratio==0.5`) read is never asserted
  with more confidence than this tier's own measured discriminative power, unlike the
  long tier's `uniformity`-based mapping (AUC ~0.94, which needed no such cap since its
  own ceiling already sits close to 1). This makes the two tiers' `retention_confidence`
  values **not directly comparable in magnitude** -- both are calibrated to their own
  tier's validated accuracy, not to a shared numeric scale. Independent of
  `shape_bias_flag`, matching the Block 2c fix for the long tier.
- **`IRratio_corrected`** -- unchanged formula (`irshape/block2.py`'s
  `apply_irratio_corrected`, reused as-is): `IRratio_classic` on `best_perf` rows, null on
  `classic` rows. No shape/two-point-driven adjustment is applied here either, for the
  same reason Block 2b removed one for the long tier -- no validated basis to differ from
  the classic ratio yet.
- **`method`** -- `"classic_ratio+two_point"` for `awkward_mid`-tier rows.

See `results/engine_block3_validation.md` for the AUC re-validation (two_point_ratio vs.
pooled long-read labels, inside the engine), the CLK1 intron 4 regression check (378bp,
now in the `awkward_mid` tier for the first time), and confirmation that the long tier and
2b/2c outputs are untouched by this block.
