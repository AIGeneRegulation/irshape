"""irshape CLI. Two subcommands:
  irshape build-reference  -- build/cache a ReferenceBundle from a GTF
  irshape run              -- run the engine using the recount3 adapter,
                               write a schema-conformant TSV

Reference resolution (Engine Packaging Part 2): both subcommands accept
`--ref-dir <bundle>` (or the `IRSHAPE_REF` env var) pointing at an external,
separately-downloaded reference bundle (GTF, genome FASTA, mappability
bigWig, PolyASite atlas, prebuilt per-tier caches -- see
`build_reference_bundle.py` / SCHEMA.md). Individual
`--gtf`/`--fasta`/... flags still work standalone (no ref-dir needed) and
override whatever a given `--ref-dir` would have supplied.

`irshape run` defaults to the FULL engine (Block 1 classic + Block 2
full-shape + Block 3 two-point, all three length tiers, one merged table --
see `irshape/pipeline.py`); pass `--classic-only` for the lighter Block-1-only
path (single unscoped bundle, no shape/GC/IPA/annotation columns populated).

The bam adapter is not exposed here yet (it is stubbed -- see
adapters/bam_adapter.py); `run --input bam` raises NotImplementedError.
"""
from __future__ import annotations

import argparse
import os
import sys


def _resolve_reference(args):
    """Merge --ref-dir (or $IRSHAPE_REF) with any explicit override flags.
    Returns a dict of resolved kwargs for reference.build_reference /
    pipeline.run_full_engine. Raises reference_dir.ReferenceDirError (caught
    in main()) if nothing usable was given."""
    from . import reference_dir as refdir_mod

    ref_dir_value = args.ref_dir or os.environ.get(refdir_mod.REF_ENV_VAR)
    ref = refdir_mod.resolve_ref_dir(ref_dir_value) if ref_dir_value else None

    gtf = args.gtf or (ref.gtf if ref else None)
    genome_build = args.genome_build or (ref.genome_build if ref else None)
    annotation_version = args.annotation_version or (ref.annotation_version if ref else None)
    fasta_path = getattr(args, "fasta", None) or (ref.fasta if ref else None)
    mappability_bw = getattr(args, "mappability_bigwig", None) or (ref.mappability_bigwig if ref else None)
    polya_bed = getattr(args, "polya_bed", None) or (ref.polya_bed if ref else None)
    cache_dir = args.reference_cache or (ref.cache_dir if ref else None)
    mappability_readlen = ref.mappability_readlen if ref else 100
    mappability_source = ref.mappability_source if ref else None

    if not gtf:
        raise refdir_mod.ReferenceDirError(
            "no GTF available: pass --gtf directly, or --ref-dir/$IRSHAPE_REF "
            "pointing at a bundle that has one under annotation/."
        )
    if not cache_dir:
        raise refdir_mod.ReferenceDirError(
            "no cache directory available: pass --reference-cache directly, or "
            "--ref-dir/$IRSHAPE_REF (caches under <ref-dir>/cache/ by default)."
        )
    if not genome_build or not annotation_version:
        raise refdir_mod.ReferenceDirError(
            "--genome-build/--annotation-version are required when not resolvable "
            "from a --ref-dir bundle's manifest.json."
        )

    if ref is not None:
        print(f"using reference bundle: {ref.describe()}", file=sys.stderr)

    return dict(
        gtf=gtf, genome_build=genome_build, annotation_version=annotation_version,
        cache_dir=cache_dir, fasta_path=fasta_path, mappability_bigwig=mappability_bw,
        mappability_readlen=mappability_readlen, mappability_source=mappability_source,
        polya_bed=polya_bed,
    )


def _cmd_build_reference(args):
    from . import reference as ref_mod

    resolved = _resolve_reference(args)
    want_genes = None
    if args.genes:
        want_genes = set(g.strip() for g in args.genes.split(",") if g.strip())

    print(f"building reference: gtf={resolved['gtf']} genome_build={resolved['genome_build']} "
          f"annotation_version={resolved['annotation_version']} "
          f"genes={'genome-wide' if want_genes is None else len(want_genes)} "
          f"context={'yes' if (resolved['fasta_path'] or resolved['mappability_bigwig']) else 'no'}",
          file=sys.stderr)
    bundle = ref_mod.build_reference(
        gtf=resolved["gtf"], genome_build=resolved["genome_build"],
        annotation_version=resolved["annotation_version"], cache_dir=resolved["cache_dir"],
        want_genes=want_genes, force=args.force, min_length=args.min_length, max_length=args.max_length,
        fasta_path=resolved["fasta_path"], mappability_bigwig=resolved["mappability_bigwig"],
        mappability_readlen=resolved["mappability_readlen"], mappability_source=resolved["mappability_source"],
    )
    path = ref_mod.cache_path(resolved["cache_dir"], resolved["genome_build"],
                               resolved["annotation_version"], bundle.scope_key)
    print(f"reference bundle: {len(bundle.introns)} introns -> {path}")


def _cmd_run(args):
    from . import schema as schema_mod

    def progress(msg):
        print(msg, file=sys.stderr, flush=True)

    resolved = _resolve_reference(args)
    want_genes = None
    if args.genes:
        want_genes = set(g.strip() for g in args.genes.split(",") if g.strip())

    if args.input == "bam":
        raise NotImplementedError(
            "the BAM adapter is stubbed in Engine Block 1 -- use --input recount3"
        )
    elif args.input != "recount3":
        raise ValueError(f"unknown --input {args.input!r}")

    import pandas as pd

    from .adapters.base import JunctionCounts
    from .adapters.recount3_adapter import Recount3Adapter

    junctions = {}
    if args.junctions:
        jdf = pd.read_csv(args.junctions, sep="\t")
        for row in jdf.itertuples(index=False):
            junctions[row.intron_id] = JunctionCounts(
                splice_left=row.splice_left, splice_right=row.splice_right,
                splice_exact=row.splice_exact,
            )
    if not args.bigwigs:
        raise ValueError("--bigwigs is required for --input recount3")
    bigwigs = [p.strip() for p in args.bigwigs.split(",") if p.strip()]
    adapter = Recount3Adapter(bigwigs, junctions)

    try:
        if args.classic_only:
            from . import reference as ref_mod
            from . import classic as classic_mod

            bundle = ref_mod.build_reference(
                gtf=resolved["gtf"], genome_build=resolved["genome_build"],
                annotation_version=resolved["annotation_version"], cache_dir=resolved["cache_dir"],
                want_genes=want_genes,
            )
            print(f"loaded reference: {len(bundle.introns)} introns", file=sys.stderr)
            table = classic_mod.compute_table(bundle, adapter, progress=progress)
        else:
            from . import pipeline as pipeline_mod

            table = pipeline_mod.run_full_engine(
                gtf=resolved["gtf"], genome_build=resolved["genome_build"],
                annotation_version=resolved["annotation_version"], cache_dir=resolved["cache_dir"],
                adapter=adapter, want_genes=want_genes, fasta_path=resolved["fasta_path"],
                mappability_bigwig=resolved["mappability_bigwig"],
                mappability_readlen=resolved["mappability_readlen"],
                mappability_source=resolved["mappability_source"], polya_bed=resolved["polya_bed"],
                progress=progress,
            )
    finally:
        adapter.close()

    table.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {len(table)} rows ({len(table) // 2} introns x 2 mask policies) -> {args.out}")
    problems = schema_mod.validate(table)
    print("contract violations:", problems if problems else "NONE", file=sys.stderr)


def _add_reference_args(p, require_genome_meta):
    p.add_argument("--ref-dir", default=None,
                    help=f"external reference bundle directory (or set $IRSHAPE_REF)")
    p.add_argument("--gtf", default=None, help="override: GTF path (else from --ref-dir)")
    p.add_argument("--genome-build", default=None, help="override: else from --ref-dir manifest")
    p.add_argument("--annotation-version", default=None, help="override: else from --ref-dir manifest")
    p.add_argument("--fasta", default=None, help="override: genome FASTA (else from --ref-dir)")
    p.add_argument("--mappability-bigwig", default=None, help="override: mappability bigWig (else from --ref-dir)")
    p.add_argument("--polya-bed", default=None, help="override: PolyASite bed.gz (else from --ref-dir)")
    p.add_argument("--reference-cache", default=None,
                    help="cache dir for built ReferenceBundles (default: <ref-dir>/cache)")
    p.add_argument("--genes", default=None, help="comma-separated gene names; omit for genome-wide")


def build_parser():
    p = argparse.ArgumentParser(prog="irshape", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    br = sub.add_parser("build-reference", help="build/cache a ReferenceBundle from a GTF/ref-dir")
    _add_reference_args(br, require_genome_meta=False)
    br.add_argument("--min-length", type=int, default=None)
    br.add_argument("--max-length", type=int, default=None)
    br.add_argument("--force", action="store_true", help="rebuild even if a cache file exists")
    br.set_defaults(func=_cmd_build_reference)

    rn = sub.add_parser("run", help="run the engine, write a schema-conformant TSV")
    _add_reference_args(rn, require_genome_meta=False)
    rn.add_argument("--input", choices=["bam", "recount3"], required=True)
    rn.add_argument("--bigwigs", default=None, help="comma-separated bigwig paths (recount3 input)")
    rn.add_argument("--junctions", default=None,
                     help="TSV with columns intron_id,splice_left,splice_right,splice_exact "
                          "(recount3 input; from parse_recount3_junction_flatfiles)")
    rn.add_argument("--classic-only", action="store_true",
                     help="Block-1 classic engine only (single unscoped bundle, no shape/GC/IPA "
                          "columns populated); default runs the full tiered engine")
    rn.add_argument("--out", required=True)
    rn.set_defaults(func=_cmd_run)

    return p


def main(argv=None):
    from . import reference_dir as refdir_mod

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except refdir_mod.ReferenceDirError as e:
        print(f"irshape: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
