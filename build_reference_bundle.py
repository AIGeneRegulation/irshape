#!/usr/bin/env python3
"""Engine Packaging Part 1 (SLURM CPU, LOCAL FILES ONLY, NO NETWORK) --
assemble the existing, already-built reference artifacts irshape needs into
ONE versioned, self-contained directory a user downloads separately from the
package/container ("external reference bundle" -- see irshape/reference_dir.py
and SCHEMA.md).

This script COLLECTS; it does not rebuild anything. The expensive artifacts
(genome-wide mask sweep + baked-in GC/mappability tracks, per length tier)
already exist on disk from Engine Blocks 1-3 (scripts 48/52/62); this just
hardlinks them (same filesystem, zero extra disk cost) into a clean,
versioned namespace with a manifest + checksums. Checksumming ~17GB of
reference data is CPU-bound work, hence SLURM, not the login node
(CLAUDE.md: "anything that ... reduces data" is compute).

Generic across genome/annotation builds: every source path is a flag; the
defaults below are read from two environment variables (documented in
README.md) rather than hardcoded, so pointing all of --gtf/--fasta/
--mappability-bigwig/--polya-bed/--long-cache/--awkward-mid-cache/
--genome-build/--annotation-version/--out-dir at a different assembly's
already-built artifacts assembles a new bundle the same way. It does NOT
invoke `irshape build-reference` itself -- build the per-tier caches for the
new genome/annotation first (that step legitimately needs to run: it's
exactly what "collect the existing on-disk ones" presupposes exists), then
point this script at them.

No EB shrinkage: files are linked as-is, at their full existing size.

Environment variables (see README.md "Rebuilding the reference bundle"):
  IRSHAPE_BUILD_ROOT       -- root containing results/ and data/external/
                              (default: current directory)
  IRSHAPE_BUILD_ANNOT_DIR  -- directory containing the GTF + genome FASTA
                              (default: IRSHAPE_BUILD_ROOT)
"""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone

ROOT = os.environ.get("IRSHAPE_BUILD_ROOT", ".")
IRLR = os.environ.get("IRSHAPE_BUILD_ANNOT_DIR", ROOT)
DEFAULT_OUT_DIR = f"{ROOT}/irshape-ref-GRCh38-gencodeV44"

DEFAULTS = dict(
    gtf=f"{IRLR}/gencode.v44.annotation.gtf",
    fasta=f"{IRLR}/GRCh38.primary_assembly.genome.fa",
    mappability_bigwig=f"{ROOT}/data/external/mappability/k100.Umap.MultiTrackMappability.bw",
    mappability_source=("https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/"
                         "k100.Umap.MultiTrackMappability.bw (UCSC mirror of Hoffman-lab Umap k100)"),
    mappability_readlen=100,
    polya_bed=f"{ROOT}/data/external/polyasite/atlas.clusters.2.0.GRCh38.96.bed.gz",
    polya_source="PolyASite 2.0 atlas.clusters.2.0.GRCh38.96 (https://polyasite.unibas.ch/)",
    intron_universe=f"{ROOT}/results/intron_universe.tsv",
    long_cache=(f"{ROOT}/data/processed/irshape_ref/"
                "irshape_ref__GRCh38__gencode_v44__genome_wide__minlen1000__ctx_k100.pkl"),
    awkward_mid_cache=(f"{ROOT}/data/processed/irshape_ref/"
                       "irshape_ref__GRCh38__gencode_v44__genome_wide__minlen100__maxlen1000__ctx_k100.pkl"),
    genome_build="GRCh38",
    annotation_version="gencode_v44",
    version="1.0.0",
)


def eprint(*a, **kw):
    print(*a, file=sys.stderr, flush=True, **kw)


def link_file(src, dst, mode="hardlink"):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        os.remove(dst)
    if mode == "copy":
        shutil.copy2(src, dst)
        return
    try:
        os.link(src, dst)
    except OSError as e:
        eprint(f"  hardlink failed ({e}), falling back to symlink")
        os.symlink(os.path.abspath(src), dst)


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    ap.add_argument("--gtf", default=DEFAULTS["gtf"])
    ap.add_argument("--fasta", default=DEFAULTS["fasta"])
    ap.add_argument("--mappability-bigwig", default=DEFAULTS["mappability_bigwig"])
    ap.add_argument("--mappability-source", default=DEFAULTS["mappability_source"])
    ap.add_argument("--mappability-readlen", type=int, default=DEFAULTS["mappability_readlen"])
    ap.add_argument("--polya-bed", default=DEFAULTS["polya_bed"])
    ap.add_argument("--polya-source", default=DEFAULTS["polya_source"])
    ap.add_argument("--intron-universe", default=DEFAULTS["intron_universe"])
    ap.add_argument("--long-cache", default=DEFAULTS["long_cache"],
                     help="prebuilt ReferenceBundle pickle, length_tier=='long' (min_length=1000, with GC/mappability context)")
    ap.add_argument("--awkward-mid-cache", default=DEFAULTS["awkward_mid_cache"],
                     help="prebuilt ReferenceBundle pickle, length_tier=='awkward_mid' (100<=length<1000, with context)")
    ap.add_argument("--genome-build", default=DEFAULTS["genome_build"])
    ap.add_argument("--annotation-version", default=DEFAULTS["annotation_version"])
    ap.add_argument("--version", default=DEFAULTS["version"], help="bundle version string")
    ap.add_argument("--link-mode", choices=["hardlink", "copy"], default="hardlink",
                     help="hardlink (default, zero extra disk on the same filesystem) or real copy")
    args = ap.parse_args()

    for label, path in [("--gtf", args.gtf), ("--fasta", args.fasta),
                         ("--mappability-bigwig", args.mappability_bigwig),
                         ("--polya-bed", args.polya_bed), ("--intron-universe", args.intron_universe),
                         ("--long-cache", args.long_cache), ("--awkward-mid-cache", args.awkward_mid_cache)]:
        if not os.path.exists(path):
            sys.exit(f"missing input for {label}: {path}")

    out = args.out_dir
    eprint(f"assembling bundle -> {out}")
    os.makedirs(out, exist_ok=True)

    placements = [
        (args.gtf, f"annotation/{os.path.basename(args.gtf)}",
         "GENCODE GTF annotation -- gene/transcript/exon models the engine masks introns against"),
        (args.fasta, f"genome/{os.path.basename(args.fasta)}",
         "Genome FASTA -- GC content + U12 splice-site motif lookup"),
        (args.mappability_bigwig, f"mappability/{os.path.basename(args.mappability_bigwig)}",
         f"Per-base multi-read mappability track (k={args.mappability_readlen}): {args.mappability_source}"),
        (args.polya_bed, f"polyasite/{os.path.basename(args.polya_bed)}",
         f"PolyASite atlas for IPA/alt-3' annotation corroboration: {args.polya_source}"),
        (args.long_cache, f"cache/{os.path.basename(args.long_cache)}",
         "Prebuilt ReferenceBundle: length_tier=='long' (>=1000bp) -- both mask policies "
         "(other/same/small-RNA per-base classification) + GC + mappability, baked in"),
        (args.awkward_mid_cache, f"cache/{os.path.basename(args.awkward_mid_cache)}",
         "Prebuilt ReferenceBundle: length_tier=='awkward_mid' (100-999bp) -- both mask policies "
         "+ GC + mappability, baked in"),
    ]

    # intron_universe.tsv is small and human-consumed (not read by the engine
    # at runtime -- the engine derives IDENTITY fields itself from the GTF);
    # gzip the collected copy rather than hardlinking the raw TSV.
    universe_dst = os.path.join(out, "intron_universe/intron_universe.tsv.gz")
    os.makedirs(os.path.dirname(universe_dst), exist_ok=True)
    with open(args.intron_universe, "rb") as fh_in, gzip.open(universe_dst, "wb") as fh_out:
        shutil.copyfileobj(fh_in, fh_out)
    eprint(f"  wrote {universe_dst}")

    files_manifest = []
    for src, rel, desc in placements:
        dst = os.path.join(out, rel)
        eprint(f"  linking {rel}")
        link_file(src, dst, mode=args.link_mode)

    # Hash everything actually placed under `out` (includes the gzipped
    # universe copy), after all linking/copying is done.
    eprint("computing checksums...")
    all_rel = []
    for dirpath, _dirnames, filenames in os.walk(out):
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            all_rel.append(os.path.relpath(full, out))
    all_rel.sort()

    checksums = {}
    for rel in all_rel:
        full = os.path.join(out, rel)
        size = os.path.getsize(full)
        digest = sha256_of(full)
        checksums[rel] = (digest, size)
        eprint(f"  sha256 {rel}: {digest} ({size} bytes)")

    desc_by_rel = {rel: desc for _src, rel, desc in placements}
    desc_by_rel["intron_universe/intron_universe.tsv.gz"] = (
        "Flat identity table (intron_id/gene/chrom/coords/strand) for reference/inspection; "
        "not read by the engine at runtime"
    )

    build_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    manifest = {
        "name": "irshape-ref-GRCh38-gencodeV44",
        "version": args.version,
        "genome_build": args.genome_build,
        "annotation_version": args.annotation_version,
        "build_date": build_date,
        "source_project": ROOT,
        "mappability": {"source": args.mappability_source, "readlen": args.mappability_readlen},
        "polyasite": {"source": args.polya_source},
        "tiering": {
            "long": {"min_length": 1000, "max_length": None, "cache": os.path.basename(args.long_cache)},
            "awkward_mid": {"min_length": 100, "max_length": 1000, "cache": os.path.basename(args.awkward_mid_cache)},
            "sub_read": {"min_length": None, "max_length": 100, "cache": None,
                         "note": "no prebuilt cache -- sub_read's classic-ratio-only path needs no "
                                 "GC/mappability context, so irshape builds it on demand in seconds "
                                 "from annotation/*.gtf on first use"},
        },
        "files": [
            {"path": rel, "bytes": size, "sha256": digest, "description": desc_by_rel.get(rel, "")}
            for rel, (digest, size) in checksums.items()
        ],
    }

    manifest_path = os.path.join(out, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    eprint(f"wrote {manifest_path}")

    checksums_path = os.path.join(out, "CHECKSUMS.sha256")
    with open(checksums_path, "w") as fh:
        for rel in all_rel:
            digest, _size = checksums[rel]
            fh.write(f"{digest}  {rel}\n")
    eprint(f"wrote {checksums_path}")

    total_bytes = sum(size for _d, size in checksums.values())
    md_path = os.path.join(out, "MANIFEST.md")
    with open(md_path, "w") as fh:
        fh.write(f"# irshape-ref-GRCh38-gencodeV44 v{args.version}\n\n")
        fh.write(f"Built {build_date} from `{ROOT}`.\n\n")
        fh.write(f"genome_build: `{args.genome_build}` / annotation_version: `{args.annotation_version}`\n\n")
        fh.write(f"Total size: {total_bytes / 1e9:.2f} GB across {len(all_rel)} files "
                 "(sizes are full, un-shrunk originals -- see manifest.json).\n\n")
        fh.write("| file | bytes | sha256 | description |\n|---|---|---|---|\n")
        for rel in all_rel:
            digest, size = checksums[rel]
            fh.write(f"| `{rel}` | {size} | `{digest[:12]}...` | {desc_by_rel.get(rel, '')} |\n")
        fh.write("\n## Verify\n\n```\ncd irshape-ref-GRCh38-gencodeV44 && sha256sum -c CHECKSUMS.sha256\n```\n")
    eprint(f"wrote {md_path}")

    eprint(f"\nbundle complete: {out}")
    eprint(f"total: {total_bytes / 1e9:.2f} GB, {len(all_rel)} files")


if __name__ == "__main__":
    main()
