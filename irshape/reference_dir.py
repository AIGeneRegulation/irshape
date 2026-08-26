"""External reference bundle resolution (Engine Packaging Part 2).

The engine's reference materials (GTF, genome FASTA, mappability bigWig,
PolyASite atlas, plus the expensive-to-build per-tier mask/GC/mappability
caches) are NOT bundled into the irshape package or its container image --
they are a separately-downloaded, versioned directory
(`irshape-ref-GRCh38-gencodeV44/`, see `build_reference_bundle.py`).
This module is the one place that knows that directory's layout, so `cli.py`
just asks a `ReferenceDir` for paths instead of hardcoding them.

Resolution order: `--ref-dir` CLI flag, then the `IRSHAPE_REF` environment
variable. Neither set (or the path not a valid bundle) is a hard error with
a message telling the user where to get one -- there is no bundled fallback.
"""
from __future__ import annotations

import json
import os

REF_ENV_VAR = "IRSHAPE_REF"

REF_BUNDLE_URL_PLACEHOLDER = "<URL/DOI placeholder -- see README.md>"


class ReferenceDirError(RuntimeError):
    pass


def _no_ref_dir_message() -> str:
    return (
        "no irshape reference bundle given.\n"
        "irshape needs an external reference bundle (GTF + genome FASTA + "
        "mappability track + PolyASite atlas + prebuilt mask/GC caches) that is "
        "downloaded separately from the package/container:\n"
        f"  1) download the bundle:  {REF_BUNDLE_URL_PLACEHOLDER}\n"
        "  2) extract it, then either:\n"
        "       irshape run --ref-dir /path/to/irshape-ref-GRCh38-gencodeV44 ...\n"
        f"     or export {REF_ENV_VAR}=/path/to/irshape-ref-GRCh38-gencodeV44"
    )


def resolve_ref_dir(cli_value: str | None):
    """cli_value: the --ref-dir argument, or None. Returns a ReferenceDir, or
    raises ReferenceDirError with a message telling the user how to fix it."""
    path = cli_value or os.environ.get(REF_ENV_VAR)
    if not path:
        raise ReferenceDirError(_no_ref_dir_message())
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        raise ReferenceDirError(
            f"reference directory not found: {path!r}\n{_no_ref_dir_message()}"
        )
    manifest_path = os.path.join(path, "manifest.json")
    if not os.path.exists(manifest_path):
        raise ReferenceDirError(
            f"{path!r} does not look like an irshape reference bundle "
            f"(no manifest.json in it).\n{_no_ref_dir_message()}"
        )
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    return ReferenceDir(path, manifest)


class ReferenceDir:
    """Resolved paths inside one irshape-ref-<genome>-<annotation>/ bundle.
    Every property returns None (not an error) if that particular material
    is absent -- e.g. a minimal bundle with no PolyASite atlas can still run
    the classic engine; callers decide what's mandatory for what they asked
    the tool to do."""

    def __init__(self, path: str, manifest: dict):
        self.path = path
        self.manifest = manifest

    def _find_one(self, subdir: str, suffixes: tuple[str, ...]):
        d = os.path.join(self.path, subdir)
        if not os.path.isdir(d):
            return None
        for name in sorted(os.listdir(d)):
            if name.endswith(suffixes):
                return os.path.join(d, name)
        return None

    @property
    def gtf(self):
        return self._find_one("annotation", (".gtf",))

    @property
    def fasta(self):
        return self._find_one("genome", (".fa", ".fasta"))

    @property
    def mappability_bigwig(self):
        return self._find_one("mappability", (".bw", ".bigWig"))

    @property
    def polya_bed(self):
        return self._find_one("polyasite", (".bed.gz", ".bed"))

    @property
    def cache_dir(self):
        return os.path.join(self.path, "cache")

    @property
    def genome_build(self):
        return self.manifest.get("genome_build", "GRCh38")

    @property
    def annotation_version(self):
        return self.manifest.get("annotation_version", "gencode_v44")

    @property
    def mappability_readlen(self):
        return self.manifest.get("mappability", {}).get("readlen", 100)

    @property
    def mappability_source(self):
        return self.manifest.get("mappability", {}).get("source")

    def describe(self) -> str:
        name = self.manifest.get("name", os.path.basename(self.path))
        version = self.manifest.get("version", "?")
        return f"{name} v{version} ({self.path})"
