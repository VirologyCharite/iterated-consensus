"""BAM inspection and read extraction, via the samtools functionality bundled in pysam.

No external `samtools` binary is required -- everything here goes through
pysam's htslib bindings, either the high-level AlignmentFile API or the
samtools-command wrappers (`pysam.view`, `pysam.index`, `pysam.merge`,
`pysam.fastq`).
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pysam

from .errors import IteratedConsensusError


class BamError(IteratedConsensusError, RuntimeError):
    """Raised for BAM-related failures: missing/ambiguous reference, bad input, etc."""


def list_references(bam_path: Path) -> list[str]:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return list(bam.references)


def resolve_reference_name(bam_path: Path, reference_id: str | None) -> str:
    """Pick the single reference to call a consensus from.

    If the BAM has exactly one reference, it is used automatically.
    Otherwise `reference_id` must be given and must match one of them.
    """
    refs = list_references(bam_path)
    if not refs:
        raise BamError(f"{bam_path} has no references in its header")
    if reference_id is not None:
        if reference_id not in refs:
            raise BamError(f"reference '{reference_id}' not found in {bam_path}; available: {refs}")
        return reference_id
    if len(refs) == 1:
        return refs[0]
    raise BamError(f"{bam_path} has multiple references {refs}; specify reference_id to pick one")


def _is_coordinate_sorted(bam_path: Path) -> bool:
    """Check the header's own SO tag -- just reads the header, not the reads.

    A file that can't even be opened as a BAM (corrupt, not actually a BAM,
    ...) isn't verifiably sorted either -- treated as False here, so the
    caller falls through to an actual sort/index attempt, which will raise
    a proper BamError explaining the real problem.
    """
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            return bam.header.get("HD", {}).get("SO") == "coordinate"
    except (ValueError, OSError):
        return False


def _sort_to(source: Path, dest: Path) -> None:
    try:
        pysam.sort("-o", str(dest), str(source))
    except pysam.SamtoolsError as exc:
        raise BamError(f"could not sort {source} for indexing: {exc}") from exc


def _sort_in_place(bam_path: Path) -> None:
    sorting_path = bam_path.with_name(bam_path.name + ".sorting.tmp")
    _sort_to(bam_path, sorting_path)
    sorting_path.replace(bam_path)


def _index(bam_path: Path) -> None:
    try:
        pysam.index(str(bam_path))
    except pysam.SamtoolsError as exc:
        raise BamError(f"could not index {bam_path}: {exc}") from exc


def ensure_indexed(bam_path: Path) -> None:
    """Create a `.bam.bai` index for `bam_path`, unless one (or a `.csi`) already exists.

    Sorts `bam_path` in place first if it isn't already coordinate-sorted --
    indexing requires that, and a mapper's own output isn't always sorted.
    The header's own sort-order tag decides whether to sort, rather than
    just trying to index and reacting to failure: for a large, genuinely
    already-sorted BAM that only costs a header read, not a failed index
    attempt (which would also print samtools' own error log to stderr even
    though nothing is actually wrong).
    """
    if Path(str(bam_path) + ".bai").exists() or Path(str(bam_path) + ".csi").exists():
        return

    if not _is_coordinate_sorted(bam_path):
        _sort_in_place(bam_path)
        _index(bam_path)
        return

    try:
        _index(bam_path)
    except BamError:
        # Header claimed coordinate-sorted but indexing failed anyway --
        # sort for real and retry.
        _sort_in_place(bam_path)
        _index(bam_path)


def ensure_indexed_readonly(bam_path: Path, work_dir: Path) -> Path:
    """Like `ensure_indexed`, but for a BAM this process doesn't own (e.g. a
    user's own BAM-start input) -- never rewrites `bam_path` itself.

    An index is still created directly beside `bam_path` if it's already
    coordinate-sorted and missing one: that's purely additive (doesn't touch
    the BAM's own bytes) and is the conventional place tools look for it
    anyway. But making an unsorted BAM indexable at all means reordering its
    contents, which would mean rewriting `bam_path` -- so in that case (or
    if indexing beside it fails for any other reason, e.g. a read-only
    source directory), a sorted copy is written under `work_dir` instead,
    indexed there, and returned. Either way, `bam_path` itself is never
    modified.
    """
    if Path(str(bam_path) + ".bai").exists() or Path(str(bam_path) + ".csi").exists():
        return bam_path

    if _is_coordinate_sorted(bam_path):
        try:
            _index(bam_path)
            return bam_path
        except BamError:
            pass  # couldn't index beside it -- fall through to a copy below

    work_dir.mkdir(parents=True, exist_ok=True)
    sorted_path = work_dir / f"{bam_path.stem}.sorted.bam"
    if not sorted_path.exists():
        _sort_to(bam_path, sorted_path)
    ensure_indexed(sorted_path)  # this copy is ours -- in-place is fine
    return sorted_path


def count_mapped_reads(bam_path: Path, *, reference_name: str | None = None) -> int:
    """Count primary mapped reads, optionally restricted to one reference."""
    if reference_name is not None:
        ensure_indexed(bam_path)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        iterator = bam.fetch(until_eof=True) if reference_name is None else bam.fetch(reference_name)
        return sum(
            1 for read in iterator if not read.is_unmapped and not read.is_secondary
            and not read.is_supplementary
        )


def get_reference_length(bam_path: Path, reference_name: str) -> int:
    """The length (LN) the BAM header records for one of its references."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return bam.get_reference_length(reference_name)


def _fastq_has_reads(path: Path) -> bool:
    """Whether a samtools-fastq gzip output actually contains any reads.

    Raw file size isn't useful here: an "empty" category still comes out as
    a valid (small) gzip stream with nonzero byte size.
    """
    if not path.exists():
        return False
    with gzip.open(path, "rb") as f:
        return f.read(1) != b""


def _collect_fastq_result(mate1: Path, mate2: Path, unpaired: Path) -> dict[str, Path]:
    result: dict[str, Path] = {}
    if _fastq_has_reads(mate1):
        result["mate1"] = mate1
        result["mate2"] = mate2
    if _fastq_has_reads(unpaired):
        result["unpaired"] = unpaired
    return result


def extract_fastq(
    bam_path: Path,
    *,
    reference_name: str,
    mode: str,
    out_dir: Path,
) -> dict[str, Path]:
    """Extract reads from a BAM into fastq.gz files, split by pairing.

    `mode` is one of:
      - "ref": reads aligned to `reference_name` only
      - "ref+unal": that, plus reads that didn't map anywhere
      - "all": every read in the BAM

    Returns a dict with whichever of "mate1"/"mate2"/"unpaired" keys ended up
    non-empty (empty categories are omitted).

    Idempotent: if `out_dir` already has extraction output from an earlier
    call (e.g. resuming a run), that's reused instead of re-extracting.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    mate1 = out_dir / "mate1.fastq.gz"
    mate2 = out_dir / "mate2.fastq.gz"
    unpaired = out_dir / "unpaired.fastq.gz"
    existing = _collect_fastq_result(mate1, mate2, unpaired)
    if existing:
        return existing

    if mode == "all":
        source_bam = bam_path
    elif mode == "ref":
        ensure_indexed(bam_path)
        source_bam = out_dir / "extraction_source.bam"
        pysam.view("-b", "-o", str(source_bam), str(bam_path), reference_name, catch_stdout=False)
    elif mode == "ref+unal":
        ensure_indexed(bam_path)
        mapped_bam = out_dir / "extraction_mapped.bam"
        unmapped_bam = out_dir / "extraction_unmapped.bam"
        pysam.view("-b", "-o", str(mapped_bam), str(bam_path), reference_name, catch_stdout=False)
        pysam.view("-b", "-f", "4", "-o", str(unmapped_bam), str(bam_path), catch_stdout=False)
        source_bam = out_dir / "extraction_source.bam"
        pysam.merge("-f", str(source_bam), str(mapped_bam), str(unmapped_bam))
    else:
        raise BamError(f"unknown bam_reads mode {mode!r}; must be 'ref', 'ref+unal', or 'all'")

    pysam.fastq(
        "-1", str(mate1),
        "-2", str(mate2),
        "-0", str(unpaired),
        str(source_bam),
    )

    result = _collect_fastq_result(mate1, mate2, unpaired)
    if not result:
        raise BamError(f"no reads extracted from {bam_path} (mode={mode!r})")
    return result
