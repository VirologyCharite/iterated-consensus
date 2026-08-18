"""FASTA I/O and fetching a starting reference sequence from NCBI."""

from __future__ import annotations

import os
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .errors import IteratedConsensusError


class ReferenceError(IteratedConsensusError, RuntimeError):
    """Raised for missing/ambiguous/malformed reference sequences."""


@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str


def parse_fasta(text: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    description = ""
    seq_lines: list[str] = []

    def flush() -> None:
        if header is not None:
            records.append(
                FastaRecord(id=header, description=description, sequence="".join(seq_lines))
            )

    for line in text.splitlines():
        if line.startswith(">"):
            flush()
            parts = line[1:].strip().split(None, 1)
            header = parts[0] if parts else ""
            description = parts[1] if len(parts) > 1 else ""
            seq_lines = []
        else:
            stripped = line.strip()
            if stripped:
                seq_lines.append(stripped)
    flush()
    return records


def _pick_record(
    records: list[FastaRecord], fasta_path: Path, reference_id: str | None
) -> FastaRecord:
    if not records:
        raise ReferenceError(f"{fasta_path} has no sequences")
    if reference_id is not None:
        for record in records:
            if record.id == reference_id:
                return record
        ids = [r.id for r in records]
        raise ReferenceError(
            f"reference_id '{reference_id}' not found in {fasta_path}; available: {ids}"
        )
    if len(records) == 1:
        return records[0]
    ids = [r.id for r in records]
    raise ReferenceError(f"{fasta_path} has multiple sequences {ids}; specify reference_id")


def load_reference(fasta_path: Path, reference_id: str | None = None) -> FastaRecord:
    records = parse_fasta(fasta_path.read_text())
    return _pick_record(records, fasta_path, reference_id)


def write_fasta(path: Path, record_id: str, sequence: str, *, line_width: int = 70) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{record_id}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i : i + line_width] + "\n")


_NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# RefSeq (NC_045512.2, NZ_CP012345.1) or GenBank (MN908947.3, U12345.1) style:
# 1-2 letters, optional underscore + up to 2 more letters (RefSeq assembly
# accessions like NZ_ embed a GenBank-style accession after the prefix),
# 5-9 digits, optional dot-version.
_ACCESSION_RE = re.compile(r"^[A-Za-z]{1,2}_?[A-Za-z]{0,2}[0-9]{5,9}(\.[0-9]+)?$")


def looks_like_ncbi_accession(name: str) -> bool:
    """Heuristic: could `name` plausibly be a RefSeq/GenBank nucleotide accession?

    Deliberately permissive -- a false positive just means an attempted
    fetch that fails cleanly with a normal ReferenceError; a false negative
    silently loses a reference the caller could have auto-fetched.
    """
    return bool(_ACCESSION_RE.match(name))


def ncbi_cache_path(accession: str, cache_dir: Path) -> Path:
    """Where `fetch_ncbi_accession` would cache (or already has cached) `accession`."""
    return cache_dir / f"{accession}.fasta"


def fetch_ncbi_accession(accession: str, out_dir: Path) -> Path:
    """Download a nucleotide FASTA record from NCBI, caching it in `out_dir`."""
    target = ncbi_cache_path(accession, out_dir)
    if target.exists():
        return target
    out_dir.mkdir(parents=True, exist_ok=True)

    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"}
    )
    url = f"{_NCBI_EFETCH_URL}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            text = response.read().decode()
    except OSError as exc:
        raise ReferenceError(f"failed to fetch accession '{accession}' from NCBI: {exc}") from exc

    if not text.startswith(">"):
        raise ReferenceError(
            f"unexpected response fetching accession '{accession}' from NCBI: {text[:200]!r}"
        )
    target.write_text(text)
    return target


@dataclass(frozen=True)
class ResolvedReference:
    record: FastaRecord
    source_path: Path
    is_whole_file: bool
    """True if `source_path` contains exactly `record` and nothing else --
    i.e. it's safe to symlink a destination straight to `source_path`
    instead of writing out a fresh copy of `record`'s sequence."""


def would_resolve_reference(*, reference_id: str | None, reference_fasta: Path | None) -> bool:
    """Whether `resolve_reference` would find something, without doing its
    (possibly network-fetching) work -- for previewing a dry run without
    actually hitting NCBI."""
    if reference_fasta is not None:
        return True
    return reference_id is not None and looks_like_ncbi_accession(reference_id)


def resolve_reference(
    *, reference_id: str | None, reference_fasta: Path | None, cache_dir: Path
) -> ResolvedReference | None:
    """Resolve "the reference sequence" from `reference_fasta` and/or `reference_id`.

    `reference_id` is just a name -- it can match a record in
    `reference_fasta` (if given), or, if `reference_fasta` isn't given and
    the name looks like an NCBI accession, be fetched.

    Returns None if nothing could be resolved from what's given here --
    callers decide whether that's acceptable (a BAM-start run can tolerate
    having no reference for iteration 0; a FASTQ-start run can't).
    """
    if not would_resolve_reference(reference_id=reference_id, reference_fasta=reference_fasta):
        return None
    if reference_fasta is not None:
        fasta_path = reference_fasta
    else:
        assert reference_id is not None  # guaranteed by would_resolve_reference above
        fasta_path = fetch_ncbi_accession(reference_id, cache_dir)

    records = parse_fasta(fasta_path.read_text())
    record = _pick_record(records, fasta_path, reference_id if reference_fasta is not None else None)
    return ResolvedReference(record=record, source_path=fasta_path, is_whole_file=len(records) == 1)


def symlink_reference(link_path: Path, target_path: Path) -> None:
    """Make `link_path` a *relative* symlink to `target_path`.

    A no-op if `link_path` already exists (as a file or a symlink, even a
    broken one) -- e.g. from an earlier call on a resumed run.
    """
    if link_path.exists() or link_path.is_symlink():
        return
    link_path.parent.mkdir(parents=True, exist_ok=True)
    relative_target = os.path.relpath(target_path.resolve(), start=link_path.parent.resolve())
    link_path.symlink_to(relative_target)
