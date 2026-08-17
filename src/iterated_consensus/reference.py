"""FASTA I/O and fetching a starting reference sequence from NCBI."""

from __future__ import annotations

import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path


class ReferenceError(RuntimeError):
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


def load_reference(fasta_path: Path, reference_id: str | None = None) -> FastaRecord:
    records = parse_fasta(fasta_path.read_text())
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


def write_fasta(path: Path, record_id: str, sequence: str, *, line_width: int = 70) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{record_id}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i : i + line_width] + "\n")


_NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_ncbi_accession(accession: str, out_dir: Path) -> Path:
    """Download a nucleotide FASTA record from NCBI, caching it in `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{accession}.fasta"
    if target.exists():
        return target

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
