"""Running the user-configured consensus pipeline and reading its result."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from .commands import CommandError, RenderedCommand, render_command, run_command
from .config import ConsensusSpec
from .errors import IteratedConsensusError
from .metrics import base_composition
from .reference import parse_fasta
from .templating import CatResolver, render


class ConsensusError(IteratedConsensusError, RuntimeError):
    """Raised when the consensus pipeline fails or produces an unusable result."""


@dataclass(frozen=True)
class ConsensusResult:
    fasta_path: Path
    record_id: str
    sequence: str
    length: int
    composition: dict[str, int]


def render_consensus_steps(
    spec: ConsensusSpec,
    values: Mapping[str, object],
    *,
    cat_resolver: CatResolver | None = None,
) -> list[RenderedCommand]:
    return [render_command(step, values, cat_resolver=cat_resolver) for step in spec.steps]


def run_consensus(
    spec: ConsensusSpec,
    values: Mapping[str, object],
    *,
    log_dir: Path | None = None,
    cat_resolver: CatResolver | None = None,
) -> ConsensusResult:
    for i, step in enumerate(spec.steps):
        rendered = render_command(step, values, cat_resolver=cat_resolver)
        log_path = None if log_dir is None else log_dir / f"consensus_step_{i:02d}.log"
        try:
            run_command(rendered, log_path=log_path)
        except CommandError as exc:
            raise ConsensusError(f"consensus step {i} failed: {exc}") from exc

    output_path = Path(render(spec.output, values, cat_resolver=cat_resolver))
    if not output_path.exists():
        raise ConsensusError(
            f"consensus steps finished but expected output '{output_path}' was not created"
        )

    try:
        records = parse_fasta(output_path.read_text())
    except OSError as exc:
        raise ConsensusError(f"could not read consensus output '{output_path}': {exc}") from exc

    if len(records) != 1:
        found = [r.id for r in records]
        raise ConsensusError(
            f"expected exactly one sequence in consensus output '{output_path}', "
            f"found {len(records)}: {found}"
        )

    record = records[0]
    if not record.sequence:
        raise ConsensusError(f"consensus output '{output_path}' has an empty sequence")

    return ConsensusResult(
        fasta_path=output_path,
        record_id=record.id,
        sequence=record.sequence,
        length=len(record.sequence),
        composition=base_composition(record.sequence),
    )
