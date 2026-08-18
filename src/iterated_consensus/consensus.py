"""Running the user-configured consensus pipeline and reading its result."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from .commands import (
    CommandError,
    CommandRun,
    RenderedCommand,
    render_command,
    run_command,
    run_logged_command,
    run_tool_version_command,
)
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
    tool_versions: dict[str, str] = field(default_factory=dict)
    commands: tuple[CommandRun, ...] = ()


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
    tool_versions: dict[str, str] = {}
    for name, step in spec.tool_versions.items():
        try:
            tool_versions[name] = run_tool_version_command(step, values)
        except CommandError as exc:
            raise ConsensusError(f"[consensus] tool-versions '{name}' failed: {exc}") from exc

    commands: list[CommandRun] = []
    for i, step in enumerate(spec.steps):
        rendered = render_command(step, values, cat_resolver=cat_resolver)
        name = f"consensus_step_{i:02d}"
        try:
            if log_dir is None:
                run_command(rendered)
            else:
                commands.append(run_logged_command(name, rendered, log_dir))
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
        tool_versions=tool_versions,
        commands=tuple(commands),
    )
