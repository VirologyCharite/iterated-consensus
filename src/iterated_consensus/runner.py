"""Orchestrates the iterate-until-convergence loop.

Iterations are numbered from 0. iter_000 is the bootstrap step: it produces
the first consensus from whatever mapping was already available at the
start, rather than one this tool built by iterating. For a FASTQ-start run
that means building an index from the given reference and mapping into it
(the same shape as every later iteration); for a BAM-start run it means
calling a consensus directly from the input BAM, with no mapping step at
all. From iter_001 on, every iteration has the same shape regardless of
start mode: build an index from the previous iteration's consensus.fasta,
remap, call a new consensus.

Directory layout written under `out_dir`:

    reads/                     extracted/cached read files, built once
    reference_initial.fasta    normalized starting reference (FASTQ-start only)
    iter_NNN/
        <mapper>_index.*       index files (absent for iter_000 of a BAM-start run)
        <mapper>.bam           that mapper's mapping output (same caveat)
        merged.bam             (only if >1 mapper) all mapper BAMs merged
        consensus.fasta        this iteration's consensus
        stats.json             this iteration's metrics
        logs/                  captured stdout+stderr of every command run
    metrics.tsv                one row per iteration
    summary.json                final outcome: iterations run, converged?, timing
    index.html                  human-readable report rendered from summary.json

`run()` resumes automatically: if `out_dir` already has a summary.json from
an earlier call, it continues from the next iteration (see
`_load_resume_state`) instead of starting over -- e.g. after raising
max_iterations in the config because a previous run stopped without
converging. index.html is (re)written every time `run()` returns, via
`report.write_report`.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path

import pysam

from .bam import (
    count_mapped_reads,
    ensure_indexed,
    ensure_indexed_readonly,
    extract_fastq,
    get_reference_length,
    resolve_reference_name,
)
from .commands import RenderedCommand, render_command, run_command
from .config import Config, ConsensusSpec, InputSpec
from .consensus import ConsensusResult, run_consensus
from .errors import IteratedConsensusError
from .metrics import (
    ConvergenceState,
    base_composition,
    check_convergence,
    sequence_identity,
)
from .reads import ReadsCatCache
from .reference import FastaRecord, parse_fasta, resolve_reference, write_fasta
from .report import write_report
from .templating import ReadsList


class RunnerError(IteratedConsensusError, RuntimeError):
    """Raised for run-level setup failures (bad/missing input, etc.)."""


@dataclass(frozen=True)
class IterationRecord:
    iteration: int
    reads_mapped: int
    consensus_length: int
    identity_to_previous: float | None
    elapsed_seconds: float


@dataclass(frozen=True)
class RunResult:
    output_dir: Path
    iterations: tuple[IterationRecord, ...]
    converged: bool
    total_elapsed_seconds: float


@dataclass(frozen=True)
class DryRunPreview:
    lines: tuple[str, ...]
    note: str


def _reads_values_from_fastq(input_spec: InputSpec) -> dict[str, object]:
    values: dict[str, object] = {}
    if input_spec.mate1:
        values["reads_1"] = ReadsList("reads_1", input_spec.mate1)
        values["reads_2"] = ReadsList("reads_2", input_spec.mate2)
    if input_spec.unpaired:
        values["reads_single"] = ReadsList("reads_single", input_spec.unpaired)
    return values


def _reads_values_from_extraction(extracted: dict[str, Path]) -> dict[str, object]:
    values: dict[str, object] = {}
    if "mate1" in extracted:
        values["reads_1"] = ReadsList("reads_1", (extracted["mate1"],))
        values["reads_2"] = ReadsList("reads_2", (extracted["mate2"],))
    if "unpaired" in extracted:
        values["reads_single"] = ReadsList("reads_single", (extracted["unpaired"],))
    return values


@dataclass(frozen=True)
class _InitialState:
    reads_values: dict[str, object]
    starting_bam: Path | None  # set only for a BAM-origin run
    reference_path: Path | None  # may still be None for a BAM-start run -- see _prepare_initial_state
    iteration_0_needs_mapping: bool


def _check_reference_matches_bam(record: FastaRecord, bam_path: Path, reference_name: str) -> None:
    if record.id != reference_name:
        raise RunnerError(
            f"reference sequence id '{record.id}' does not match BAM reference '{reference_name}'"
        )
    bam_length = get_reference_length(bam_path, reference_name)
    if len(record.sequence) != bam_length:
        raise RunnerError(
            f"reference '{reference_name}' is {len(record.sequence)} bp, but the BAM header says "
            f"'{reference_name}' is {bam_length} bp -- these must match exactly"
        )


def _consensus_uses_reference_placeholder(consensus: ConsensusSpec) -> bool:
    for step in consensus.steps:
        texts = step if isinstance(step, list) else (step,)
        if any("{reference}" in text for text in texts):
            return True
    return False


def _prepare_initial_state(config: Config, out_dir: Path) -> _InitialState:
    input_spec = config.input
    if input_spec is None:
        raise RunnerError("no input specified: config needs an [input] section")
    input_spec.validate()
    reads_dir = out_dir / "reads"
    cache_dir = out_dir / "reference_cache"

    if input_spec.bam is not None:
        reference_name = resolve_reference_name(input_spec.bam, input_spec.reference_id)
        record = resolve_reference(
            reference_id=reference_name,
            reference_fasta=input_spec.reference_fasta,
            cache_dir=cache_dir,
        )
        reference_path: Path | None = None
        if record is not None:
            _check_reference_matches_bam(record, input_spec.bam, reference_name)
            reference_path = out_dir / "reference_initial.fasta"
            write_fasta(reference_path, record.id, record.sequence)
        elif _consensus_uses_reference_placeholder(config.consensus):
            raise RunnerError(
                "no reference is available for iteration 0 (BAM-start with no "
                "[input].reference_fasta given, and contig "
                f"'{reference_name}' doesn't look like an NCBI accession), but [consensus] steps "
                "use {reference} -- this can never succeed, since iteration 0 always runs first "
                "and has no previously-computed consensus to fall back on. Give "
                "[input].reference_fasta, or set [input].reference_id to an NCBI accession."
            )

        # Never rewrite the user's own input BAM: region queries below need
        # it indexed (and, if it isn't already sorted, sorted) -- do that
        # against a copy under reads_dir instead, if needed.
        safe_bam = ensure_indexed_readonly(input_spec.bam, reads_dir)

        extracted = extract_fastq(
            safe_bam,
            reference_name=reference_name,
            mode=input_spec.bam_reads,
            out_dir=reads_dir,
        )
        starting_bam = reads_dir / "iteration_0_source.bam"
        if not starting_bam.exists():
            pysam.view(
                "-b", "-o", str(starting_bam), str(safe_bam), reference_name, catch_stdout=False
            )
        return _InitialState(
            reads_values=_reads_values_from_extraction(extracted),
            starting_bam=starting_bam,
            reference_path=reference_path,
            iteration_0_needs_mapping=False,
        )

    record = resolve_reference(
        reference_id=input_spec.reference_id,
        reference_fasta=input_spec.reference_fasta,
        cache_dir=cache_dir,
    )
    if record is None:
        raise RunnerError(
            "FASTQ-start needs a starting reference: give [input].reference_fasta, or set "
            "[input].reference_id to an NCBI accession"
        )
    reference_path = out_dir / "reference_initial.fasta"
    write_fasta(reference_path, record.id, record.sequence)
    return _InitialState(
        reads_values=_reads_values_from_fastq(input_spec),
        starting_bam=None,
        reference_path=reference_path,
        iteration_0_needs_mapping=True,
    )


@dataclass(frozen=True)
class _MappingStep:
    mapper_name: str
    index_rendered: RenderedCommand
    map_rendered: RenderedCommand
    bam_path: Path


def _plan_mapping(
    config: Config,
    iter_dir: Path,
    base_values: dict[str, object],
    cat_cache: ReadsCatCache,
) -> tuple[list[_MappingStep], Path]:
    """Render (but do not run) every configured mapper's index/map commands.

    Returns (one _MappingStep per mapper, the bam path the consensus step
    would end up reading from -- a single mapper's bam, or merged.bam).
    """
    steps: list[_MappingStep] = []
    for mapper in config.mappers:
        index_prefix = iter_dir / f"{mapper.name}_index"
        values = {**base_values, "index_prefix": str(index_prefix)}
        index_rendered = render_command(mapper.index_cmd, values, cat_resolver=cat_cache.resolve)

        bam_out = iter_dir / f"{mapper.name}.bam"
        map_values = {**values, "bam": str(bam_out)}
        map_rendered = render_command(mapper.map_cmd, map_values, cat_resolver=cat_cache.resolve)

        steps.append(_MappingStep(mapper.name, index_rendered, map_rendered, bam_out))

    bam_path = steps[0].bam_path if len(steps) == 1 else iter_dir / "merged.bam"
    return steps, bam_path


def _run_mapping(
    config: Config,
    iter_dir: Path,
    base_values: dict[str, object],
    cat_cache: ReadsCatCache,
    *,
    log: bool,
) -> Path:
    steps, bam_path = _plan_mapping(config, iter_dir, base_values, cat_cache)
    for step in steps:
        run_command(
            step.index_rendered,
            log_path=(iter_dir / "logs" / f"{step.mapper_name}_index.log") if log else None,
        )
        run_command(
            step.map_rendered,
            log_path=(iter_dir / "logs" / f"{step.mapper_name}_map.log") if log else None,
        )

    if len(steps) > 1:
        pysam.merge("-f", str(bam_path), *(str(s.bam_path) for s in steps))

    # Mapper commands don't always leave a .bam.bai behind; several
    # consensus tools (e.g. `samtools consensus`) require one.
    ensure_indexed(bam_path)
    return bam_path


def _write_iteration_stats(
    iter_dir: Path,
    consensus_result: ConsensusResult,
    reads_mapped: int,
    identity: float | None,
    elapsed: float,
) -> None:
    stats = {
        "reads_mapped": reads_mapped,
        "consensus_length": consensus_result.length,
        "identity_to_previous": identity,
        "elapsed_seconds": elapsed,
        "composition": base_composition(consensus_result.sequence),
    }
    (iter_dir / "stats.json").write_text(json.dumps(stats, indent=2))


def _write_metrics_tsv(out_dir: Path, records: list[IterationRecord]) -> None:
    lines = ["iteration\treads_mapped\tconsensus_length\tidentity_to_previous\telapsed_seconds"]
    for r in records:
        identity = "" if r.identity_to_previous is None else f"{r.identity_to_previous:.4f}"
        lines.append(
            f"{r.iteration}\t{r.reads_mapped}\t{r.consensus_length}\t{identity}\t"
            f"{r.elapsed_seconds:.3f}"
        )
    (out_dir / "metrics.tsv").write_text("\n".join(lines) + "\n")


def _write_summary_json(
    out_dir: Path, records: list[IterationRecord], converged: bool, total_elapsed: float
) -> None:
    summary = {
        "iterations_run": len(records),
        "converged": converged,
        "total_elapsed_seconds": total_elapsed,
        "iterations": [asdict(r) for r in records],
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))


@dataclass(frozen=True)
class _ResumeState:
    next_iteration: int
    last_consensus_path: Path
    last_consensus_sequence: str
    convergence_state: ConvergenceState
    converged: bool
    records: list[IterationRecord]
    total_elapsed_seconds: float


def _load_resume_state(out_dir: Path, config: Config) -> _ResumeState | None:
    """Reconstruct where a previous run under `out_dir` left off, if any.

    Resume is detected purely from a prior `summary.json` -- there's no
    separate flag. A missing/absent summary.json means a fresh run. Replays
    the historical identity values through `check_convergence` using the
    *current* config's thresholds, so a loosened/tightened convergence
    setting takes effect on resume too.
    """
    summary_path = out_dir / "summary.json"
    if not summary_path.exists():
        return None

    try:
        summary = json.loads(summary_path.read_text())
        old_records = [IterationRecord(**r) for r in summary["iterations"]]
    except (json.JSONDecodeError, KeyError, TypeError) as exc:
        raise RunnerError(f"cannot resume: {summary_path} is not a valid run summary: {exc}") from exc
    if not old_records:
        return None

    last = old_records[-1]
    consensus_path = out_dir / f"iter_{last.iteration:03d}" / "consensus.fasta"
    if not consensus_path.exists():
        raise RunnerError(f"cannot resume: {consensus_path} referenced by {summary_path} is missing")
    fasta_records = parse_fasta(consensus_path.read_text())
    if len(fasta_records) != 1:
        raise RunnerError(f"cannot resume: {consensus_path} does not contain exactly one sequence")

    convergence_state = ConvergenceState()
    converged = False
    for record in old_records:
        if record.identity_to_previous is not None:
            converged, convergence_state = check_convergence(
                record.identity_to_previous,
                threshold=config.convergence_identity,
                required_streak=config.convergence_streak,
                state=convergence_state,
            )

    return _ResumeState(
        next_iteration=last.iteration + 1,
        last_consensus_path=consensus_path,
        last_consensus_sequence=fasta_records[0].sequence,
        convergence_state=convergence_state,
        converged=converged,
        records=list(old_records),
        total_elapsed_seconds=summary["total_elapsed_seconds"],
    )


def run(
    config: Config,
    out_dir: Path,
    *,
    on_iteration: Callable[[IterationRecord], None] | None = None,
) -> RunResult:
    """Run the iterate-until-convergence loop (or resume one, see module docstring).

    If given, `on_iteration` is called once per iteration actually executed
    in this call -- not for iterations loaded from a resumed prior run --
    right after that iteration's IterationRecord is finalized. Intended for
    progress reporting (e.g. the CLI's `--progress`); errors it raises
    propagate and abort the run.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    resume = _load_resume_state(out_dir, config)
    if resume is not None and resume.converged:
        write_report(out_dir)
        return RunResult(
            output_dir=out_dir,
            iterations=tuple(resume.records),
            converged=True,
            total_elapsed_seconds=resume.total_elapsed_seconds,
        )

    cat_cache = ReadsCatCache(out_dir / "reads")
    state = _prepare_initial_state(config, out_dir)

    if resume is not None:
        start_iteration = resume.next_iteration
        current_reference_path = resume.last_consensus_path
        current_bam_path: Path | None = None
        previous_sequence: str | None = resume.last_consensus_sequence
        convergence_state = resume.convergence_state
        records: list[IterationRecord] = list(resume.records)
        prior_elapsed = resume.total_elapsed_seconds
    else:
        start_iteration = 0
        current_reference_path = state.reference_path
        current_bam_path = state.starting_bam
        previous_sequence = None
        convergence_state = ConvergenceState()
        records = []
        prior_elapsed = 0.0

    converged = False
    run_start = time.monotonic()

    for iteration in range(start_iteration, config.max_iterations + 1):
        iter_dir = out_dir / f"iter_{iteration:03d}"
        t0 = time.monotonic()
        needs_mapping = state.iteration_0_needs_mapping or iteration >= 1

        if needs_mapping:
            base_values: dict[str, object] = {
                "reference": str(current_reference_path),
                "threads": config.threads,
                **state.reads_values,
            }
            current_bam_path = _run_mapping(config, iter_dir, base_values, cat_cache, log=True)

        assert current_bam_path is not None
        reads_mapped = count_mapped_reads(current_bam_path)

        consensus_values: dict[str, object] = {
            "bam": str(current_bam_path),
            "consensus_prefix": str(iter_dir / "consensus"),
            "threads": config.threads,
        }
        if current_reference_path is not None:
            consensus_values["reference"] = str(current_reference_path)

        consensus_result = run_consensus(
            config.consensus,
            consensus_values,
            log_dir=iter_dir / "logs",
            cat_resolver=cat_cache.resolve,
        )
        consensus_copy = iter_dir / "consensus.fasta"
        write_fasta(consensus_copy, consensus_result.record_id, consensus_result.sequence)

        identity: float | None = None
        if previous_sequence is not None:
            identity = sequence_identity(previous_sequence, consensus_result.sequence).identity

        elapsed = time.monotonic() - t0
        record = IterationRecord(
            iteration=iteration,
            reads_mapped=reads_mapped,
            consensus_length=consensus_result.length,
            identity_to_previous=identity,
            elapsed_seconds=elapsed,
        )
        records.append(record)
        _write_iteration_stats(iter_dir, consensus_result, reads_mapped, identity, elapsed)
        if on_iteration is not None:
            on_iteration(record)

        if identity is not None:
            converged, convergence_state = check_convergence(
                identity,
                threshold=config.convergence_identity,
                required_streak=config.convergence_streak,
                state=convergence_state,
            )

        previous_sequence = consensus_result.sequence
        current_reference_path = consensus_copy

        if converged:
            break

    total_elapsed = prior_elapsed + (time.monotonic() - run_start)
    _write_metrics_tsv(out_dir, records)
    _write_summary_json(out_dir, records, converged, total_elapsed)
    write_report(out_dir)
    return RunResult(
        output_dir=out_dir,
        iterations=tuple(records),
        converged=converged,
        total_elapsed_seconds=total_elapsed,
    )


def _preview_iteration(
    config: Config,
    iter_dir: Path,
    *,
    needs_mapping: bool,
    reference_path: Path | None,
    bam_path: Path | None,
    reads_values: dict[str, object],
    cat_cache: ReadsCatCache,
) -> list[str]:
    """Render (but do not run) one iteration's commands."""
    lines: list[str] = []

    if needs_mapping:
        base_values: dict[str, object] = {
            "reference": str(reference_path),
            "threads": config.threads,
            **reads_values,
        }
        steps, bam_path = _plan_mapping(config, iter_dir, base_values, cat_cache)
        for step in steps:
            lines.append(step.index_rendered.display)
            lines.append(step.map_rendered.display)

    assert bam_path is not None
    consensus_values: dict[str, object] = {
        "bam": str(bam_path),
        "consensus_prefix": str(iter_dir / "consensus"),
        "threads": config.threads,
    }
    if reference_path is not None:
        consensus_values["reference"] = str(reference_path)
    for step in config.consensus.steps:
        rendered = render_command(step, consensus_values, cat_resolver=cat_cache.resolve)
        lines.append(rendered.display)

    return lines


def preview(config: Config, out_dir: Path) -> DryRunPreview:
    """Render (but do not run) the commands for iterations 0 and 1, for `--dry-run`.

    Both always run -- config validation requires max_iterations >= 1 -- so
    both can be shown in full, even though iteration 1's reference
    (iteration 0's consensus.fasta) doesn't exist yet: its *path* is fixed
    by the config, only its *content* depends on actually running iteration
    0. Whether the loop continues past iteration 1 depends on convergence,
    which can only be determined by actually running the pipeline.
    """
    cat_cache = ReadsCatCache(out_dir / "reads")
    state = _prepare_initial_state(config, out_dir)

    iter0_dir = out_dir / "iter_000"
    lines0 = _preview_iteration(
        config,
        iter0_dir,
        needs_mapping=state.iteration_0_needs_mapping,
        reference_path=state.reference_path,
        bam_path=state.starting_bam,
        reads_values=state.reads_values,
        cat_cache=cat_cache,
    )

    iter1_dir = out_dir / "iter_001"
    lines1 = _preview_iteration(
        config,
        iter1_dir,
        needs_mapping=True,
        reference_path=iter0_dir / "consensus.fasta",
        bam_path=None,
        reads_values=state.reads_values,
        cat_cache=cat_cache,
    )

    if state.iteration_0_needs_mapping:
        iteration_0_summary = "maps against the reference you gave"
    else:
        iteration_0_summary = "calls a consensus directly from the input BAM, with no mapping step"
    note = (
        f"Iterations 0 and 1 both always run. Iteration 0 {iteration_0_summary}; "
        "iteration 1 maps against iteration 0's consensus for the first time. "
        "Whether the loop continues past iteration 1 depends on convergence, "
        "which can only be determined by actually running the pipeline -- "
        "iteration 2 onward would repeat iteration 1's shape against each new "
        "consensus."
    )
    return DryRunPreview(lines=tuple(lines0 + lines1), note=note)
