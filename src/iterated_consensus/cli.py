"""Thin CLI: parses arguments and calls into the library. See `runner.py`."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Annotated

import typer

from .config import InputOverrides, apply_input_overrides, load_config, parse_file_list
from .errors import IteratedConsensusError
from .presets import PresetError, get_preset_text, list_presets
from .report import format_ambiguous_count, format_elapsed, format_identity
from .runner import IterationRecord, preview, run

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=(
        "Iteratively call a consensus sequence from a BAM or FASTQ file: "
        "map, call a consensus, remap against it, and repeat until convergence."
    ),
)


_PROGRESS_COLUMNS = (
    ("iter", 4, ">"),
    ("reads_mapped", 13, ">"),
    ("consensus_length", 18, ">"),
    ("ambiguous", 9, ">"),
    ("identity_to_previous", 22, ">"),
    ("elapsed", 8, ">"),
    ("consensus_md5", 32, ">"),
)


def _print_progress_header() -> None:
    typer.echo("  ".join(f"{name:{align}{width}}" for name, width, align in _PROGRESS_COLUMNS))


def _print_iteration_progress(record: IterationRecord) -> None:
    values = (
        str(record.iteration),
        str(record.reads_mapped),
        str(record.consensus_length),
        format_ambiguous_count(record.ambiguous_count),
        format_identity(record.identity_to_previous),
        format_elapsed(record.elapsed_seconds),
        record.consensus_md5 or "",
    )
    typer.echo(
        "  ".join(
            f"{value:{align}{width}}"
            for value, (_, width, align) in zip(values, _PROGRESS_COLUMNS, strict=True)
        )
    )


@app.command("run")
def run_cmd(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", exists=True, readable=True, help="Pipeline TOML config.")
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output-dir",
            "-o",
            help=(
                "Output directory (created if missing). Falls back to [run].output_dir "
                "in the config if not given here."
            ),
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run", "-n", help="Print the commands iteration 0 would run, without running them."
        ),
    ] = False,
    progress: Annotated[
        bool,
        typer.Option("--progress", help="Print a one-line summary after each iteration completes."),
    ] = False,
    traceback: Annotated[
        bool,
        typer.Option(
            "--traceback",
            help="Let errors crash with a full Python traceback, instead of a short message.",
        ),
    ] = False,
    bam: Annotated[
        Path | None, typer.Option(help="Starting BAM file (overrides [input].bam in the config).")
    ] = None,
    bam_reads: Annotated[
        str | None,
        typer.Option(help="Reads to extract from --bam for remapping: ref, ref+unal, or all."),
    ] = None,
    reads_1: Annotated[
        str | None, typer.Option("-1", help="Comma-separated mate-1 FASTQ file(s).")
    ] = None,
    reads_2: Annotated[
        str | None, typer.Option("-2", help="Comma-separated mate-2 FASTQ file(s).")
    ] = None,
    reads_single: Annotated[
        str | None, typer.Option("-U", help="Comma-separated unpaired FASTQ file(s).")
    ] = None,
    reference_fasta: Annotated[
        Path | None, typer.Option(help="A local FASTA file containing the starting reference.")
    ] = None,
    reference_id: Annotated[
        str | None,
        typer.Option(
            help=(
                "The name of a sequence: a record ID within --reference-fasta (needed if it "
                "has multiple sequences), a contig name within --bam (needed if it has more "
                "than one), and/or an NCBI accession -- fetched automatically if no "
                "--reference-fasta is given and this looks like one."
            )
        ),
    ] = None,
) -> None:
    """Run the iterative consensus loop."""
    try:
        config = load_config(config_path)
        overrides = InputOverrides(
            reads_1=parse_file_list(reads_1) if reads_1 is not None else None,
            reads_2=parse_file_list(reads_2) if reads_2 is not None else None,
            reads_single=parse_file_list(reads_single) if reads_single is not None else None,
            reference_fasta=reference_fasta,
            reference_id=reference_id,
            bam=bam,
            bam_reads=bam_reads,  # type: ignore[arg-type]
        )
        config = replace(config, input=apply_input_overrides(config.input, overrides))
    except IteratedConsensusError as exc:
        if traceback:
            raise
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(1) from None

    out_dir = output_dir if output_dir is not None else config.output_dir
    if out_dir is None:
        typer.echo(
            "error: no output directory given: pass --output-dir or set [run].output_dir "
            "in the config",
            err=True,
        )
        raise typer.Exit(1)

    try:
        if dry_run:
            result = preview(config, out_dir)
            for line in result.lines:
                typer.echo(f"$ {line}")
            typer.echo(f"\n# {result.note}")
        else:
            if progress:
                if (out_dir / "summary.json").exists():
                    typer.echo(f"Resuming from existing output in {out_dir}")
                _print_progress_header()
            on_iteration = _print_iteration_progress if progress else None
            outcome = run(config, out_dir, on_iteration=on_iteration)
            if outcome.cycle is not None:
                status = f"Cycle detected (period {outcome.cycle.period})"
            elif outcome.converged:
                status = "Converged"
            else:
                status = "Stopped (max_iterations reached)"
            typer.echo(
                f"{status} after {len(outcome.iterations)} iteration(s) "
                f"in {format_elapsed(outcome.total_elapsed_seconds)}"
            )
            if outcome.cycle is not None:
                c = outcome.cycle
                typer.echo(
                    f"Iteration {c.repeat_iteration}'s consensus matches iteration "
                    f"{c.first_iteration}'s (MD5 {c.consensus_md5}) -- the run is oscillating "
                    f"rather than converging. Using iteration {c.first_iteration} as the final consensus."
                )
            typer.echo(f"Output written to {outcome.output_dir}")
            typer.echo(f"Report: {outcome.output_dir / 'index.html'}")
            if outcome.final_consensus_path is not None:
                typer.echo(f"Final consensus: {outcome.final_consensus_path}")
            if outcome.final_reference_fasta_path is not None:
                typer.echo(f"Final reference: {outcome.final_reference_fasta_path}")
            if outcome.final_reference_bam_path is not None:
                typer.echo(f"Final reference BAM: {outcome.final_reference_bam_path}")
    except IteratedConsensusError as exc:
        if traceback:
            raise
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(1) from None


@app.command("config-template")
def config_template_cmd(
    name: Annotated[
        str | None, typer.Argument(help="Preset name to print. Omit to list available presets.")
    ] = None,
) -> None:
    """Print a bundled example pipeline config, as a starting point for your own."""
    if name is None:
        for preset_name, description in list_presets():
            typer.echo(f"{preset_name}\t{description}")
        return
    try:
        # get_preset_text() already ends with its own trailing newline;
        # echo's default nl=True would add a second one, printing a stray
        # blank line at the end.
        typer.echo(get_preset_text(name), nl=False)
    except PresetError as exc:
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(1) from None


def main() -> None:
    app()


if __name__ == "__main__":
    main()
