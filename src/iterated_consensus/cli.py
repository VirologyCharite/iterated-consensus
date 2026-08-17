"""Thin CLI: parses arguments and calls into the library. See `runner.py`."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Annotated

import typer

from .config import (
    ConfigError,
    InputOverrides,
    apply_input_overrides,
    load_config,
    parse_file_list,
)
from .presets import PresetError, get_preset_text, list_presets
from .report import format_elapsed, format_identity
from .runner import IterationRecord, RunnerError, preview, run

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help=(
        "Iteratively call a consensus sequence from a BAM or FASTQ file: "
        "map, call a consensus, remap against it, and repeat until convergence."
    ),
)


def _print_iteration_progress(record: IterationRecord) -> None:
    typer.echo(
        f"iter {record.iteration:>3}  "
        f"reads_mapped={record.reads_mapped:<8} "
        f"consensus_length={record.consensus_length:<8} "
        f"identity_to_previous={format_identity(record.identity_to_previous):<9} "
        f"elapsed={format_elapsed(record.elapsed_seconds)}"
    )


@app.command("run")
def run_cmd(
    config_path: Annotated[
        Path, typer.Option("--config", "-c", exists=True, readable=True, help="Pipeline TOML config.")
    ],
    out_dir: Annotated[
        Path, typer.Option("--out-dir", "-o", help="Output directory (created if missing).")
    ],
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
    bam: Annotated[
        Path | None, typer.Option(help="Starting BAM file (overrides [input].bam in the config).")
    ] = None,
    reference_name: Annotated[
        str | None,
        typer.Option(help="Reference/contig in --bam to call a consensus from (needed if it has >1)."),
    ] = None,
    bam_reads: Annotated[
        str | None,
        typer.Option(help="Reads to extract from --bam for remapping: ref, ref+unal, or all."),
    ] = None,
    mate1: Annotated[
        str | None, typer.Option("-1", help="Comma-separated mate-1 FASTQ file(s).")
    ] = None,
    mate2: Annotated[
        str | None, typer.Option("-2", help="Comma-separated mate-2 FASTQ file(s).")
    ] = None,
    unpaired: Annotated[
        str | None, typer.Option("-U", help="Comma-separated unpaired FASTQ file(s).")
    ] = None,
    reference: Annotated[
        Path | None, typer.Option(help="Starting reference FASTA (FASTQ-start only).")
    ] = None,
    reference_id: Annotated[
        str | None, typer.Option(help="Record ID within --reference, if it has multiple sequences.")
    ] = None,
    accession: Annotated[
        str | None,
        typer.Option(help="NCBI accession to fetch as the starting reference (FASTQ-start only)."),
    ] = None,
) -> None:
    """Run the iterative consensus loop."""
    try:
        config = load_config(config_path)
        overrides = InputOverrides(
            mate1=parse_file_list(mate1) if mate1 is not None else None,
            mate2=parse_file_list(mate2) if mate2 is not None else None,
            unpaired=parse_file_list(unpaired) if unpaired is not None else None,
            reference=reference,
            accession=accession,
            reference_id=reference_id,
            bam=bam,
            reference_name=reference_name,
            bam_reads=bam_reads,  # type: ignore[arg-type]
        )
        config = replace(config, input=apply_input_overrides(config.input, overrides))
    except ConfigError as exc:
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(1) from None

    try:
        if dry_run:
            result = preview(config, out_dir)
            for line in result.lines:
                typer.echo(f"$ {line}")
            typer.echo(f"\n# {result.note}")
        else:
            if progress and (out_dir / "summary.json").exists():
                typer.echo(f"Resuming from existing output in {out_dir}")
            on_iteration = _print_iteration_progress if progress else None
            outcome = run(config, out_dir, on_iteration=on_iteration)
            status = "Converged" if outcome.converged else "Stopped (max_iterations reached)"
            typer.echo(
                f"{status} after {len(outcome.iterations)} iteration(s) "
                f"in {format_elapsed(outcome.total_elapsed_seconds)}"
            )
            typer.echo(f"Output written to {outcome.output_dir}")
            typer.echo(f"Report: {outcome.output_dir / 'index.html'}")
    except (RunnerError, ConfigError) as exc:
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
        typer.echo(get_preset_text(name))
    except PresetError as exc:
        typer.echo(f"error: {exc}", err=True)
        raise typer.Exit(1) from None


def main() -> None:
    app()


if __name__ == "__main__":
    main()
