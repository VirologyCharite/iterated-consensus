"""TOML configuration schema for mapper/consensus pipelines and run inputs."""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from .errors import IteratedConsensusError

CommandStep = str | list[str]
BamReadsMode = Literal["ref", "ref+unal", "all"]
_VALID_BAM_READS_MODES = ("ref", "ref+unal", "all")

# [run] keys with dedicated Config fields -- anything else under [run] becomes
# a custom {name} placeholder usable in mapper/consensus commands (see
# Config.extra_vars).
_KNOWN_RUN_KEYS = frozenset(
    {
        "threads", "threads_reserve", "max_iterations", "convergence_identity", "convergence_streak",
        "output_dir",
    }
)

# Placeholder names iterated-consensus fills in itself, once per iteration --
# a custom [run] variable can't reuse one of these (see Config.validate). This
# is _KNOWN_RUN_KEYS (all of which double as placeholders -- see runner.py)
# plus the placeholders that aren't [run] keys at all.
_RESERVED_PLACEHOLDER_NAMES = _KNOWN_RUN_KEYS | frozenset(
    {
        "reference", "index_prefix", "bam", "consensus_prefix", "reads_1", "reads_2", "reads_single",
        "consensus_fasta", "consensus_id",
    }
)


def available_cpu_count() -> int:
    """CPUs actually usable by this process, not just installed.

    Prefers the scheduler affinity mask where the platform exposes one
    (Linux -- also respects cgroup/container CPU limits there), since that's
    a tighter bound than the physical core count when running inside a
    container or under `taskset`. Falls back to the total installed count
    elsewhere (e.g. macOS, which has no `sched_getaffinity`).
    """
    try:
        return len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except AttributeError:
        return os.cpu_count() or 1


class ConfigError(IteratedConsensusError, ValueError):
    """Raised for malformed or inconsistent configuration."""


@dataclass(frozen=True)
class Mapper:
    name: str
    index_cmd: CommandStep
    map_cmd: CommandStep
    tool_versions: dict[str, CommandStep] = field(default_factory=dict)
    """From this mapper's [tool-versions] sub-table, if given: tool name ->
    a command whose stdout is that tool's version string. Run once per
    iteration, immediately before index_cmd/map_cmd -- see runner.py."""


@dataclass(frozen=True)
class ConsensusSpec:
    steps: tuple[CommandStep, ...]
    output: str
    tool_versions: dict[str, CommandStep] = field(default_factory=dict)
    """From [consensus]'s [tool-versions] sub-table, if given -- see
    Mapper.tool_versions above; run immediately before the consensus steps."""


@dataclass(frozen=True)
class InputSpec:
    reads_1: tuple[Path, ...] = ()
    reads_2: tuple[Path, ...] = ()
    reads_single: tuple[Path, ...] = ()
    reference_fasta: Path | None = None
    reference_id: str | None = None
    bam: Path | None = None
    bam_reads: BamReadsMode = "ref"

    def validate(self) -> None:
        have_fastq = bool(self.reads_1 or self.reads_2 or self.reads_single)
        have_bam = self.bam is not None
        if have_fastq and have_bam:
            raise ConfigError("specify FASTQ input (reads_1/reads_2/reads_single) or bam, not both")
        if not have_fastq and not have_bam:
            raise ConfigError("no input given: specify FASTQ input or bam")

        if have_fastq:
            if len(self.reads_1) != len(self.reads_2):
                raise ConfigError(
                    f"reads_1 has {len(self.reads_1)} file(s) but reads_2 has "
                    f"{len(self.reads_2)}; they must pair up 1:1"
                )
            if not self.reads_1 and not self.reads_single:
                raise ConfigError("FASTQ input needs reads_1/reads_2 pairs and/or reads_single reads")
            if self.reference_fasta is None and self.reference_id is None:
                raise ConfigError(
                    "FASTQ input needs a starting reference (reference_fasta and/or reference_id)"
                )

        # For have_bam, reference_fasta/reference_id are both optional: if
        # neither resolves to a usable reference, iteration 0 just has no
        # {reference} available (see runner.py).

        if self.bam_reads not in _VALID_BAM_READS_MODES:
            raise ConfigError(
                f"invalid bam_reads {self.bam_reads!r}; must be one of {_VALID_BAM_READS_MODES}"
            )


@dataclass(frozen=True)
class OutputSpec:
    consensus_fasta: str | None = None
    """A template string, not a plain path -- rendered against the usual
    [run] placeholders (plus {consensus_id}) before use, so it can reference
    {output_dir} (or anything else) to control where it lands. A relative
    result is *not* auto-placed under the output directory: it's just a
    normal relative path, resolved against the current working directory
    like any other file argument. Use "{output_dir}/..." explicitly to put
    it under --output-dir."""
    consensus_id: str | None = None
    commands: tuple[CommandStep, ...] = ()
    """Run once, after consensus_fasta is written, in order. {consensus_fasta}
    and {consensus_id} are available as placeholders (alongside every other
    placeholder [run] commands can use), resolved to the final values --
    the actual id used, even if consensus_id itself was left unset and the
    consensus tool's own id was kept."""

    def validate(self) -> None:
        if self.consensus_id is not None and self.consensus_fasta is None:
            raise ConfigError("[output] consensus_id needs consensus_fasta to also be given")
        if self.commands and self.consensus_fasta is None:
            raise ConfigError("[output] commands needs consensus_fasta to also be given")


@dataclass(frozen=True)
class InputOverrides:
    """CLI-supplied input fields; each `None` means "leave the config value as-is"."""

    reads_1: tuple[Path, ...] | None = None
    reads_2: tuple[Path, ...] | None = None
    reads_single: tuple[Path, ...] | None = None
    reference_fasta: Path | None = None
    reference_id: str | None = None
    bam: Path | None = None
    bam_reads: BamReadsMode | None = None


def parse_file_list(value: str) -> tuple[Path, ...]:
    """Split a bowtie2-style comma-separated file list into paths."""
    return tuple(Path(v) for v in value.split(",") if v)


def apply_input_overrides(base: InputSpec | None, overrides: InputOverrides) -> InputSpec:
    """Layer CLI-supplied input fields over a config's [input] section (or a blank one)."""
    base = base if base is not None else InputSpec()
    merged = InputSpec(
        reads_1=overrides.reads_1 if overrides.reads_1 is not None else base.reads_1,
        reads_2=overrides.reads_2 if overrides.reads_2 is not None else base.reads_2,
        reads_single=overrides.reads_single if overrides.reads_single is not None else base.reads_single,
        reference_fasta=(
            overrides.reference_fasta if overrides.reference_fasta is not None else base.reference_fasta
        ),
        reference_id=overrides.reference_id if overrides.reference_id is not None else base.reference_id,
        bam=overrides.bam if overrides.bam is not None else base.bam,
        bam_reads=overrides.bam_reads if overrides.bam_reads is not None else base.bam_reads,
    )
    merged.validate()
    return merged


@dataclass(frozen=True)
class Config:
    mappers: tuple[Mapper, ...]
    consensus: ConsensusSpec
    input: InputSpec | None = None
    output: OutputSpec | None = None
    output_dir: Path | None = None
    """Fallback for --output-dir: used only when the CLI flag isn't given.
    Also available as the {output_dir} placeholder in every command, set to
    whichever of the two actually applies."""
    threads: int = 1
    threads_reserve: int = 0
    """Only meaningful (and only settable) alongside threads = "auto" -- see
    _resolve_threads. Kept here (rather than discarded once resolved) purely
    so it's available as a {threads_reserve} placeholder like the rest of
    [run], e.g. for a logging step."""
    extra_vars: dict[str, str | int | float | bool] = field(default_factory=dict)
    """Custom [run] variables (anything under [run] besides the dedicated
    fields on this class), usable as {name} placeholders in every
    mapper/consensus command -- same as {threads}, {max_iterations}, etc."""
    max_iterations: int = 20
    convergence_identity: float = 100.0
    convergence_streak: int = 1

    def validate(self) -> None:
        if not self.mappers:
            raise ConfigError("at least one [[mapper]] must be given")
        names = [m.name for m in self.mappers]
        if len(names) != len(set(names)):
            raise ConfigError(f"mapper names must be unique, got {names}")
        if not self.consensus.steps:
            raise ConfigError("[consensus] must have at least one step")
        if not self.consensus.output:
            raise ConfigError("[consensus] must set 'output'")
        if self.threads < 1:
            raise ConfigError(f"threads must be >= 1, got {self.threads}")
        reserved_used = sorted(set(self.extra_vars) & _RESERVED_PLACEHOLDER_NAMES)
        if reserved_used:
            raise ConfigError(
                f"[run] variable(s) {reserved_used} collide with placeholder(s) "
                "iterated-consensus sets automatically each iteration; choose different name(s)"
            )
        if self.max_iterations < 1:
            raise ConfigError(f"max_iterations must be >= 1, got {self.max_iterations}")
        if not (0 < self.convergence_identity <= 100):
            raise ConfigError(
                f"convergence_identity must be in (0, 100], got {self.convergence_identity}"
            )
        if self.convergence_streak < 1:
            raise ConfigError(f"convergence_streak must be >= 1, got {self.convergence_streak}")
        if self.input is not None:
            self.input.validate()
        if self.output is not None:
            self.output.validate()


def _command_step(value: object, where: str) -> CommandStep:
    if isinstance(value, str):
        return value
    if isinstance(value, list) and all(isinstance(v, str) for v in value):
        return list(value)
    raise ConfigError(f"{where} must be a string or a list of strings, got {value!r}")


def _parse_tool_versions(raw: dict, where: str) -> dict[str, CommandStep]:
    tv_raw = raw.get("tool-versions", {})
    if not isinstance(tv_raw, dict):
        raise ConfigError(f"{where} [tool-versions] must be a table")
    return {
        name: _command_step(cmd, f"{where} [tool-versions] '{name}'") for name, cmd in tv_raw.items()
    }


def _paths(value: object, where: str) -> tuple[Path, ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise ConfigError(f"{where} must be a list of strings, got {value!r}")
    return tuple(Path(v) for v in value)


def _parse_mapper(raw: dict, index: int) -> Mapper:
    try:
        name = raw["name"]
    except KeyError as exc:
        raise ConfigError(f"[[mapper]] #{index} is missing 'name'") from exc
    if not isinstance(name, str):
        raise ConfigError(f"[[mapper]] #{index} 'name' must be a string")
    try:
        index_cmd = raw["index_cmd"]
        map_cmd = raw["map_cmd"]
    except KeyError as exc:
        raise ConfigError(f"mapper '{name}' is missing {exc.args[0]!r}") from exc
    return Mapper(
        name=name,
        index_cmd=_command_step(index_cmd, f"mapper '{name}' index_cmd"),
        map_cmd=_command_step(map_cmd, f"mapper '{name}' map_cmd"),
        tool_versions=_parse_tool_versions(raw, f"mapper '{name}'"),
    )


def _parse_consensus(raw: dict) -> ConsensusSpec:
    steps_raw = raw.get("steps", [])
    if not isinstance(steps_raw, list):
        raise ConfigError("[consensus] steps must be a list")
    steps = tuple(_command_step(s, f"[consensus] steps[{i}]") for i, s in enumerate(steps_raw))
    output = raw.get("output", "")
    if not isinstance(output, str):
        raise ConfigError("[consensus] output must be a string")
    return ConsensusSpec(steps=steps, output=output, tool_versions=_parse_tool_versions(raw, "[consensus]"))


def _parse_input(raw: dict) -> InputSpec:
    reference_fasta = raw.get("reference_fasta")
    bam = raw.get("bam")
    return InputSpec(
        reads_1=_paths(raw.get("reads_1"), "[input] reads_1"),
        reads_2=_paths(raw.get("reads_2"), "[input] reads_2"),
        reads_single=_paths(raw.get("reads_single"), "[input] reads_single"),
        reference_fasta=Path(reference_fasta) if reference_fasta is not None else None,
        reference_id=raw.get("reference_id"),
        bam=Path(bam) if bam is not None else None,
        bam_reads=raw.get("bam_reads", "ref"),
    )


def _parse_output(raw: dict) -> OutputSpec:
    consensus_fasta = raw.get("consensus_fasta")
    if consensus_fasta is not None and not isinstance(consensus_fasta, str):
        raise ConfigError(f"[output] consensus_fasta must be a string, got {consensus_fasta!r}")
    commands_raw = raw.get("commands", [])
    if not isinstance(commands_raw, list):
        raise ConfigError("[output] commands must be a list")
    commands = tuple(_command_step(c, f"[output] commands[{i}]") for i, c in enumerate(commands_raw))
    return OutputSpec(
        consensus_fasta=consensus_fasta,
        consensus_id=raw.get("consensus_id"),
        commands=commands,
    )


def _resolve_threads(run_raw: dict) -> int:
    value = run_raw.get("threads", 1)
    if value == "auto":
        reserve = run_raw.get("threads_reserve", 0)
        if not isinstance(reserve, int) or isinstance(reserve, bool) or reserve < 0:
            raise ConfigError(
                f"[run] threads_reserve must be a non-negative integer, got {reserve!r}"
            )
        # Leave at least 1 thread even if threads_reserve meets or exceeds
        # the detected count -- a run that can't use any threads can't run.
        return max(1, available_cpu_count() - reserve)
    if "threads_reserve" in run_raw:
        raise ConfigError('[run] threads_reserve only applies when threads = "auto"')
    if not isinstance(value, int) or isinstance(value, bool):
        raise ConfigError(f'[run] threads must be an integer or "auto", got {value!r}')
    return value


def _parse_output_dir(run_raw: dict) -> Path | None:
    value = run_raw.get("output_dir")
    if value is None:
        return None
    if not isinstance(value, str):
        raise ConfigError(f"[run] output_dir must be a string, got {value!r}")
    return Path(value)


def _parse_run_extras(run_raw: dict) -> dict[str, str | int | float | bool]:
    extras: dict[str, str | int | float | bool] = {}
    for key, value in run_raw.items():
        if key in _KNOWN_RUN_KEYS:
            continue
        if not isinstance(value, str | int | float | bool):
            raise ConfigError(
                f"[run] '{key}' must be a string, number, or boolean to be used as a "
                f"command placeholder, got {value!r}"
            )
        extras[key] = value
    return extras


def parse_config(text: str) -> Config:
    """Parse and validate config TOML source (does not touch the filesystem)."""
    try:
        raw = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise ConfigError(f"invalid TOML: {exc}") from exc

    mappers_raw = raw.get("mapper", [])
    if not isinstance(mappers_raw, list):
        raise ConfigError("[[mapper]] must be a list of tables")
    mappers = tuple(_parse_mapper(m, i) for i, m in enumerate(mappers_raw))

    consensus_raw = raw.get("consensus")
    if consensus_raw is None:
        raise ConfigError("missing [consensus] section")
    consensus = _parse_consensus(consensus_raw)

    input_raw = raw.get("input")
    input_spec = _parse_input(input_raw) if input_raw is not None else None

    output_raw = raw.get("output")
    output_spec = _parse_output(output_raw) if output_raw is not None else None

    run_raw = raw.get("run", {})
    if not isinstance(run_raw, dict):
        raise ConfigError("[run] must be a table")
    config = Config(
        mappers=mappers,
        consensus=consensus,
        input=input_spec,
        output=output_spec,
        output_dir=_parse_output_dir(run_raw),
        threads=_resolve_threads(run_raw),
        threads_reserve=run_raw.get("threads_reserve", 0),  # validated inside _resolve_threads above
        extra_vars=_parse_run_extras(run_raw),
        max_iterations=run_raw.get("max_iterations", 20),
        convergence_identity=run_raw.get("convergence_identity", 100.0),
        convergence_streak=run_raw.get("convergence_streak", 1),
    )
    config.validate()
    return config


def load_config(path: Path) -> Config:
    return parse_config(path.read_text())
