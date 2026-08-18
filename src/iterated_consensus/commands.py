"""Rendering and running a single pipeline step (mapper index/map, consensus step)."""

from __future__ import annotations

import shlex
import subprocess
import time
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from .config import CommandStep
from .errors import IteratedConsensusError
from .templating import CatResolver, render


class CommandError(IteratedConsensusError, RuntimeError):
    """Raised when a rendered command exits non-zero."""


@dataclass(frozen=True)
class RenderedCommand:
    argv_or_shell: list[str] | str
    display: str

    @property
    def is_shell(self) -> bool:
        return isinstance(self.argv_or_shell, str)


def render_command(
    step: CommandStep,
    values: Mapping[str, object],
    *,
    cat_resolver: CatResolver | None = None,
) -> RenderedCommand:
    """Render a config CommandStep (list or shell string) against `values`."""
    if isinstance(step, list):
        rendered = [render(token, values, cat_resolver=cat_resolver) for token in step]
        return RenderedCommand(argv_or_shell=rendered, display=shlex.join(rendered))
    rendered_str = render(step, values, cat_resolver=cat_resolver)
    return RenderedCommand(argv_or_shell=rendered_str, display=rendered_str)


def run_command(
    cmd: RenderedCommand,
    *,
    log_path: Path | None = None,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess:
    """Execute a rendered command, optionally teeing stdout+stderr to a log file."""
    if log_path is None:
        result = subprocess.run(cmd.argv_or_shell, shell=cmd.is_shell, cwd=cwd, check=False)
    else:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("wb") as log_file:
            log_file.write(f"$ {cmd.display}\n".encode())
            log_file.flush()
            result = subprocess.run(
                cmd.argv_or_shell,
                shell=cmd.is_shell,
                cwd=cwd,
                check=False,
                stdout=log_file,
                stderr=subprocess.STDOUT,
            )

    if result.returncode != 0:
        where = f" (see {log_path})" if log_path is not None else ""
        raise CommandError(
            f"command failed (exit {result.returncode}): {cmd.display}{where}"
        )
    return result


@dataclass(frozen=True)
class CommandRun:
    """A record of one executed command, for surfacing in stats.json and the
    report -- unlike RenderedCommand, this is about what actually happened,
    not what was about to be run."""

    name: str
    display: str
    elapsed_seconds: float
    log: str


def run_logged_command(name: str, cmd: RenderedCommand, log_dir: Path) -> CommandRun:
    """Like run_command, but times the run and reads the log file's own text
    back in (not just its path), for embedding in stats.json / the report."""
    log_path = log_dir / f"{name}.log"
    t0 = time.monotonic()
    run_command(cmd, log_path=log_path)
    elapsed = time.monotonic() - t0
    return CommandRun(name=name, display=cmd.display, elapsed_seconds=elapsed, log=log_path.read_text())


def run_tool_version_command(step: CommandStep, values: Mapping[str, object]) -> str:
    """Run a [tool-versions] command and return its stdout, stripped -- may
    legitimately span multiple lines. Not written to a log file (the result
    goes straight into stats.json instead); only stdout is captured, so a
    tool that prints its version to stderr needs its own `2>&1` in the
    command, same as it would need `head`/`cut`/etc. to trim the output."""
    rendered = render_command(step, values)
    try:
        result = subprocess.run(
            rendered.argv_or_shell, shell=rendered.is_shell, capture_output=True, text=True, check=False
        )
    except OSError as exc:
        raise CommandError(f"could not run tool-versions command '{rendered.display}': {exc}") from exc
    if result.returncode != 0:
        detail = f"\n{result.stderr}" if result.stderr else ""
        raise CommandError(
            f"tool-versions command failed (exit {result.returncode}): {rendered.display}{detail}"
        )
    return result.stdout.strip()
