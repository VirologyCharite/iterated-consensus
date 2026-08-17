"""Rendering and running a single pipeline step (mapper index/map, consensus step)."""

from __future__ import annotations

import shlex
import subprocess
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from .config import CommandStep
from .templating import CatResolver, render


class CommandError(RuntimeError):
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
