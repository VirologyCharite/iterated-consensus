"""The `{cat:...}` materialization cache: concatenate a ReadsList on first use."""

from __future__ import annotations

import shutil
from pathlib import Path

from .templating import ReadsList


class ReadsCatCache:
    """Lazily concatenates a ReadsList's files into one, then reuses the result.

    Used as the `cat_resolver` passed to `templating.render`/`commands.render_command`
    so that `{cat:name}` only pays the concatenation I/O cost for mappers that
    actually need a single-file input, and only once per run (not once per
    iteration -- the original reads never change between iterations).
    """

    def __init__(self, work_dir: Path, *, dry_run: bool = False):
        self._work_dir = work_dir
        self._cache: dict[str, Path] = {}
        self._dry_run = dry_run

    def resolve(self, reads_list: ReadsList) -> Path:
        cached = self._cache.get(reads_list.name)
        if cached is not None:
            return cached

        if len(reads_list.paths) == 1:
            target = reads_list.paths[0]
        else:
            target = self._concatenate(reads_list)

        self._cache[reads_list.name] = target
        return target

    def _concatenate(self, reads_list: ReadsList) -> Path:
        is_gz = {p.suffix == ".gz" for p in reads_list.paths}
        if len(is_gz) != 1:
            raise ValueError(
                f"reads list '{reads_list.name}' mixes gzipped and non-gzipped files; "
                "cannot concatenate"
            )
        suffix = ".fastq.gz" if is_gz.pop() else ".fastq"
        target = self._work_dir / f"{reads_list.name}.cat{suffix}"
        if self._dry_run:
            # Report the path this would end up at, without actually
            # concatenating anything.
            return target
        if not target.exists():
            self._work_dir.mkdir(parents=True, exist_ok=True)
            with target.open("wb") as out:
                for path in reads_list.paths:
                    with path.open("rb") as src:
                        shutil.copyfileobj(src, out)
        return target
