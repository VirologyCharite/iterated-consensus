"""Shared base for this project's own recognized errors.

Every error type raised by this library's own code (bad config, a failed
command, a malformed template, ...) inherits from this, so the CLI can catch
all of them uniformly with one `except IteratedConsensusError` -- see
cli.py's `--traceback` handling. Each still also keeps its original base
(ValueError/RuntimeError/...) for anyone catching that directly.
"""

from __future__ import annotations


class IteratedConsensusError(Exception):
    """Base class for all errors raised by this library's own code."""
