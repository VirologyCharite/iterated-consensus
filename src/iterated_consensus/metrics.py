"""Per-iteration metrics: sequence identity, convergence, base composition."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass

import edlib

_CIGAR_OP_RE = re.compile(r"(\d+)([=XIDM])")


@dataclass(frozen=True)
class IdentityResult:
    identity: float  # percent, 0-100
    alignment_length: int
    edit_distance: int


def sequence_identity(a: str, b: str) -> IdentityResult:
    """Global-alignment percent identity between two sequences, via edlib."""
    if not a or not b:
        raise ValueError("cannot compute identity of an empty sequence")
    result = edlib.align(a, b, mode="NW", task="path")
    edit_distance: int = result["editDistance"]
    alignment_length = sum(int(n) for n, _op in _CIGAR_OP_RE.findall(result["cigar"]))
    identity = 100.0 * (alignment_length - edit_distance) / alignment_length
    return IdentityResult(
        identity=identity,
        alignment_length=alignment_length,
        edit_distance=edit_distance,
    )


@dataclass(frozen=True)
class ConvergenceState:
    streak: int = 0


def check_convergence(
    identity: float,
    *,
    threshold: float,
    required_streak: int,
    state: ConvergenceState,
) -> tuple[bool, ConvergenceState]:
    """Update a convergence streak with the latest identity value.

    Returns (converged, new_state). `converged` is True once `identity` has
    been >= `threshold` for `required_streak` consecutive calls.
    """
    streak = state.streak + 1 if identity >= threshold else 0
    new_state = ConvergenceState(streak=streak)
    return streak >= required_streak, new_state


def base_composition(seq: str) -> dict[str, int]:
    """Count of each character (base/ambiguity code/gap) in a sequence."""
    return dict(Counter(seq.upper()))
