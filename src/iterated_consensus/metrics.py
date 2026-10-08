"""Per-iteration metrics: sequence identity, convergence, base composition."""

from __future__ import annotations

import hashlib
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


UNAMBIGUOUS_BASES = ("A", "C", "G", "T")


def ambiguous_count(composition: dict[str, int]) -> int:
    """Count of characters in a base_composition() result that aren't plain
    A/C/G/T -- IUPAC ambiguity codes, N, gaps, anything else."""
    unambiguous = sum(composition.get(b, 0) for b in UNAMBIGUOUS_BASES)
    return sum(composition.values()) - unambiguous


@dataclass(frozen=True)
class AmbiguityTransitions:
    formerly_ambiguous: int
    """Sites that were ambiguous (not plain A/C/G/T) in the previous
    sequence but are unambiguous in the current one."""
    newly_ambiguous: int
    """Sites that were unambiguous in the previous sequence but are
    ambiguous in the current one."""


def ambiguity_transitions(previous: str, current: str) -> AmbiguityTransitions:
    """Count sites whose ambiguity status flipped between two sequences.

    Aligns `previous` to `current` (global alignment via edlib, same method
    as sequence_identity) and compares ambiguity status -- plain A/C/G/T vs.
    anything else -- at each aligned (match/substitution) column. Inserted
    or deleted sites have no counterpart in the other sequence, so they
    aren't counted toward either total.
    """
    if not previous or not current:
        raise ValueError("cannot compare ambiguity transitions of an empty sequence")
    cigar = edlib.align(previous, current, mode="NW", task="path")["cigar"]
    formerly_ambiguous = newly_ambiguous = 0
    i = j = 0
    for n, op in _CIGAR_OP_RE.findall(cigar):
        n = int(n)
        if op in ("=", "X", "M"):
            for _ in range(n):
                prev_ambiguous = previous[i].upper() not in UNAMBIGUOUS_BASES
                curr_ambiguous = current[j].upper() not in UNAMBIGUOUS_BASES
                if prev_ambiguous and not curr_ambiguous:
                    formerly_ambiguous += 1
                elif curr_ambiguous and not prev_ambiguous:
                    newly_ambiguous += 1
                i += 1
                j += 1
        elif op == "D":
            # A site present in `current` with no counterpart in `previous`.
            j += n
        else:  # "I"
            # A site present in `previous` with no counterpart in `current`.
            i += n
    return AmbiguityTransitions(formerly_ambiguous=formerly_ambiguous, newly_ambiguous=newly_ambiguous)


def sequence_md5(seq: str) -> str:
    """Hex-digest MD5 of a sequence, case-normalized and with no whitespace --
    the convention used by NCBI/ENA for sequence checksums, so it's directly
    comparable to hashes computed elsewhere, not just within this tool."""
    normalized = "".join(seq.split()).upper()
    return hashlib.md5(normalized.encode("ascii"), usedforsecurity=False).hexdigest()
