import pytest

from iterated_consensus.metrics import (
    ConvergenceState,
    base_composition,
    check_convergence,
    sequence_identity,
)


def test_identical_sequences_are_100_percent() -> None:
    result = sequence_identity("ACGTACGT", "ACGTACGT")
    assert result.identity == 100.0
    assert result.edit_distance == 0
    assert result.alignment_length == 8


def test_single_substitution() -> None:
    result = sequence_identity("ACGTACGT", "ACGTACGA")
    assert result.edit_distance == 1
    assert result.identity == pytest.approx(87.5)


def test_insertion_changes_alignment_length() -> None:
    result = sequence_identity("ACGT", "ACGGT")
    assert result.edit_distance == 1
    assert result.alignment_length == 5
    assert result.identity == pytest.approx(80.0)


def test_completely_different_sequences_low_identity() -> None:
    result = sequence_identity("AAAA", "TTTT")
    assert result.edit_distance == 4
    assert result.identity == 0.0


def test_empty_sequence_raises() -> None:
    with pytest.raises(ValueError, match="empty"):
        sequence_identity("", "ACGT")
    with pytest.raises(ValueError, match="empty"):
        sequence_identity("ACGT", "")


def test_convergence_streak_builds_and_resets() -> None:
    state = ConvergenceState()
    converged, state = check_convergence(100.0, threshold=100.0, required_streak=2, state=state)
    assert not converged
    assert state.streak == 1

    converged, state = check_convergence(100.0, threshold=100.0, required_streak=2, state=state)
    assert converged
    assert state.streak == 2


def test_convergence_streak_resets_on_drop() -> None:
    state = ConvergenceState(streak=1)
    converged, state = check_convergence(90.0, threshold=100.0, required_streak=2, state=state)
    assert not converged
    assert state.streak == 0


def test_convergence_single_streak_default() -> None:
    state = ConvergenceState()
    converged, state = check_convergence(99.9, threshold=99.5, required_streak=1, state=state)
    assert converged
    assert state.streak == 1


def test_base_composition_counts_and_uppercases() -> None:
    assert base_composition("acgtACGTN") == {
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 1,
    }


def test_base_composition_empty_sequence() -> None:
    assert base_composition("") == {}
