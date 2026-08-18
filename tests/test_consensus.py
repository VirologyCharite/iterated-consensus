from pathlib import Path

import pytest

from iterated_consensus.config import ConsensusSpec
from iterated_consensus.consensus import ConsensusError, run_consensus


def test_run_consensus_single_step(tmp_path: Path) -> None:
    output = tmp_path / "consensus.fa"
    spec = ConsensusSpec(
        steps=('printf ">seq1\\nACGTACGT\\n" > {output}',),
        output="{output}",
    )
    result = run_consensus(spec, {"output": str(output)})
    assert result.record_id == "seq1"
    assert result.sequence == "ACGTACGT"
    assert result.length == 8
    assert result.composition == {"A": 2, "C": 2, "G": 2, "T": 2}
    assert result.fasta_path == output


def test_run_consensus_multiple_steps_run_in_order(tmp_path: Path) -> None:
    intermediate = tmp_path / "pileup.txt"
    output = tmp_path / "out.fa"
    spec = ConsensusSpec(
        steps=(
            f'printf "pileup-data" > {intermediate}',
            'printf ">c1\\nGGGG\\n" > {output}',
        ),
        output="{output}",
    )
    result = run_consensus(spec, {"output": str(output)})
    assert intermediate.read_text() == "pileup-data"
    assert result.sequence == "GGGG"


def test_run_consensus_missing_output_raises(tmp_path: Path) -> None:
    spec = ConsensusSpec(steps=("true",), output=str(tmp_path / "nope.fa"))
    with pytest.raises(ConsensusError, match="was not created"):
        run_consensus(spec, {})


def test_run_consensus_failing_step_raises(tmp_path: Path) -> None:
    spec = ConsensusSpec(steps=("exit 1",), output=str(tmp_path / "out.fa"))
    with pytest.raises(ConsensusError, match="consensus step 0 failed"):
        run_consensus(spec, {})


def test_run_consensus_multi_record_output_raises(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    spec = ConsensusSpec(
        steps=('printf ">a\\nAAAA\\n>b\\nCCCC\\n" > {output}',),
        output="{output}",
    )
    with pytest.raises(ConsensusError, match="expected exactly one sequence"):
        run_consensus(spec, {"output": str(output)})


def test_run_consensus_writes_logs(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    log_dir = tmp_path / "logs"
    spec = ConsensusSpec(
        steps=('printf ">a\\nAAAA\\n" > {output}',),
        output="{output}",
    )
    run_consensus(spec, {"output": str(output)}, log_dir=log_dir)
    assert (log_dir / "consensus_step_00.log").exists()


def test_run_consensus_no_log_dir_leaves_commands_empty(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    spec = ConsensusSpec(steps=('printf ">a\\nAAAA\\n" > {output}',), output="{output}")
    result = run_consensus(spec, {"output": str(output)})
    assert result.commands == ()


def test_run_consensus_with_log_dir_records_commands(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    log_dir = tmp_path / "logs"
    spec = ConsensusSpec(
        steps=('printf ">a\\nAAAA\\n" > {output}', 'echo "second step"'),
        output="{output}",
    )
    result = run_consensus(spec, {"output": str(output)}, log_dir=log_dir)
    assert [c.name for c in result.commands] == ["consensus_step_00", "consensus_step_01"]
    assert "second step" in result.commands[1].log


def test_run_consensus_runs_tool_versions_and_records_them(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    spec = ConsensusSpec(
        steps=('printf ">a\\nAAAA\\n" > {output}',),
        output="{output}",
        tool_versions={"ivar": "echo 'ivar 1.4.2'"},
    )
    result = run_consensus(spec, {"output": str(output)})
    assert result.tool_versions == {"ivar": "ivar 1.4.2"}


def test_run_consensus_tool_versions_failure_raises_consensus_error(tmp_path: Path) -> None:
    output = tmp_path / "out.fa"
    spec = ConsensusSpec(
        steps=('printf ">a\\nAAAA\\n" > {output}',),
        output="{output}",
        tool_versions={"ivar": "exit 1"},
    )
    with pytest.raises(ConsensusError, match="tool-versions 'ivar' failed"):
        run_consensus(spec, {"output": str(output)})
