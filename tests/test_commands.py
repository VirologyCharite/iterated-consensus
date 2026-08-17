from pathlib import Path

import pytest

from iterated_consensus.commands import CommandError, render_command, run_command


def test_render_list_step() -> None:
    rendered = render_command(["bowtie2-build", "{reference}", "{index_prefix}"], {
        "reference": "ref.fa",
        "index_prefix": "idx",
    })
    assert rendered.argv_or_shell == ["bowtie2-build", "ref.fa", "idx"]
    assert not rendered.is_shell


def test_render_shell_step() -> None:
    rendered = render_command(
        "bowtie2 -x {index_prefix} -p {threads} | samtools sort -o {bam}",
        {"index_prefix": "idx", "threads": 4, "bam": "out.bam"},
    )
    assert rendered.argv_or_shell == "bowtie2 -x idx -p 4 | samtools sort -o out.bam"
    assert rendered.is_shell


def test_run_list_command_succeeds(tmp_path: Path) -> None:
    rendered = render_command(["true"], {})
    run_command(rendered)


def test_run_shell_command_writes_log(tmp_path: Path) -> None:
    log_path = tmp_path / "logs" / "step.log"
    rendered = render_command("echo hello world", {})
    run_command(rendered, log_path=log_path)
    content = log_path.read_text()
    assert "$ echo hello world" in content
    assert "hello world" in content


def test_failing_command_raises_command_error(tmp_path: Path) -> None:
    rendered = render_command(["false"], {})
    with pytest.raises(CommandError, match="exit 1"):
        run_command(rendered)


def test_failing_command_error_mentions_log_path(tmp_path: Path) -> None:
    log_path = tmp_path / "step.log"
    rendered = render_command("exit 3", {})
    with pytest.raises(CommandError, match=str(log_path)):
        run_command(rendered, log_path=log_path)
