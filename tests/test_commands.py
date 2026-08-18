from pathlib import Path

import pytest

from iterated_consensus.commands import (
    CommandError,
    render_command,
    run_command,
    run_logged_command,
    run_tool_version_command,
)


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


def test_run_logged_command_captures_display_elapsed_and_log(tmp_path: Path) -> None:
    rendered = render_command("echo hello world", {})
    result = run_logged_command("greeting", rendered, tmp_path / "logs")
    assert result.name == "greeting"
    assert result.display == "echo hello world"
    assert result.elapsed_seconds >= 0
    assert "$ echo hello world" in result.log
    assert "hello world" in result.log
    assert (tmp_path / "logs" / "greeting.log").exists()


def test_run_logged_command_failure_raises_command_error(tmp_path: Path) -> None:
    rendered = render_command(["false"], {})
    with pytest.raises(CommandError, match="exit 1"):
        run_logged_command("fails", rendered, tmp_path / "logs")


def test_run_tool_version_command_captures_stdout() -> None:
    version = run_tool_version_command("echo 'MyTool 1.2.3'", {})
    assert version == "MyTool 1.2.3"


def test_run_tool_version_command_supports_placeholders() -> None:
    version = run_tool_version_command("echo tool-{name}", {"name": "x"})
    assert version == "tool-x"


def test_run_tool_version_command_strips_and_keeps_multiple_lines() -> None:
    version = run_tool_version_command("printf '  v1\\nbuild abc  \\n'", {})
    assert version == "v1\nbuild abc"


def test_run_tool_version_command_failure_raises() -> None:
    with pytest.raises(CommandError, match="exit 1"):
        run_tool_version_command(["false"], {})


def test_run_tool_version_command_ignores_stderr() -> None:
    version = run_tool_version_command("echo out; echo err >&2", {})
    assert version == "out"
