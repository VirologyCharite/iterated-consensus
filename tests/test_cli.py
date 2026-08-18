from typer.testing import CliRunner

from iterated_consensus.cli import app

runner = CliRunner()


def test_config_template_output_has_no_trailing_blank_line() -> None:
    result = runner.invoke(app, ["config-template", "bowtie2-ivar"])
    assert result.exit_code == 0
    # Exactly one trailing newline (the file's own), not two (which would
    # print as a stray blank line at the end).
    assert result.output.endswith("convergence_streak = 1\n")
    assert not result.output.endswith("\n\n")


def test_config_template_list_names_presets() -> None:
    result = runner.invoke(app, ["config-template"])
    assert result.exit_code == 0
    assert "bowtie2-ivar" in result.output
    assert "bwa-samtools" in result.output


def test_config_template_unknown_preset_errors() -> None:
    result = runner.invoke(app, ["config-template", "does-not-exist"])
    assert result.exit_code == 1
    assert "unknown preset" in result.output
