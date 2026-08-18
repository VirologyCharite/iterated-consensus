import sys
from pathlib import Path

import pysam
import pytest
from typer.testing import CliRunner

from iterated_consensus.cli import app

runner = CliRunner()
FIXTURES = Path(__file__).parent / "fixtures"


def _write_synthetic_bam(bam_path: Path) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 20, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        segment = pysam.AlignedSegment()
        segment.query_name = "r1"
        segment.query_sequence = "ACGT"
        segment.flag = 0
        segment.reference_id = 0
        segment.reference_start = 0
        segment.mapping_quality = 60
        segment.cigar = [(0, 4)]
        segment.query_qualities = pysam.qualitystring_to_array("IIII")
        f.write(segment)


def test_config_template_output_has_no_trailing_blank_line() -> None:
    result = runner.invoke(app, ["config-template", "bowtie2-ivar"])
    assert result.exit_code == 0
    # Exactly one trailing newline (the file's own), not two (which would
    # print as a stray blank line at the end).
    assert result.output.endswith("sort_threads = 1\n")
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


def test_run_prints_final_consensus_path_when_output_configured(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    final_path = tmp_path / "final.fasta"
    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(f"""
[[mapper]]
name = "fake"
index_cmd = ["true"]
map_cmd = ["{sys.executable}", "{FIXTURES / "make_fake_bam.py"}", "{{bam}}"]

[consensus]
steps = [["{sys.executable}", "{FIXTURES / "write_fixed_fasta.py"}", "{{consensus_prefix}}.fa"]]
output = "{{consensus_prefix}}.fa"

[input]
bam = "{bam_path}"

[output]
consensus_fasta = "{final_path}"
consensus_id = "final-id"

[run]
max_iterations = 1
""")
    result = runner.invoke(
        app, ["run", "--config", str(config_path), "--output-dir", str(tmp_path / "out")]
    )
    assert result.exit_code == 0, result.output
    assert f"Final consensus: {final_path}" in result.output
    assert final_path.read_text().startswith(">final-id\n")


def _minimal_config_text(bam_path: Path, *, output_section: str = "", run_extra: str = "") -> str:
    return f"""
[[mapper]]
name = "fake"
index_cmd = ["true"]
map_cmd = ["{sys.executable}", "{FIXTURES / "make_fake_bam.py"}", "{{bam}}"]

[consensus]
steps = [["{sys.executable}", "{FIXTURES / "write_fixed_fasta.py"}", "{{consensus_prefix}}.fa"]]
output = "{{consensus_prefix}}.fa"

[input]
bam = "{bam_path}"

{output_section}

[run]
max_iterations = 1
{run_extra}
"""


def test_run_output_dir_placeholder_places_consensus_fasta_under_out_dir(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    out_dir = tmp_path / "out"
    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(
        _minimal_config_text(bam_path, output_section='[output]\nconsensus_fasta = "{output_dir}/final.fasta"')
    )
    result = runner.invoke(app, ["run", "--config", str(config_path), "--output-dir", str(out_dir)])
    assert result.exit_code == 0, result.output
    expected_path = out_dir / "final.fasta"
    assert f"Final consensus: {expected_path}" in result.output
    assert expected_path.exists()


def test_run_relative_consensus_fasta_is_relative_to_cwd_not_out_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    out_dir = tmp_path / "out"
    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(
        _minimal_config_text(bam_path, output_section='[output]\nconsensus_fasta = "final.fasta"')
    )
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)

    result = runner.invoke(app, ["run", "--config", str(config_path), "--output-dir", str(out_dir)])
    assert result.exit_code == 0, result.output
    assert "Final consensus: final.fasta" in result.output  # left relative, not resolved against out_dir
    assert (cwd / "final.fasta").exists()
    assert not (out_dir / "final.fasta").exists()


def test_run_output_dir_falls_back_to_run_section_config(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    out_dir = tmp_path / "out"
    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(_minimal_config_text(bam_path, run_extra=f'output_dir = "{out_dir}"'))

    result = runner.invoke(app, ["run", "--config", str(config_path)])  # no --output-dir given
    assert result.exit_code == 0, result.output
    assert (out_dir / "iter_000" / "consensus.fasta").exists()


def test_run_cli_output_dir_overrides_run_section_config(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    cli_out_dir = tmp_path / "cli-out"
    config_out_dir = tmp_path / "config-out"
    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(_minimal_config_text(bam_path, run_extra=f'output_dir = "{config_out_dir}"'))

    result = runner.invoke(
        app, ["run", "--config", str(config_path), "--output-dir", str(cli_out_dir)]
    )
    assert result.exit_code == 0, result.output
    assert (cli_out_dir / "iter_000" / "consensus.fasta").exists()
    assert not config_out_dir.exists()


def test_run_missing_output_dir_errors_clearly(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(_minimal_config_text(bam_path))  # no [run].output_dir either

    result = runner.invoke(app, ["run", "--config", str(config_path)])
    assert result.exit_code == 1
    assert "output directory" in result.output
    assert "--output-dir" in result.output


def test_progress_prints_header_once_not_per_line(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(f"""
[[mapper]]
name = "fake"
index_cmd = ["true"]
map_cmd = ["{sys.executable}", "{FIXTURES / "make_fake_bam.py"}", "{{bam}}"]

[consensus]
steps = [["{sys.executable}", "{FIXTURES / "copy_fasta.py"}", "{{reference}}", "{{consensus_prefix}}.fa"]]
output = "{{consensus_prefix}}.fa"

[input]
reads_single = ["{unpaired}"]
reference_fasta = "{reference}"

[run]
max_iterations = 5
""")
    result = runner.invoke(
        app, ["run", "--config", str(config_path), "--output-dir", str(tmp_path / "out"), "--progress"]
    )
    assert result.exit_code == 0, result.output
    lines = result.output.splitlines()
    header_lines = [line for line in lines if line.startswith("iter")]
    assert len(header_lines) == 1  # column names appear once, not per iteration
    assert "reads_mapped" in header_lines[0]
    assert "consensus_length" in header_lines[0]
    assert "identity_to_previous" in header_lines[0]
    assert "elapsed" in header_lines[0]
    assert "consensus_md5" in header_lines[0]
    # iteration rows don't repeat the column names
    iteration_lines = [line for line in lines if line.strip().startswith(("0", "1"))]
    assert iteration_lines
    assert not any("reads_mapped=" in line for line in iteration_lines)
    # each iteration row ends with a 32-char hex MD5
    for line in iteration_lines:
        md5 = line.split()[-1]
        assert len(md5) == 32
        int(md5, 16)  # raises ValueError if not valid hex


def test_run_prints_cycle_detected_status_and_summary(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _write_synthetic_bam(bam_path)

    config_path = tmp_path / "pipelines.toml"
    config_path.write_text(f"""
[[mapper]]
name = "fake"
index_cmd = ["true"]
map_cmd = ["{sys.executable}", "{FIXTURES / "make_fake_bam.py"}", "{{bam}}"]

[consensus]
steps = [["{sys.executable}", "{FIXTURES / "write_alternating_sequence.py"}", "{{consensus_prefix}}.fa"]]
output = "{{consensus_prefix}}.fa"

[input]
bam = "{bam_path}"

[run]
max_iterations = 20
""")
    result = runner.invoke(
        app, ["run", "--config", str(config_path), "--output-dir", str(tmp_path / "out")]
    )
    assert result.exit_code == 0, result.output
    assert "Cycle detected (period 2) after 3 iteration(s)" in result.output
    assert "Iteration 2's consensus matches iteration 0's" in result.output
    assert "Using iteration 0 as the final consensus" in result.output
