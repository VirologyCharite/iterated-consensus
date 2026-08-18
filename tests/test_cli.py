import sys
from pathlib import Path

import pysam
from typer.testing import CliRunner

from iterated_consensus.cli import app

runner = CliRunner()
FIXTURES = Path(__file__).parent / "fixtures"


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
        app, ["run", "--config", str(config_path), "--out-dir", str(tmp_path / "out")]
    )
    assert result.exit_code == 0, result.output
    assert f"Final consensus: {final_path}" in result.output
    assert final_path.read_text().startswith(">final-id\n")
