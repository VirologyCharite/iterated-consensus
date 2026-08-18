from pathlib import Path

import pytest

from iterated_consensus.config import (
    ConfigError,
    InputOverrides,
    InputSpec,
    apply_input_overrides,
    parse_config,
    parse_file_list,
)

MAPPER = """
[[mapper]]
name = "bowtie2"
index_cmd = ["bowtie2-build", "{reference}", "{index_prefix}"]
map_cmd = "bowtie2 -x {index_prefix} -1 {reads_1:,} -2 {reads_2:,} -p {threads} | samtools sort -o {bam}"
"""

CONSENSUS = """
[consensus]
steps = [
    "samtools mpileup -aa -A -d 0 -Q 0 -f {reference} {bam} | ivar consensus -p {consensus_prefix}",
]
output = "{consensus_prefix}.fa"
"""


def test_minimal_valid_config_with_fastq_input() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_1 = ["a_R1.fq", "b_R1.fq"]
reads_2 = ["a_R2.fq", "b_R2.fq"]
reference_fasta = "ref.fasta"
"""
    config = parse_config(text)
    assert config.mappers[0].name == "bowtie2"
    assert config.consensus.output == "{consensus_prefix}.fa"
    assert config.input is not None
    assert config.input.reads_1 == (Path("a_R1.fq"), Path("b_R1.fq"))
    assert config.input.reference_fasta == Path("ref.fasta")
    assert config.input.bam_reads == "ref"


def test_minimal_valid_config_with_bam_input() -> None:
    text = MAPPER + CONSENSUS + """
[input]
bam = "input.bam"
reference_id = "chr2"
bam_reads = "ref+unal"
"""
    config = parse_config(text)
    assert config.input is not None
    assert config.input.bam == Path("input.bam")
    assert config.input.reference_id == "chr2"
    assert config.input.bam_reads == "ref+unal"


def test_run_section_overrides_defaults() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
threads = 8
max_iterations = 5
convergence_identity = 99.5
convergence_streak = 2
"""
    config = parse_config(text)
    assert config.threads == 8
    assert config.max_iterations == 5
    assert config.convergence_identity == 99.5
    assert config.convergence_streak == 2


def test_no_mappers_raises() -> None:
    with pytest.raises(ConfigError, match="at least one"):
        parse_config(CONSENSUS)


def test_duplicate_mapper_names_raise() -> None:
    text = MAPPER + MAPPER + CONSENSUS
    with pytest.raises(ConfigError, match="unique"):
        parse_config(text)


def test_missing_consensus_section_raises() -> None:
    with pytest.raises(ConfigError, match="consensus"):
        parse_config(MAPPER)


def test_consensus_missing_output_raises() -> None:
    text = MAPPER + """
[consensus]
steps = ["samtools mpileup -f {reference} {bam}"]
"""
    with pytest.raises(ConfigError, match="output"):
        parse_config(text)


def test_input_both_fastq_and_bam_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_fasta = "ref.fasta"
bam = "input.bam"
"""
    with pytest.raises(ConfigError, match="not both"):
        parse_config(text)


def test_input_neither_fastq_nor_bam_raises() -> None:
    text = MAPPER + CONSENSUS + "[input]\n"
    with pytest.raises(ConfigError, match="no input given"):
        parse_config(text)


def test_reads_1_reads_2_length_mismatch_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_1 = ["a_R1.fq", "b_R1.fq"]
reads_2 = ["a_R2.fq"]
reference_fasta = "ref.fasta"
"""
    with pytest.raises(ConfigError, match="pair up"):
        parse_config(text)


def test_fastq_input_missing_reference_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
"""
    with pytest.raises(ConfigError, match="starting reference"):
        parse_config(text)


def test_fastq_input_reference_id_alone_is_valid_at_config_level() -> None:
    # Whether "NC_045512.2" actually resolves to something fetchable is a
    # runtime concern (see runner.py) -- config validation only checks that
    # *something* was given.
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"
"""
    config = parse_config(text)
    assert config.input is not None
    assert config.input.reference_fasta is None
    assert config.input.reference_id == "NC_045512.2"


def test_bam_input_with_reference_fasta_is_allowed() -> None:
    text = MAPPER + CONSENSUS + """
[input]
bam = "input.bam"
reference_fasta = "ref.fasta"
"""
    config = parse_config(text)
    assert config.input is not None
    assert config.input.bam == Path("input.bam")
    assert config.input.reference_fasta == Path("ref.fasta")


def test_bam_input_with_reference_fasta_and_reference_id_together() -> None:
    # Not mutually exclusive: reference_id picks a record out of a
    # multi-sequence reference_fasta (and doubles as the BAM contig to use).
    text = MAPPER + CONSENSUS + """
[input]
bam = "input.bam"
reference_fasta = "panel.fasta"
reference_id = "chr2"
"""
    config = parse_config(text)
    assert config.input is not None
    assert config.input.reference_fasta == Path("panel.fasta")
    assert config.input.reference_id == "chr2"


def test_invalid_bam_reads_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
bam = "input.bam"
bam_reads = "everything"
"""
    with pytest.raises(ConfigError, match="bam_reads"):
        parse_config(text)


def test_invalid_convergence_identity_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_fasta = "ref.fasta"

[run]
convergence_identity = 150
"""
    with pytest.raises(ConfigError, match="convergence_identity"):
        parse_config(text)


def test_list_command_step_preserved() -> None:
    config = parse_config(MAPPER + CONSENSUS)
    assert config.mappers[0].index_cmd == ["bowtie2-build", "{reference}", "{index_prefix}"]
    assert isinstance(config.mappers[0].map_cmd, str)


def test_invalid_toml_raises() -> None:
    with pytest.raises(ConfigError, match="invalid TOML"):
        parse_config("this is not [ valid toml")


def test_parse_file_list_splits_on_comma() -> None:
    assert parse_file_list("a.fq,b.fq,c.fq") == (Path("a.fq"), Path("b.fq"), Path("c.fq"))


def test_parse_file_list_single_file() -> None:
    assert parse_file_list("a.fq") == (Path("a.fq"),)


def test_apply_input_overrides_on_blank_base() -> None:
    overrides = InputOverrides(reads_single=(Path("s.fq"),), reference_fasta=Path("ref.fa"))
    merged = apply_input_overrides(None, overrides)
    assert merged.reads_single == (Path("s.fq"),)
    assert merged.reference_fasta == Path("ref.fa")


def test_apply_input_overrides_cli_wins_over_config() -> None:
    base = InputSpec(reads_single=(Path("config_s.fq"),), reference_fasta=Path("config_ref.fa"))
    overrides = InputOverrides(reads_single=(Path("cli_s.fq"),))
    merged = apply_input_overrides(base, overrides)
    assert merged.reads_single == (Path("cli_s.fq"),)
    assert merged.reference_fasta == Path("config_ref.fa")  # untouched


def test_apply_input_overrides_validates_result() -> None:
    base = InputSpec(bam=Path("in.bam"))
    overrides = InputOverrides(reads_single=(Path("s.fq"),), reference_fasta=Path("ref.fa"))
    with pytest.raises(ConfigError, match="not both"):
        apply_input_overrides(base, overrides)
