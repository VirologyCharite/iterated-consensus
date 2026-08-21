from pathlib import Path

import pytest

from iterated_consensus.config import (
    ConfigError,
    InputOverrides,
    InputSpec,
    apply_input_overrides,
    available_cpu_count,
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


def test_run_section_arbitrary_variables_become_extra_vars() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
sort_threads = 6
min_depth = 10
sample_name = "patient-42"
strict = true
"""
    config = parse_config(text)
    assert config.extra_vars == {
        "sort_threads": 6,
        "min_depth": 10,
        "sample_name": "patient-42",
        "strict": True,
    }
    # known [run] keys never leak into extra_vars
    assert "threads" not in config.extra_vars
    assert "max_iterations" not in config.extra_vars


def test_run_section_extra_var_colliding_with_reserved_placeholder_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
reference = "oops"
"""
    with pytest.raises(ConfigError, match="collide"):
        parse_config(text)


def test_run_section_extra_var_bad_type_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
bad = [1, 2, 3]
"""
    with pytest.raises(ConfigError, match="string, number, or boolean"):
        parse_config(text)


def test_run_section_threads_auto_uses_available_cpu_count() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
threads = "auto"
"""
    config = parse_config(text)
    assert config.threads == available_cpu_count()


def test_run_section_threads_auto_with_reserve() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
threads = "auto"
threads_reserve = 2
"""
    config = parse_config(text)
    assert config.threads == max(1, available_cpu_count() - 2)


def test_run_section_threads_reserve_never_drops_below_one() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
threads = "auto"
threads_reserve = 999999
"""
    config = parse_config(text)
    assert config.threads == 1


def test_run_section_threads_reserve_without_auto_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
threads = 4
threads_reserve = 2
"""
    with pytest.raises(ConfigError, match="threads_reserve"):
        parse_config(text)


def test_run_section_threads_bad_value_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
threads = "fast"
"""
    with pytest.raises(ConfigError, match="threads"):
        parse_config(text)


def test_run_section_output_dir_parsed() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
output_dir = "results/run1"
"""
    config = parse_config(text)
    assert config.output_dir == Path("results/run1")


def test_no_run_section_output_dir_is_none() -> None:
    config = parse_config(MAPPER + CONSENSUS)
    assert config.output_dir is None


def test_run_section_output_dir_bad_type_raises() -> None:
    text = MAPPER + CONSENSUS + """
[run]
output_dir = 5
"""
    with pytest.raises(ConfigError, match="output_dir"):
        parse_config(text)


def test_run_section_output_dir_not_leaked_into_extra_vars() -> None:
    text = MAPPER + CONSENSUS + """
[run]
output_dir = "results/run1"
"""
    config = parse_config(text)
    assert "output_dir" not in config.extra_vars


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


def test_output_section_parsed() -> None:
    text = MAPPER + CONSENSUS + """
[output]
consensus_fasta = "final/result.fasta"
consensus_id = "my-sample"
"""
    config = parse_config(text)
    assert config.output is not None
    assert config.output.consensus_fasta == "final/result.fasta"
    assert config.output.consensus_id == "my-sample"


def test_output_section_consensus_id_alone_is_optional() -> None:
    text = MAPPER + CONSENSUS + """
[output]
consensus_fasta = "final/result.fasta"
"""
    config = parse_config(text)
    assert config.output is not None
    assert config.output.consensus_fasta == "final/result.fasta"
    assert config.output.consensus_id is None


def test_output_section_consensus_id_without_consensus_fasta_raises() -> None:
    text = MAPPER + CONSENSUS + """
[output]
consensus_id = "my-sample"
"""
    with pytest.raises(ConfigError, match="consensus_id needs consensus_fasta"):
        parse_config(text)


def test_no_output_section_is_fine() -> None:
    config = parse_config(MAPPER + CONSENSUS)
    assert config.output is None


def test_output_section_commands_parsed() -> None:
    text = MAPPER + CONSENSUS + """
[output]
consensus_fasta = "final/result.fasta"
commands = [
    ["samtools", "faidx", "{consensus_fasta}"],
    "gzip -k {consensus_fasta}",
]
"""
    config = parse_config(text)
    assert config.output is not None
    assert config.output.commands == (
        ["samtools", "faidx", "{consensus_fasta}"],
        "gzip -k {consensus_fasta}",
    )


def test_output_section_commands_without_consensus_fasta_raises() -> None:
    text = MAPPER + CONSENSUS + """
[output]
commands = [["samtools", "faidx", "{consensus_fasta}"]]
"""
    with pytest.raises(ConfigError, match="commands needs consensus_fasta"):
        parse_config(text)


def test_output_section_commands_default_empty() -> None:
    text = MAPPER + CONSENSUS + """
[output]
consensus_fasta = "final/result.fasta"
"""
    config = parse_config(text)
    assert config.output is not None
    assert config.output.commands == ()


def test_output_section_commands_not_a_list_raises() -> None:
    text = MAPPER + CONSENSUS + """
[output]
consensus_fasta = "final/result.fasta"
commands = "not a list"
"""
    with pytest.raises(ConfigError, match="commands must be a list"):
        parse_config(text)


def test_output_section_final_reference_fields_parsed() -> None:
    text = MAPPER + CONSENSUS + """
[output]
final_reference_fasta = "final/reference.fasta"
final_reference_bam = "final/reference.bam"
"""
    config = parse_config(text)
    assert config.output is not None
    assert config.output.final_reference_fasta == "final/reference.fasta"
    assert config.output.final_reference_bam == "final/reference.bam"


def test_output_section_final_reference_fields_independent_of_consensus_fasta() -> None:
    text = MAPPER + CONSENSUS + """
[output]
final_reference_fasta = "final/reference.fasta"
"""
    config = parse_config(text)
    assert config.output is not None
    assert config.output.consensus_fasta is None
    assert config.output.final_reference_fasta == "final/reference.fasta"


def test_output_section_final_reference_fasta_not_a_string_raises() -> None:
    text = MAPPER + CONSENSUS + """
[output]
final_reference_fasta = 123
"""
    with pytest.raises(ConfigError, match="final_reference_fasta must be a string"):
        parse_config(text)


def test_output_section_final_reference_bam_not_a_string_raises() -> None:
    text = MAPPER + CONSENSUS + """
[output]
final_reference_bam = 123
"""
    with pytest.raises(ConfigError, match="final_reference_bam must be a string"):
        parse_config(text)


def test_run_section_extra_var_colliding_with_final_reference_fasta_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
final_reference_fasta = "oops"
"""
    with pytest.raises(ConfigError, match="collide"):
        parse_config(text)


def test_run_section_extra_var_colliding_with_final_reference_bam_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
final_reference_bam = "oops"
"""
    with pytest.raises(ConfigError, match="collide"):
        parse_config(text)


def test_run_section_extra_var_colliding_with_consensus_fasta_raises() -> None:
    text = MAPPER + CONSENSUS + """
[input]
reads_single = ["s.fq"]
reference_id = "NC_045512.2"

[run]
consensus_fasta = "oops"
"""
    with pytest.raises(ConfigError, match="collide"):
        parse_config(text)


def test_list_command_step_preserved() -> None:
    config = parse_config(MAPPER + CONSENSUS)
    assert config.mappers[0].index_cmd == ["bowtie2-build", "{reference}", "{index_prefix}"]
    assert isinstance(config.mappers[0].map_cmd, str)


def test_no_tool_versions_is_fine() -> None:
    config = parse_config(MAPPER + CONSENSUS)
    assert config.mappers[0].tool_versions == {}
    assert config.consensus.tool_versions == {}


def test_mapper_tool_versions_parsed() -> None:
    text = MAPPER + """
[mapper.tool-versions]
bowtie2 = "bowtie2 --version | head -n 1 | cut -f3 -d' '"
samtools = ["samtools", "--version"]
""" + CONSENSUS
    config = parse_config(text)
    assert config.mappers[0].tool_versions == {
        "bowtie2": "bowtie2 --version | head -n 1 | cut -f3 -d' '",
        "samtools": ["samtools", "--version"],
    }


def test_consensus_tool_versions_parsed() -> None:
    text = MAPPER + CONSENSUS + """
[consensus.tool-versions]
ivar = "ivar version | head -n 1 | cut -f3 -d' '"
"""
    config = parse_config(text)
    assert config.consensus.tool_versions == {"ivar": "ivar version | head -n 1 | cut -f3 -d' '"}


def test_tool_versions_attach_to_correct_mapper_in_array() -> None:
    text = MAPPER + """
[mapper.tool-versions]
bowtie2 = "echo bowtie2-version"

[[mapper]]
name = "bwa"
index_cmd = ["bwa", "index"]
map_cmd = "bwa mem"
""" + CONSENSUS
    config = parse_config(text)
    assert config.mappers[0].name == "bowtie2"
    assert config.mappers[0].tool_versions == {"bowtie2": "echo bowtie2-version"}
    assert config.mappers[1].name == "bwa"
    assert config.mappers[1].tool_versions == {}


def test_mapper_tool_versions_not_a_table_raises() -> None:
    text = """
[[mapper]]
name = "bowtie2"
index_cmd = ["bowtie2-build"]
map_cmd = "bowtie2"
"tool-versions" = "not a table"
""" + CONSENSUS
    with pytest.raises(ConfigError, match="must be a table"):
        parse_config(text)


def test_tool_versions_bad_command_type_raises() -> None:
    text = MAPPER + CONSENSUS + """
[consensus.tool-versions]
ivar = 5
"""
    with pytest.raises(ConfigError, match="string or a list of strings"):
        parse_config(text)


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
