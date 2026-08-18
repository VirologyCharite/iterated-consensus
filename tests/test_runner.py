import sys
from pathlib import Path

import pysam
import pytest

from iterated_consensus.config import (
    Config,
    ConsensusSpec,
    InputSpec,
    Mapper,
    OutputSpec,
)
from iterated_consensus.reference import ReferenceError
from iterated_consensus.runner import RunnerError, preview, run

FIXTURES = Path(__file__).parent / "fixtures"


def _fake_bam_mapper(name: str = "fake") -> Mapper:
    return Mapper(
        name=name,
        index_cmd=["true"],
        map_cmd=[sys.executable, str(FIXTURES / "make_fake_bam.py"), "{bam}"],
    )


def test_fastq_start_converges_and_writes_outputs(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGTACGTACGTACGTACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "copy_fasta.py"), "{reference}", "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
        threads=1,
        max_iterations=10,
        convergence_identity=100.0,
        convergence_streak=1,
    )

    out_dir = tmp_path / "out"
    result = run(config, out_dir)

    assert result.converged
    assert len(result.iterations) == 2  # iter 0 has no previous consensus to compare
    assert result.iterations[0].identity_to_previous is None
    assert result.iterations[1].identity_to_previous == 100.0
    assert result.iterations[0].consensus_length == 20
    assert result.iterations[1].reads_mapped == 1

    assert (out_dir / "iter_000" / "consensus.fasta").exists()
    assert (out_dir / "iter_001" / "consensus.fasta").exists()
    assert (out_dir / "iter_000" / "logs" / "fake_index.log").exists()
    assert (out_dir / "iter_000" / "logs" / "fake_map.log").exists()
    # The fake mapper doesn't index its own output; the runner does it for us.
    assert (out_dir / "iter_000" / "fake.bam.bai").exists()
    assert (out_dir / "metrics.tsv").exists()
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "index.html").exists()
    assert "Converged" in (out_dir / "index.html").read_text()

    metrics_lines = (out_dir / "metrics.tsv").read_text().splitlines()
    assert len(metrics_lines) == 3  # header + 2 iterations


def test_on_iteration_callback_fires_once_per_executed_iteration(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGTACGTACGTACGTACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "copy_fasta.py"), "{reference}", "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
    )

    seen: list[int] = []
    result = run(config, tmp_path / "out", on_iteration=lambda record: seen.append(record.iteration))

    assert seen == [r.iteration for r in result.iterations]
    assert seen == [0, 1]


def test_bam_start_iteration_0_has_no_mapping_step(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": 50, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        segment = pysam.AlignedSegment()
        segment.query_name = "r1"
        segment.query_sequence = "ACGTACGTAC"
        segment.flag = 0
        segment.reference_id = 0
        segment.reference_start = 0
        segment.mapping_quality = 60
        segment.cigar = [(0, 10)]
        segment.query_qualities = pysam.qualitystring_to_array("I" * 10)
        f.write(segment)

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(bam=bam_path),
        max_iterations=10,
        convergence_identity=100.0,
        convergence_streak=1,
    )

    out_dir = tmp_path / "out"
    result = run(config, out_dir)

    assert result.converged
    assert len(result.iterations) == 2
    # No mapper ran in iteration 0: no index/map logs, only consensus logs.
    assert not (out_dir / "iter_000" / "logs" / "fake_index.log").exists()
    assert not (out_dir / "iter_000" / "logs" / "fake_map.log").exists()
    assert (out_dir / "iter_000" / "logs" / "consensus_step_00.log").exists()
    # Iteration 1 does map, against iteration 0's consensus.
    assert (out_dir / "iter_001" / "logs" / "fake_index.log").exists()
    assert (out_dir / "iter_001" / "logs" / "fake_map.log").exists()


def test_bam_start_never_modifies_the_original_unsorted_input_bam(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 50, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        # coordinate order 10, then 0 -- not sorted, and no SO tag either.
        for name, pos in (("r10", 10), ("r0", 0)):
            segment = pysam.AlignedSegment()
            segment.query_name = name
            segment.query_sequence = "ACGTACGTAC"
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = pos
            segment.mapping_quality = 60
            segment.cigar = [(0, 10)]
            segment.query_qualities = pysam.qualitystring_to_array("I" * 10)
            f.write(segment)
    bytes_before = bam_path.read_bytes()

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(bam=bam_path),
        max_iterations=1,
    )
    run(config, tmp_path / "out")

    assert bam_path.read_bytes() == bytes_before
    assert not Path(str(bam_path) + ".bai").exists()


def test_run_sorts_and_indexes_unsorted_mapper_output(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    unsorted_mapper = Mapper(
        name="unsorted",
        index_cmd=["true"],
        map_cmd=[sys.executable, str(FIXTURES / "make_unsorted_bam.py"), "{bam}"],
    )
    config = Config(
        mappers=(unsorted_mapper,),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
        max_iterations=1,
    )
    out_dir = tmp_path / "out"
    run(config, out_dir)

    bam_path = out_dir / "iter_000" / "unsorted.bam"
    assert (out_dir / "iter_000" / "unsorted.bam.bai").exists()
    with pysam.AlignmentFile(str(bam_path)) as f:
        positions = [r.reference_start for r in f]
    assert positions == sorted(positions)


def test_run_without_input_raises(tmp_path: Path) -> None:
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(steps=(["true"],), output=str(tmp_path / "x.fa")),
        input=None,
    )
    with pytest.raises(RunnerError, match="no input specified"):
        run(config, tmp_path / "out")


def test_preview_fastq_start_renders_iterations_0_and_1(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "copy_fasta.py"), "{reference}", "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
    )
    result = preview(config, tmp_path / "out")
    assert len(result.lines) == 6  # iter 0: index, map, consensus; iter 1: same shape
    assert "true" in result.lines[0]
    assert "make_fake_bam.py" in result.lines[1]
    assert "copy_fasta.py" in result.lines[2]
    assert "true" in result.lines[3]
    assert "make_fake_bam.py" in result.lines[4]
    assert "copy_fasta.py" in result.lines[5]
    # iter 1's mapping output lives under iter_001, but its {reference} is iter_000's consensus.
    assert "iter_001" in result.lines[4]
    assert "iter_000" in result.lines[5] and "consensus.fasta" in result.lines[5]
    assert "always run" in result.note
    assert not (tmp_path / "out").exists()  # nothing actually run, not even out_dir itself


def test_run_section_extra_vars_available_in_mapper_and_consensus_commands(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    mapper = Mapper(
        name="fake",
        index_cmd=["true"],
        map_cmd=[
            sys.executable, str(FIXTURES / "make_fake_bam.py"), "{bam}", "--sort-threads={sort_threads}"
        ],
    )
    config = Config(
        mappers=(mapper,),
        consensus=ConsensusSpec(
            steps=(
                [
                    sys.executable, str(FIXTURES / "copy_fasta.py"), "{reference}", "{consensus_prefix}.fa",
                    "--min-depth={min_depth}",
                ],
            ),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
        extra_vars={"sort_threads": 6, "min_depth": 10},
    )
    result = preview(config, tmp_path / "out")
    assert any("--sort-threads=6" in line for line in result.lines)
    assert any("--min-depth=10" in line for line in result.lines)


def test_run_section_builtin_fields_are_also_placeholders(tmp_path: Path) -> None:
    """threads/threads_reserve/max_iterations/convergence_identity/convergence_streak
    are usable as {name} placeholders too, not just custom [run] variables --
    e.g. for a logging step that records what a run was configured with."""
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=(
                [
                    "echo",
                    "threads={threads}",
                    "threads_reserve={threads_reserve}",
                    "max_iterations={max_iterations}",
                    "convergence_identity={convergence_identity}",
                    "convergence_streak={convergence_streak}",
                ],
            ),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
        threads=6,
        threads_reserve=2,
        max_iterations=15,
        convergence_identity=99.5,
        convergence_streak=3,
    )
    result = preview(config, tmp_path / "out")
    consensus_line = next(line for line in result.lines if "echo" in line)
    assert "threads=6" in consensus_line
    assert "threads_reserve=2" in consensus_line
    assert "max_iterations=15" in consensus_line
    assert "convergence_identity=99.5" in consensus_line
    assert "convergence_streak=3" in consensus_line


def test_preview_bam_start_iteration_0_has_no_mapping_lines(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 50, "SN": "chr1"}]}
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

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=ConsensusSpec(
            steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
            output="{consensus_prefix}.fa",
        ),
        input=InputSpec(bam=bam_path),
    )
    result = preview(config, tmp_path / "out")
    assert len(result.lines) == 4  # iter 0: consensus only; iter 1: index, map, consensus
    assert "write_fixed_fasta.py" in result.lines[0]
    assert "true" in result.lines[1]
    assert "make_fake_bam.py" in result.lines[2]
    assert "write_fixed_fasta.py" in result.lines[3]
    assert "no mapping step" in result.note
    # Nothing actually run: no reads extracted, no BAM sorted/indexed, no
    # out_dir (or reads/ under it) created at all.
    assert not (tmp_path / "out").exists()


def test_preview_bam_start_does_not_extract_reads_or_create_out_dir(tmp_path: Path) -> None:
    """A dry run must not create out_dir, extract FASTQs into reads/, or sort/index
    the input BAM -- it should only print the commands a real run would execute."""
    bam_path = tmp_path / "input.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 50, "SN": "chr1"}]}
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
    original_bam_bytes = bam_path.read_bytes()

    reads_aware_mapper = Mapper(
        name="fake",
        index_cmd=["true"],
        map_cmd=[
            sys.executable, str(FIXTURES / "make_fake_bam.py"), "{bam}", "{reads_1}", "{reads_2}"
        ],
    )
    config = Config(
        mappers=(reads_aware_mapper,),
        consensus=_consensus_ignoring_reference(),
        input=InputSpec(bam=bam_path),
    )
    out_dir = tmp_path / "out"
    result = preview(config, out_dir)

    assert not out_dir.exists()
    assert not bam_path.with_suffix(".bam.bai").exists()
    assert bam_path.read_bytes() == original_bam_bytes  # untouched
    # iter 1's map_cmd still references plausible extracted-reads paths.
    assert any("reads/mate1.fastq.gz" in line and "reads/mate2.fastq.gz" in line for line in result.lines)


def _make_synthetic_bam(bam_path: Path, *, ref_name: str = "chr1", ref_length: int = 50) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": ref_length, "SN": ref_name}]}
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


def _consensus_using_reference() -> ConsensusSpec:
    return ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "copy_fasta.py"), "{reference}", "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )


def _consensus_ignoring_reference() -> ConsensusSpec:
    return ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )


def test_preview_bam_start_uses_explicit_reference(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="chr1", ref_length=20)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">chr1\n" + "ACGT" * 5 + "\n")  # 20bp, matches the BAM header

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path, reference_fasta=reference),
    )
    result = preview(config, tmp_path / "out")
    assert any("reference_initial.fasta" in line for line in result.lines)


def test_preview_bam_start_explicit_reference_unknown_id_raises(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="chr1", ref_length=20)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">wrong_name\n" + "ACGT" * 5 + "\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path, reference_fasta=reference),
    )
    with pytest.raises(ReferenceError, match="not found"):
        preview(config, tmp_path / "out")


def test_preview_bam_start_explicit_reference_length_mismatch_raises(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="chr1", ref_length=20)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">chr1\n" + "ACGT" * 3 + "\n")  # 12bp, BAM header says 20bp

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path, reference_fasta=reference),
    )
    with pytest.raises(RunnerError, match=r"12 bp.*20 bp"):
        preview(config, tmp_path / "out")


def test_preview_bam_start_auto_fetches_accession_looking_contig(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="NC_045512.2", ref_length=20)
    out_dir = tmp_path / "out"
    # Exploit fetch_ncbi_accession's own cache check so no real network call happens.
    cache_dir = out_dir / "reference_cache"
    cache_dir.mkdir(parents=True)
    (cache_dir / "NC_045512.2.fasta").write_text(">NC_045512.2\n" + "ACGT" * 5 + "\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path),
    )
    result = preview(config, out_dir)
    assert any("reference_initial.fasta" in line for line in result.lines)


def test_preview_bam_start_auto_fetch_id_mismatch_raises(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="NC_045512.2", ref_length=20)
    out_dir = tmp_path / "out"
    cache_dir = out_dir / "reference_cache"
    cache_dir.mkdir(parents=True)
    # A cached file that doesn't actually match what we asked for (id differs).
    (cache_dir / "NC_045512.2.fasta").write_text(">NC_045512.1\n" + "ACGT" * 5 + "\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path),
    )
    with pytest.raises(RunnerError, match="does not match BAM reference"):
        preview(config, out_dir)


def test_preview_bam_start_no_reference_ok_if_consensus_does_not_need_it(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="chr1", ref_length=20)  # not accession-like

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_ignoring_reference(),
        input=InputSpec(bam=bam_path),
    )
    result = preview(config, tmp_path / "out")
    assert len(result.lines) == 4  # iter 0: consensus only; iter 1: index, map, consensus


def test_preview_bam_start_no_reference_but_consensus_needs_it_raises(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="chr1", ref_length=20)  # not accession-like

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path),
    )
    with pytest.raises(RunnerError, match="can never succeed"):
        preview(config, tmp_path / "out")


def test_run_symlinks_reference_initial_for_single_record_fastq_start(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fasta"
    reference.write_text(">ref1\nACGTACGTACGTACGTACGT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference),
        max_iterations=1,
    )
    out_dir = tmp_path / "out"
    run(config, out_dir)

    link = out_dir / "reference_initial.fasta"
    assert link.is_symlink()
    assert not link.readlink().is_absolute()
    assert link.read_text() == reference.read_text()


def test_run_symlinks_reference_initial_for_bam_start(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path, ref_name="chr1", ref_length=20)
    reference = tmp_path / "ref.fasta"
    reference.write_text(">chr1\n" + "ACGT" * 5 + "\n")  # 20bp, matches the BAM header

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(bam=bam_path, reference_fasta=reference),
        max_iterations=1,
    )
    out_dir = tmp_path / "out"
    run(config, out_dir)

    link = out_dir / "reference_initial.fasta"
    assert link.is_symlink()
    assert not link.readlink().is_absolute()
    assert link.read_text() == reference.read_text()


def test_run_does_not_symlink_multi_record_reference(tmp_path: Path) -> None:
    reference = tmp_path / "panel.fasta"
    reference.write_text(">ref1\nACGTACGTACGTACGTACGT\n>ref2\nTTTT\n")
    unpaired = tmp_path / "s.fastq"
    unpaired.write_text("@r1\nACGT\n+\nIIII\n")

    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=_consensus_using_reference(),
        input=InputSpec(reads_single=(unpaired,), reference_fasta=reference, reference_id="ref1"),
        max_iterations=1,
    )
    out_dir = tmp_path / "out"
    run(config, out_dir)

    link = out_dir / "reference_initial.fasta"
    assert not link.is_symlink()
    assert link.read_text() == ">ref1\nACGTACGTACGTACGTACGT\n"


def test_resume_continues_past_a_raised_max_iterations(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_sequence_by_iteration.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"

    capped_config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        max_iterations=2,
    )
    first = run(capped_config, out_dir)
    assert not first.converged
    assert [r.iteration for r in first.iterations] == [0, 1, 2]

    index_log = out_dir / "iter_001" / "logs" / "fake_index.log"
    mtime_before_resume = index_log.stat().st_mtime

    resumed_config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        max_iterations=3,
    )
    second = run(resumed_config, out_dir)

    assert second.converged
    assert [r.iteration for r in second.iterations] == [0, 1, 2, 3]
    # iter_001 wasn't redone.
    assert index_log.stat().st_mtime == mtime_before_resume
    assert second.total_elapsed_seconds >= first.total_elapsed_seconds
    assert (out_dir / "index.html").exists()


def test_on_iteration_not_called_for_resumed_records(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_sequence_by_iteration.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"
    capped_config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        max_iterations=2,
    )
    run(capped_config, out_dir)

    resumed_config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        max_iterations=3,
    )
    seen: list[int] = []
    run(resumed_config, out_dir, on_iteration=lambda record: seen.append(record.iteration))
    assert seen == [3]  # not [0, 1, 2, 3] -- those were loaded, not executed, this call


def test_resume_is_a_noop_once_already_converged(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        max_iterations=10,
    )
    first = run(config, out_dir)
    assert first.converged

    second = run(config, out_dir)
    assert second.converged
    assert second.iterations == first.iterations
    assert second.total_elapsed_seconds == first.total_elapsed_seconds


def test_output_section_copies_final_consensus(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"
    final_path = tmp_path / "delivered" / "final.fasta"
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        output=OutputSpec(consensus_fasta=final_path),
        max_iterations=10,
    )
    result = run(config, out_dir)

    assert result.converged
    assert final_path.exists()
    assert not final_path.is_symlink()  # a real, standalone copy
    assert final_path.read_text() == f">seed\n{'ACGTACGTACGTACGTACGT'}\n"


def test_output_section_renames_consensus_id(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"
    final_path = tmp_path / "final.fasta"
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        output=OutputSpec(consensus_fasta=final_path, consensus_id="my-sample"),
        max_iterations=10,
    )
    run(config, out_dir)

    assert final_path.read_text() == f">my-sample\n{'ACGTACGTACGTACGTACGT'}\n"


def test_output_section_written_on_resumed_already_converged_run(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"
    final_path = tmp_path / "final.fasta"
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        output=OutputSpec(consensus_fasta=final_path),
        max_iterations=10,
    )
    first = run(config, out_dir)
    assert first.converged
    final_path.unlink()  # prove the resumed (already-converged) call writes it again

    second = run(config, out_dir)
    assert second.converged
    assert final_path.exists()


def test_no_output_section_does_not_write_anything_extra(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    out_dir = tmp_path / "out"
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
        max_iterations=10,
    )
    run(config, out_dir)  # output is None -- should not raise or write anything odd


def test_resume_with_no_previous_run_behaves_like_fresh_run(tmp_path: Path) -> None:
    bam_path = tmp_path / "input.bam"
    _make_synthetic_bam(bam_path)

    consensus = ConsensusSpec(
        steps=([sys.executable, str(FIXTURES / "write_fixed_fasta.py"), "{consensus_prefix}.fa"],),
        output="{consensus_prefix}.fa",
    )
    config = Config(
        mappers=(_fake_bam_mapper(),),
        consensus=consensus,
        input=InputSpec(bam=bam_path),
    )
    result = run(config, tmp_path / "out")
    assert result.converged
    assert result.iterations[0].iteration == 0
