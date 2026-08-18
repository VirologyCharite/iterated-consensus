import sys
from pathlib import Path

import pysam
import pytest

from iterated_consensus.config import Config, ConsensusSpec, InputSpec, Mapper
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
        input=InputSpec(unpaired=(unpaired,), reference_fasta=reference),
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
        input=InputSpec(unpaired=(unpaired,), reference_fasta=reference),
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
        input=InputSpec(unpaired=(unpaired,), reference_fasta=reference),
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
        input=InputSpec(unpaired=(unpaired,), reference_fasta=reference),
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
    assert not (tmp_path / "out" / "iter_000").exists()  # nothing actually run
    assert not (tmp_path / "out" / "iter_001").exists()


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
