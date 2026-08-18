import gzip
from pathlib import Path

import pysam
import pytest

from iterated_consensus.bam import (
    BamError,
    count_mapped_reads,
    ensure_indexed,
    ensure_indexed_readonly,
    extract_fastq,
    resolve_reference_name,
)


def _segment(
    name: str,
    seq: str,
    *,
    flag: int,
    reference_id: int,
    reference_start: int,
    next_reference_id: int = -1,
    next_reference_start: int = -1,
) -> pysam.AlignedSegment:
    s = pysam.AlignedSegment()
    s.query_name = name
    s.query_sequence = seq
    s.flag = flag
    s.reference_id = reference_id
    s.reference_start = reference_start
    s.next_reference_id = next_reference_id
    s.next_reference_start = next_reference_start
    s.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    if not flag & 4:
        s.mapping_quality = 60
        s.cigar = [(0, len(seq))]
    return s


@pytest.fixture
def two_ref_bam(tmp_path: Path) -> Path:
    bam_path = tmp_path / "input.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": 50, "SN": "chr1"}, {"LN": 50, "SN": "chr2"}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        # a properly paired pair on chr1
        f.write(_segment("pairA", "ACGTACGTAC", flag=99, reference_id=0, reference_start=0,
                          next_reference_id=0, next_reference_start=0))
        f.write(_segment("pairA", "GTACGTACGT", flag=147, reference_id=0, reference_start=0,
                          next_reference_id=0, next_reference_start=0))
        # a properly paired pair on chr2
        f.write(_segment("pairB", "TTTTAAAAGG", flag=99, reference_id=1, reference_start=5,
                          next_reference_id=1, next_reference_start=5))
        f.write(_segment("pairB", "CCGGTTAACC", flag=147, reference_id=1, reference_start=5,
                          next_reference_id=1, next_reference_start=5))
        # an unmapped single read
        f.write(_segment("single_unmapped", "AAAA", flag=4, reference_id=-1, reference_start=-1))
    return bam_path


def test_resolve_reference_name_requires_choice_for_multi_ref(two_ref_bam: Path) -> None:
    with pytest.raises(BamError, match="multiple references"):
        resolve_reference_name(two_ref_bam, None)


def test_resolve_reference_name_accepts_valid_choice(two_ref_bam: Path) -> None:
    assert resolve_reference_name(two_ref_bam, "chr2") == "chr2"


def test_resolve_reference_name_rejects_unknown_choice(two_ref_bam: Path) -> None:
    with pytest.raises(BamError, match="not found"):
        resolve_reference_name(two_ref_bam, "chr9")


def test_resolve_reference_name_auto_picks_single_ref(tmp_path: Path) -> None:
    bam_path = tmp_path / "single.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 50, "SN": "only"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        f.write(_segment("r", "ACGT", flag=0, reference_id=0, reference_start=0))
    assert resolve_reference_name(bam_path, None) == "only"


def test_count_mapped_reads_total_and_per_reference(two_ref_bam: Path) -> None:
    assert count_mapped_reads(two_ref_bam) == 4
    assert count_mapped_reads(two_ref_bam, reference_name="chr1") == 2
    assert count_mapped_reads(two_ref_bam, reference_name="chr2") == 2


def test_extract_fastq_ref_mode_excludes_other_contig_and_unmapped(
    two_ref_bam: Path, tmp_path: Path
) -> None:
    result = extract_fastq(two_ref_bam, reference_name="chr1", mode="ref", out_dir=tmp_path / "out")
    assert set(result) == {"mate1", "mate2"}
    m1 = gzip.decompress(result["mate1"].read_bytes()).decode()
    m2 = gzip.decompress(result["mate2"].read_bytes()).decode()
    assert "pairA" in m1
    assert "pairA" in m2
    assert "pairB" not in m1 and "pairB" not in m2


def test_extract_fastq_ref_plus_unal_includes_unmapped(two_ref_bam: Path, tmp_path: Path) -> None:
    result = extract_fastq(
        two_ref_bam, reference_name="chr1", mode="ref+unal", out_dir=tmp_path / "out"
    )
    assert "unpaired" in result
    unpaired = gzip.decompress(result["unpaired"].read_bytes()).decode()
    assert "single_unmapped" in unpaired
    m1 = gzip.decompress(result["mate1"].read_bytes()).decode()
    assert "pairB" not in m1


def test_extract_fastq_all_mode_includes_everything(two_ref_bam: Path, tmp_path: Path) -> None:
    result = extract_fastq(two_ref_bam, reference_name="chr1", mode="all", out_dir=tmp_path / "out")
    m1 = gzip.decompress(result["mate1"].read_bytes()).decode()
    assert "pairA" in m1 and "pairB" in m1
    unpaired = gzip.decompress(result["unpaired"].read_bytes()).decode()
    assert "single_unmapped" in unpaired


def test_extract_fastq_unknown_mode_raises(two_ref_bam: Path, tmp_path: Path) -> None:
    with pytest.raises(BamError, match="unknown bam_reads mode"):
        extract_fastq(two_ref_bam, reference_name="chr1", mode="bogus", out_dir=tmp_path / "out")


def test_extract_fastq_is_idempotent(two_ref_bam: Path, tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    first = extract_fastq(two_ref_bam, reference_name="chr1", mode="ref", out_dir=out_dir)

    two_ref_bam.unlink()  # prove the second call doesn't need the source bam at all

    second = extract_fastq(two_ref_bam, reference_name="chr1", mode="ref", out_dir=out_dir)
    assert second == first


def test_ensure_indexed_creates_bai_when_missing(two_ref_bam: Path) -> None:
    bai_path = Path(str(two_ref_bam) + ".bai")
    assert not bai_path.exists()
    ensure_indexed(two_ref_bam)
    assert bai_path.exists()


def test_ensure_indexed_skips_sorting_when_header_says_sorted(two_ref_bam: Path) -> None:
    # two_ref_bam's header already declares SO=coordinate (and really is
    # sorted) -- indexing should use that directly rather than sorting
    # first, which would rewrite the file even though it doesn't need to.
    bytes_before = two_ref_bam.read_bytes()
    ensure_indexed(two_ref_bam)
    assert two_ref_bam.read_bytes() == bytes_before


def test_ensure_indexed_falls_back_to_sorting_if_header_lied(tmp_path: Path) -> None:
    bam_path = tmp_path / "mislabeled.bam"
    # Header claims coordinate-sorted, but the reads aren't actually.
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": 20, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        f.write(_segment("r10", "ACGT", flag=0, reference_id=0, reference_start=10))
        f.write(_segment("r0", "ACGT", flag=0, reference_id=0, reference_start=0))

    ensure_indexed(bam_path)

    assert Path(str(bam_path) + ".bai").exists()
    with pysam.AlignmentFile(str(bam_path)) as f:
        positions = [r.reference_start for r in f]
    assert positions == [0, 10]  # actually sorted despite the misleading header


def test_ensure_indexed_is_a_noop_if_already_indexed(two_ref_bam: Path) -> None:
    bai_path = Path(str(two_ref_bam) + ".bai")
    ensure_indexed(two_ref_bam)
    mtime_before = bai_path.stat().st_mtime

    two_ref_bam.unlink()  # prove the second call doesn't need the source bam at all
    ensure_indexed(two_ref_bam)

    assert bai_path.stat().st_mtime == mtime_before


def test_ensure_indexed_sorts_unsorted_bam_then_indexes(tmp_path: Path) -> None:
    bam_path = tmp_path / "unsorted.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 20, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        # coordinate order 10, then 0 -- not sorted
        f.write(_segment("r10", "ACGT", flag=0, reference_id=0, reference_start=10))
        f.write(_segment("r0", "ACGT", flag=0, reference_id=0, reference_start=0))

    ensure_indexed(bam_path)

    assert Path(str(bam_path) + ".bai").exists()
    with pysam.AlignmentFile(str(bam_path)) as f:
        positions = [r.reference_start for r in f]
    assert positions == [0, 10]  # now actually sorted


def test_ensure_indexed_unreadable_bam_raises_bam_error(tmp_path: Path) -> None:
    bam_path = tmp_path / "garbage.bam"
    bam_path.write_bytes(b"this is not a bam file")

    with pytest.raises(BamError, match="could not sort"):
        ensure_indexed(bam_path)


def test_ensure_indexed_readonly_already_indexed_returns_original_untouched(
    two_ref_bam: Path, tmp_path: Path
) -> None:
    ensure_indexed(two_ref_bam)  # pre-index it, as if from an earlier run
    bytes_before = two_ref_bam.read_bytes()

    work_dir = tmp_path / "work"
    result = ensure_indexed_readonly(two_ref_bam, work_dir)

    assert result == two_ref_bam
    assert two_ref_bam.read_bytes() == bytes_before
    assert not work_dir.exists()  # never even created -- nothing needed it


def test_ensure_indexed_readonly_sorted_indexes_beside_original(
    two_ref_bam: Path, tmp_path: Path
) -> None:
    # two_ref_bam's header already declares SO=coordinate and really is sorted.
    bytes_before = two_ref_bam.read_bytes()
    work_dir = tmp_path / "work"

    result = ensure_indexed_readonly(two_ref_bam, work_dir)

    assert result == two_ref_bam
    assert two_ref_bam.read_bytes() == bytes_before  # content untouched
    assert Path(str(two_ref_bam) + ".bai").exists()  # index added beside it
    assert not work_dir.exists()  # no copy was needed


def test_ensure_indexed_readonly_unsorted_never_modifies_original(tmp_path: Path) -> None:
    bam_path = tmp_path / "unsorted.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 20, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        f.write(_segment("r10", "ACGT", flag=0, reference_id=0, reference_start=10))
        f.write(_segment("r0", "ACGT", flag=0, reference_id=0, reference_start=0))
    bytes_before = bam_path.read_bytes()
    work_dir = tmp_path / "work"

    result = ensure_indexed_readonly(bam_path, work_dir)

    # The original is completely untouched: same bytes, no sidecar index.
    assert bam_path.read_bytes() == bytes_before
    assert not Path(str(bam_path) + ".bai").exists()

    # A separate, sorted, indexed copy was made instead.
    assert result != bam_path
    assert result.parent == work_dir
    assert Path(str(result) + ".bai").exists()
    with pysam.AlignmentFile(str(result)) as f:
        positions = [r.reference_start for r in f]
    assert positions == [0, 10]


def test_ensure_indexed_readonly_is_idempotent(tmp_path: Path) -> None:
    bam_path = tmp_path / "unsorted.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 20, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as f:
        f.write(_segment("r10", "ACGT", flag=0, reference_id=0, reference_start=10))
        f.write(_segment("r0", "ACGT", flag=0, reference_id=0, reference_start=0))
    work_dir = tmp_path / "work"

    first = ensure_indexed_readonly(bam_path, work_dir)
    mtime_before = first.stat().st_mtime
    second = ensure_indexed_readonly(bam_path, work_dir)

    assert second == first
    assert second.stat().st_mtime == mtime_before  # not re-sorted
