import gzip
from pathlib import Path

import pysam
import pytest

from iterated_consensus.bam import (
    BamError,
    count_mapped_reads,
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
