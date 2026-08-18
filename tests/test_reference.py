from pathlib import Path

import pytest

from iterated_consensus.reference import (
    ReferenceError,
    load_reference,
    looks_like_ncbi_accession,
    parse_fasta,
    resolve_reference,
    symlink_reference,
    write_fasta,
)


def test_parse_single_record() -> None:
    text = ">chr1 some description\nACGT\nACGT\n"
    records = parse_fasta(text)
    assert len(records) == 1
    assert records[0].id == "chr1"
    assert records[0].description == "some description"
    assert records[0].sequence == "ACGTACGT"


def test_parse_multi_record() -> None:
    text = ">a\nAAAA\n>b desc\nCCCC\nGGGG\n"
    records = parse_fasta(text)
    assert [r.id for r in records] == ["a", "b"]
    assert records[1].sequence == "CCCCGGGG"


def test_parse_empty_text() -> None:
    assert parse_fasta("") == []


def test_load_reference_single_record(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">only\nACGT\n")
    record = load_reference(fasta)
    assert record.id == "only"


def test_load_reference_multi_record_requires_id(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">a\nAAAA\n>b\nCCCC\n")
    with pytest.raises(ReferenceError, match="multiple sequences"):
        load_reference(fasta)


def test_load_reference_multi_record_with_id(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">a\nAAAA\n>b\nCCCC\n")
    record = load_reference(fasta, reference_id="b")
    assert record.sequence == "CCCC"


def test_load_reference_unknown_id_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">a\nAAAA\n")
    with pytest.raises(ReferenceError, match="not found"):
        load_reference(fasta, reference_id="nope")


def test_load_reference_no_sequences_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text("")
    with pytest.raises(ReferenceError, match="no sequences"):
        load_reference(fasta)


def test_write_fasta_wraps_lines(tmp_path: Path) -> None:
    out = tmp_path / "out.fasta"
    write_fasta(out, "seq1", "A" * 150, line_width=70)
    lines = out.read_text().splitlines()
    assert lines[0] == ">seq1"
    assert lines[1] == "A" * 70
    assert lines[2] == "A" * 70
    assert lines[3] == "A" * 10


def test_write_fasta_creates_parent_dirs(tmp_path: Path) -> None:
    out = tmp_path / "nested" / "dir" / "out.fasta"
    write_fasta(out, "seq1", "ACGT")
    assert out.read_text() == ">seq1\nACGT\n"


def test_write_then_load_round_trip(tmp_path: Path) -> None:
    out = tmp_path / "out.fasta"
    write_fasta(out, "seq1", "ACGTACGTACGT")
    record = load_reference(out)
    assert record.id == "seq1"
    assert record.sequence == "ACGTACGTACGT"


@pytest.mark.parametrize(
    "name",
    [
        "NC_045512.2",  # RefSeq, versioned
        "NC_045512",  # RefSeq, unversioned
        "NZ_CP012345.1",
        "MN908947.3",  # GenBank
        "U12345.1",
        "mn908947.3",  # lowercase should still count
    ],
)
def test_looks_like_ncbi_accession_true(name: str) -> None:
    assert looks_like_ncbi_accession(name)


@pytest.mark.parametrize(
    "name",
    [
        "chr1",
        "chr2",
        "MT",
        "1",
        "scaffold_12",
        "",
    ],
)
def test_looks_like_ncbi_accession_false(name: str) -> None:
    assert not looks_like_ncbi_accession(name)


def test_resolve_reference_single_record_file_is_whole_file(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">only\nACGT\n")
    resolved = resolve_reference(reference_id=None, reference_fasta=fasta, cache_dir=tmp_path / "cache")
    assert resolved is not None
    assert resolved.record.id == "only"
    assert resolved.source_path == fasta
    assert resolved.is_whole_file


def test_resolve_reference_multi_record_file_is_not_whole_file(tmp_path: Path) -> None:
    fasta = tmp_path / "panel.fasta"
    fasta.write_text(">a\nAAAA\n>b\nCCCC\n")
    resolved = resolve_reference(reference_id="b", reference_fasta=fasta, cache_dir=tmp_path / "cache")
    assert resolved is not None
    assert resolved.record.id == "b"
    assert resolved.source_path == fasta
    assert not resolved.is_whole_file


def test_resolve_reference_none_when_nothing_given(tmp_path: Path) -> None:
    resolved = resolve_reference(reference_id=None, reference_fasta=None, cache_dir=tmp_path / "cache")
    assert resolved is None


def test_resolve_reference_id_not_accession_shaped_returns_none(tmp_path: Path) -> None:
    resolved = resolve_reference(
        reference_id="chr1", reference_fasta=None, cache_dir=tmp_path / "cache"
    )
    assert resolved is None


def test_symlink_reference_creates_relative_link(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    target = data_dir / "ref.fasta"
    target.write_text(">chr1\nACGT\n")

    out_dir = tmp_path / "results"
    out_dir.mkdir()
    link = out_dir / "reference_initial.fasta"
    symlink_reference(link, target)

    assert link.is_symlink()
    assert not Path(link.readlink()).is_absolute()
    assert link.read_text() == ">chr1\nACGT\n"


def test_symlink_reference_is_noop_if_already_present(tmp_path: Path) -> None:
    target = tmp_path / "ref.fasta"
    target.write_text(">chr1\nACGT\n")
    link = tmp_path / "out" / "reference_initial.fasta"
    symlink_reference(link, target)
    original_link_target = link.readlink()

    target.unlink()  # prove the second call doesn't need the target at all
    symlink_reference(link, target)  # should not raise

    assert link.readlink() == original_link_target


def test_symlink_reference_noop_for_existing_broken_symlink(tmp_path: Path) -> None:
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    link = out_dir / "reference_initial.fasta"
    link.symlink_to("nonexistent.fasta")  # pre-existing, broken

    target = tmp_path / "ref.fasta"
    target.write_text(">chr1\nACGT\n")
    symlink_reference(link, target)  # should not raise, and not touch it

    assert link.readlink() == Path("nonexistent.fasta")
