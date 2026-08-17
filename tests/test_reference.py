from pathlib import Path

import pytest

from iterated_consensus.reference import (
    ReferenceError,
    load_reference,
    parse_fasta,
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
