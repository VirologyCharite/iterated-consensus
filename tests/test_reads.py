import gzip
from pathlib import Path

import pytest

from iterated_consensus.reads import ReadsCatCache
from iterated_consensus.templating import ReadsList


def test_single_file_returned_as_is(tmp_path: Path) -> None:
    f = tmp_path / "a.fastq"
    f.write_text("@r1\nACGT\n+\nIIII\n")
    cache = ReadsCatCache(tmp_path / "work")
    result = cache.resolve(ReadsList("reads_1", (f,)))
    assert result == f


def test_concatenates_multiple_plain_files(tmp_path: Path) -> None:
    a = tmp_path / "a.fastq"
    b = tmp_path / "b.fastq"
    a.write_text("@r1\nAAAA\n+\nIIII\n")
    b.write_text("@r2\nCCCC\n+\nIIII\n")
    cache = ReadsCatCache(tmp_path / "work")
    result = cache.resolve(ReadsList("reads_1", (a, b)))
    assert result.read_text() == a.read_text() + b.read_text()
    assert result.suffix == ".fastq"


def test_concatenates_multiple_gz_files(tmp_path: Path) -> None:
    a = tmp_path / "a.fastq.gz"
    b = tmp_path / "b.fastq.gz"
    with gzip.open(a, "wt") as f:
        f.write("@r1\nAAAA\n+\nIIII\n")
    with gzip.open(b, "wt") as f:
        f.write("@r2\nCCCC\n+\nIIII\n")
    cache = ReadsCatCache(tmp_path / "work")
    result = cache.resolve(ReadsList("reads_1", (a, b)))
    assert result.name.endswith(".fastq.gz")
    with gzip.open(result, "rt") as f:
        content = f.read()
    assert content == "@r1\nAAAA\n+\nIIII\n@r2\nCCCC\n+\nIIII\n"


def test_mixed_gz_and_plain_raises(tmp_path: Path) -> None:
    a = tmp_path / "a.fastq.gz"
    b = tmp_path / "b.fastq"
    a.write_bytes(b"x")
    b.write_text("x")
    cache = ReadsCatCache(tmp_path / "work")
    with pytest.raises(ValueError, match="mixes gzipped"):
        cache.resolve(ReadsList("reads_1", (a, b)))


def test_result_is_cached_across_calls(tmp_path: Path) -> None:
    a = tmp_path / "a.fastq"
    b = tmp_path / "b.fastq"
    a.write_text("A")
    b.write_text("B")
    cache = ReadsCatCache(tmp_path / "work")
    reads_list = ReadsList("reads_1", (a, b))
    first = cache.resolve(reads_list)
    a.write_text("CHANGED")  # mutate source after first concat
    second = cache.resolve(reads_list)
    assert first == second
    assert second.read_text() == "AB"  # not re-concatenated


def test_different_names_are_independent(tmp_path: Path) -> None:
    a = tmp_path / "a.fastq"
    b = tmp_path / "b.fastq"
    a.write_text("A")
    b.write_text("B")
    cache = ReadsCatCache(tmp_path / "work")
    r1 = cache.resolve(ReadsList("reads_1", (a, b)))
    r2 = cache.resolve(ReadsList("reads_2", (a, b)))
    assert r1 != r2
