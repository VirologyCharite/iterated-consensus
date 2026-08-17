from pathlib import Path

import pytest

from iterated_consensus.templating import ReadsList, TemplateError, render


@pytest.fixture
def reads_1() -> ReadsList:
    return ReadsList("reads_1", (Path("a.fq"), Path("b.fq"), Path("c.fq")))


@pytest.fixture
def values(reads_1: ReadsList) -> dict[str, object]:
    return {
        "reference": "ref.fa",
        "threads": 4,
        "reads_1": reads_1,
    }


def test_bare_scalar(values: dict[str, object]) -> None:
    assert render("{reference}", values) == "ref.fa"


def test_scalar_non_string_value(values: dict[str, object]) -> None:
    assert render("-p {threads}", values) == "-p 4"


def test_bare_list_is_space_joined(values: dict[str, object]) -> None:
    assert render("{reads_1}", values) == "a.fq b.fq c.fq"


def test_list_with_separator(values: dict[str, object]) -> None:
    assert render("{reads_1:,}", values) == "a.fq,b.fq,c.fq"


def test_list_with_prefix_no_space(values: dict[str, object]) -> None:
    assert render("{-1:reads_1}", values) == "-1a.fq -1b.fq -1c.fq"


def test_list_with_prefix_and_trailing_space(values: dict[str, object]) -> None:
    assert render("{-1 :reads_1}", values) == "-1 a.fq -1 b.fq -1 c.fq"


def test_list_with_prefix_and_separator(values: dict[str, object]) -> None:
    assert render("{-1:reads_1:,}", values) == "-1a.fq,-1b.fq,-1c.fq"


def test_full_command_line(values: dict[str, object]) -> None:
    rendered = render("bowtie2 -x idx {reads_1:,} -p {threads}", values)
    assert rendered == "bowtie2 -x idx a.fq,b.fq,c.fq -p 4"


def test_cat_resolves_via_callback(values: dict[str, object]) -> None:
    calls = []

    def resolver(reads_list: ReadsList) -> Path:
        calls.append(reads_list)
        return Path("/tmp/merged.fq")

    result = render("bwa mem ref {cat:reads_1}", values, cat_resolver=resolver)
    assert result == "bwa mem ref /tmp/merged.fq"
    assert calls == [values["reads_1"]]


def test_cat_without_resolver_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match="cat_resolver"):
        render("{cat:reads_1}", values)


def test_cat_with_separator_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match="cannot be combined with a separator"):
        render("{cat:reads_1:,}", values, cat_resolver=lambda r: Path("x"))


def test_unknown_placeholder_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match="unknown placeholder"):
        render("{nope}", values)


def test_modifier_on_scalar_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match="not a file list"):
        render("{-p:reference}", values)


def test_too_many_colons_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match="too many ':'"):
        render("{a:b:c:d}", values)


def test_empty_reads_list_raises() -> None:
    with pytest.raises(TemplateError, match="no files"):
        ReadsList("reads_1", ())


def test_literal_braces_are_escaped(values: dict[str, object]) -> None:
    assert render("{{not a placeholder}}", values) == "{not a placeholder}"


def test_unmatched_open_brace_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match=r"unmatched '\{'"):
        render("{reference", values)


def test_unmatched_close_brace_raises(values: dict[str, object]) -> None:
    with pytest.raises(TemplateError, match=r"unmatched '\}'"):
        render("reference}", values)


def test_second_list_disambiguates_independently() -> None:
    reads_2 = ReadsList("reads_2", (Path("x.fq"), Path("y.fq")))
    values = {"reads_2": reads_2}
    # 'reads_2' is a known list name, so this is name:sep, not prefix:name.
    assert render("{reads_2:+}", values) == "x.fq+y.fq"
