import json
from pathlib import Path

import pytest

from iterated_consensus.report import (
    format_elapsed,
    format_identity,
    render_report_html,
    write_report,
)


def test_format_elapsed_seconds() -> None:
    assert format_elapsed(4.2) == "4.2s"
    assert format_elapsed(59.9) == "59.9s"


def test_format_elapsed_minutes() -> None:
    assert format_elapsed(125) == "2m 5s"


def test_format_elapsed_hours() -> None:
    assert format_elapsed(3725) == "1h 2m"


def test_format_identity_none() -> None:
    assert format_identity(None) == "—"


def test_format_identity_value() -> None:
    assert format_identity(98.7134) == "98.71%"


SAMPLE_SUMMARY = {
    "iterations_run": 3,
    "converged": True,
    "total_elapsed_seconds": 12.4,
    "iterations": [
        {
            "iteration": 0,
            "reads_mapped": 100,
            "consensus_length": 1500,
            "identity_to_previous": None,
            "elapsed_seconds": 4.1,
        },
        {
            "iteration": 1,
            "reads_mapped": 142,
            "consensus_length": 1502,
            "identity_to_previous": 95.5,
            "elapsed_seconds": 4.2,
        },
        {
            "iteration": 2,
            "reads_mapped": 150,
            "consensus_length": 1502,
            "identity_to_previous": 100.0,
            "elapsed_seconds": 4.1,
        },
    ],
}


def test_render_report_html_contains_key_facts() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "<!doctype html>" in html.lower()
    assert "Converged" in html
    assert "1,502" in html  # final consensus length, thousands-separated
    assert "100.00%" in html  # final identity
    assert "12.4s" in html  # total time
    assert "<table" in html
    assert "95.50%" in html  # per-row identity


def test_render_report_html_stopped_status() -> None:
    summary = {**SAMPLE_SUMMARY, "converged": False}
    html = render_report_html(summary)
    assert "Stopped (max iterations reached)" in html
    assert 'class="badge stopped"' in html


def test_render_report_html_includes_chart_with_enough_points() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "<svg" in html
    assert "<path" in html


def test_render_report_html_omits_chart_with_too_few_points() -> None:
    summary = {
        "iterations_run": 1,
        "converged": False,
        "total_elapsed_seconds": 4.1,
        "iterations": [SAMPLE_SUMMARY["iterations"][0]],
    }
    html = render_report_html(summary)
    assert "<svg" not in html


def test_render_report_html_empty_iterations() -> None:
    summary = {"iterations_run": 0, "converged": False, "total_elapsed_seconds": 0.0, "iterations": []}
    html = render_report_html(summary)
    assert "<!doctype html>" in html.lower()
    assert "—" in html  # final identity placeholder


def test_write_report_reads_summary_and_writes_index_html(tmp_path: Path) -> None:
    (tmp_path / "summary.json").write_text(json.dumps(SAMPLE_SUMMARY))
    report_path = write_report(tmp_path)
    assert report_path == tmp_path / "index.html"
    assert report_path.exists()
    assert "Converged" in report_path.read_text()


def test_write_report_missing_summary_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        write_report(tmp_path)
