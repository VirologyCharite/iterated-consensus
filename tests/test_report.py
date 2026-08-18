import json
import re
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
    assert "Final mapped reads" in html
    assert "150" in html  # final iteration's reads_mapped
    assert "12.4s" in html  # total time
    assert "<table" in html
    assert "95.50%" in html  # per-row identity
    assert "100.00%" in html  # per-row identity (final iteration)


def test_render_report_html_table_includes_consensus_md5() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**it, "consensus_md5": f"md5-for-iter-{it['iteration']}"}
            for it in SAMPLE_SUMMARY["iterations"]
        ],
    }
    html = render_report_html(summary)
    assert "Consensus MD5" in html
    assert "<code>md5-for-iter-0</code>" in html
    assert "<code>md5-for-iter-2</code>" in html


def test_render_report_html_table_tolerates_missing_consensus_md5() -> None:
    # A resumed run's summary.json may predate this field.
    html = render_report_html(SAMPLE_SUMMARY)
    assert "<code>—</code>" in html


def test_render_report_html_table_includes_ambiguous_count() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**SAMPLE_SUMMARY["iterations"][0]},
            {**SAMPLE_SUMMARY["iterations"][1]},
            {
                **SAMPLE_SUMMARY["iterations"][2],
                "composition": {"A": 30, "C": 20, "G": 20, "T": 25, "N": 4, "R": 1},
            },
        ],
    }
    html = render_report_html(summary)
    assert "Ambiguous" in html
    assert "<td>5</td>" in html  # N(4) + R(1)


def test_render_report_html_table_ambiguous_count_dash_when_no_composition() -> None:
    # A resumed run's summary.json may predate per-iteration composition.
    html = render_report_html(SAMPLE_SUMMARY)
    assert "<td>—</td>" in html


def test_render_report_html_stopped_status() -> None:
    summary = {**SAMPLE_SUMMARY, "converged": False}
    html = render_report_html(summary)
    assert "Stopped (max iterations reached)" in html
    assert 'class="badge stopped"' in html


def test_render_report_html_no_cycle_notice_when_no_cycle() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "Cycle detected" not in html
    assert '<section class="card cycle-notice">' not in html


def test_render_report_html_cycle_badge_and_notice() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "converged": False,
        "cycle": {"first_iteration": 0, "repeat_iteration": 2, "consensus_md5": "abc123"},
        "final_iteration": 0,
    }
    html = render_report_html(summary)
    assert 'class="badge cycle"' in html
    assert "Cycle detected (period 2)" in html
    assert 'class="card cycle-notice"' in html
    assert "iteration 2" in html.lower() or "Iteration 2" in html
    assert "iteration 0" in html.lower() or "Iteration 0" in html
    assert "abc123" in html


def test_render_report_html_cycle_uses_first_occurrence_for_final_stats() -> None:
    """The "Final ..." stat tiles reflect the *first* occurrence of the
    repeated consensus, not the last iteration actually run."""
    summary = {
        "iterations_run": 3,
        "converged": False,
        "total_elapsed_seconds": 9.0,
        "cycle": {"first_iteration": 0, "repeat_iteration": 2, "consensus_md5": "abc123"},
        "final_iteration": 0,
        "iterations": [
            {
                "iteration": 0, "reads_mapped": 111, "consensus_length": 999,
                "identity_to_previous": None, "elapsed_seconds": 3.0,
            },
            {
                "iteration": 1, "reads_mapped": 222, "consensus_length": 888,
                "identity_to_previous": 0.0, "elapsed_seconds": 3.0,
            },
            {
                "iteration": 2, "reads_mapped": 333, "consensus_length": 999,
                "identity_to_previous": 0.0, "elapsed_seconds": 3.0,
            },
        ],
    }
    html = render_report_html(summary)
    # The "Final mapped reads" stat tile uses iteration 0's value (111), not
    # iteration 2's (333) -- even though 333 legitimately still shows up in
    # the per-iteration table below, since that lists every iteration run.
    assert '<span class="stat-value">111</span><span class="stat-label">Final mapped reads</span>' in html
    assert '<span class="stat-value">999</span><span class="stat-label">Final consensus length</span>' in html


def test_render_report_html_includes_chart_with_enough_points() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "<svg" in html
    assert "<path" in html


def test_render_report_html_chart_has_axis_labels_and_ticks() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "Iteration</text>" in html  # x-axis title
    assert "Identity to previous (%)</text>" in html  # y-axis title
    assert 'class="chart-axis"' in html
    assert 'class="chart-tick"' in html
    assert 'class="chart-tick-label"' in html
    # x tick labels for the two non-null-identity iterations
    assert '>1</text>' in html
    assert '>2</text>' in html
    # y tick labels include the data's actual min/max identity, without a
    # redundant "%" (the axis title already says "(%)")
    assert "95.50</text>" in html
    assert "100.00</text>" in html
    assert "95.50%</text>" not in html


def test_render_report_html_chart_ticks_stay_distinct_for_a_tight_identity_range() -> None:
    """A run converging within a fraction of a percent shouldn't render
    duplicate-looking y tick labels just because they're all rounded to a
    fixed 2 decimal places."""
    summary = {
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
                "identity_to_previous": 99.9877,
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
    html = render_report_html(summary)
    tick_labels = re.findall(r'class="chart-tick-label"[^>]*>([^<]+)</text>', html)
    y_tick_labels = [label for label in tick_labels if label not in {"1", "2"}]
    assert len(y_tick_labels) == len(set(y_tick_labels))  # no two identical


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
    assert "Final mapped reads" in html
    assert "Consensus composition" not in html  # nothing to show


def test_write_report_reads_summary_and_writes_index_html(tmp_path: Path) -> None:
    (tmp_path / "summary.json").write_text(json.dumps(SAMPLE_SUMMARY))
    report_path = write_report(tmp_path)
    assert report_path == tmp_path / "index.html"
    assert report_path.exists()
    assert "Converged" in report_path.read_text()


def test_write_report_missing_summary_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        write_report(tmp_path)


def test_render_report_html_no_composition_section_when_nothing_to_show() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "Consensus composition" not in html


def test_render_report_html_composition_section_shows_counts_and_percentages() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**SAMPLE_SUMMARY["iterations"][0]},
            {**SAMPLE_SUMMARY["iterations"][1]},
            {
                **SAMPLE_SUMMARY["iterations"][2],
                "composition": {"A": 30, "C": 20, "G": 20, "T": 25, "N": 5},
            },
        ],
    }
    html = render_report_html(summary)
    assert "Consensus composition" in html
    # unambiguous = 95/100 = 95.00%; GC = (20+20)/95 = 42.11%
    assert "95.00%" in html
    assert "42.11%" in html
    assert "<td>A</td><td>30</td><td>30.00%</td>" in html
    assert "<td>N</td><td>5</td><td>5.00%</td>" in html
    # canonical order: A, C, G, T before N
    a_pos = html.index("<td>A</td>")
    t_pos = html.index("<td>T</td>")
    n_pos = html.index("<td>N</td>")
    assert a_pos < t_pos < n_pos


def test_render_report_html_composition_all_ambiguous_gc_is_zero() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**SAMPLE_SUMMARY["iterations"][0]},
            {**SAMPLE_SUMMARY["iterations"][1]},
            {**SAMPLE_SUMMARY["iterations"][2], "composition": {"N": 10}},
        ],
    }
    html = render_report_html(summary)
    assert "0.00%" in html  # unambiguous and GC both 0%, no division by zero


def test_render_report_html_no_logs_section_when_nothing_to_show() -> None:
    html = render_report_html(SAMPLE_SUMMARY)
    assert "<h2>Logs</h2>" not in html


def test_render_report_html_shows_command_log_with_elapsed_and_output() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {
                **SAMPLE_SUMMARY["iterations"][0],
                "commands": [
                    {"name": "bowtie2_index", "display": "bowtie2-build ref idx", "elapsed_seconds": 1.5, "log": "$ bowtie2-build ref idx\nbuilding index...\n"},
                ],
            },
        ],
    }
    html = render_report_html(summary)
    assert "<h2>Logs</h2>" in html
    assert "Iteration 0" in html
    assert "bowtie2_index" in html
    assert "1.5s" in html
    assert "building index..." in html


def test_render_report_html_command_with_no_output_says_so() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {
                **SAMPLE_SUMMARY["iterations"][0],
                "commands": [
                    {"name": "true_step", "display": "true", "elapsed_seconds": 0.0, "log": "$ true\n"},
                ],
            },
        ],
    }
    html = render_report_html(summary)
    assert "(no output)" in html


def test_render_report_html_escapes_command_log_content() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {
                **SAMPLE_SUMMARY["iterations"][0],
                "commands": [
                    {
                        "name": "evil",
                        "display": "echo <script>",
                        "elapsed_seconds": 0.0,
                        "log": "$ echo <script>\n<script>alert(1)</script>\n",
                    },
                ],
            },
        ],
    }
    html = render_report_html(summary)
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html


def test_render_report_html_tool_version_shown_once_when_stable() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**SAMPLE_SUMMARY["iterations"][0], "tool_versions": {"ivar": "1.4.2"}},
            {**SAMPLE_SUMMARY["iterations"][1], "tool_versions": {"ivar": "1.4.2"}},
        ],
    }
    html = render_report_html(summary)
    assert html.count("1.4.2") == 1
    assert "Version changed" not in html


def test_render_report_html_flags_tool_version_change_prominently() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**SAMPLE_SUMMARY["iterations"][0], "tool_versions": {"ivar": "1.4.1"}},
            {**SAMPLE_SUMMARY["iterations"][1], "tool_versions": {"ivar": "1.4.2"}},
        ],
    }
    html = render_report_html(summary)
    assert "Version changed during this run" in html
    assert "tool-version-changed" in html
    assert "1.4.1" in html
    assert "1.4.2" in html
    assert "iter 0" in html
    assert "iter 1" in html


def test_render_report_html_tool_version_can_span_multiple_lines() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "iterations": [
            {**SAMPLE_SUMMARY["iterations"][0], "tool_versions": {"ivar": "ivar 1.4.2\nbuild: abc123"}},
        ],
    }
    html = render_report_html(summary)
    assert "ivar 1.4.2\nbuild: abc123" in html


def test_render_report_html_shows_final_output_commands() -> None:
    summary = {
        **SAMPLE_SUMMARY,
        "output_commands": [
            {"name": "output_command_00", "display": "samtools faidx final.fasta", "elapsed_seconds": 0.3, "log": "$ samtools faidx final.fasta\n"},
        ],
    }
    html = render_report_html(summary)
    assert "Final output" in html
    assert "samtools faidx final.fasta" in html


def test_write_report_merges_per_iteration_stats_json(tmp_path: Path) -> None:
    (tmp_path / "summary.json").write_text(json.dumps(SAMPLE_SUMMARY))
    iter_dir = tmp_path / "iter_000"
    iter_dir.mkdir()
    (iter_dir / "stats.json").write_text(json.dumps({
        "tool_versions": {"ivar": "1.4.2"},
        "commands": [{"name": "step", "display": "true", "elapsed_seconds": 0.1, "log": "$ true\n"}],
    }))
    report_path = write_report(tmp_path)
    html = report_path.read_text()
    assert "<h2>Logs</h2>" in html
    assert "1.4.2" in html
    assert "Iteration 0" in html


def test_write_report_merges_final_iteration_composition(tmp_path: Path) -> None:
    (tmp_path / "summary.json").write_text(json.dumps(SAMPLE_SUMMARY))
    final_iter_dir = tmp_path / "iter_002"  # SAMPLE_SUMMARY's last iteration
    final_iter_dir.mkdir()
    (final_iter_dir / "stats.json").write_text(json.dumps({"composition": {"A": 10, "T": 10}}))
    report_path = write_report(tmp_path)
    html = report_path.read_text()
    assert "Consensus composition" in html
    assert "100.00%" in html  # all A/T -- fully unambiguous, 0% GC
