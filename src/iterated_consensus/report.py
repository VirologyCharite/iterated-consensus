"""Renders summary.json into a self-contained, human-readable index.html."""

from __future__ import annotations

import html
import json
from pathlib import Path

from .metrics import ambiguous_count as _composition_ambiguous_count


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m {secs:.0f}s"
    hours, minutes = divmod(minutes, 60)
    return f"{int(hours)}h {int(minutes)}m"


def format_identity(identity: float | None) -> str:
    return "—" if identity is None else f"{identity:.2f}%"


def _linear_ticks(min_v: float, max_v: float, count: int) -> list[float]:
    """`count` evenly spaced values from `min_v` to `max_v`, inclusive."""
    if min_v == max_v:
        return [min_v]
    return [min_v + i * (max_v - min_v) / (count - 1) for i in range(count)]


def _integer_ticks(min_v: int, max_v: int, max_count: int) -> list[int]:
    """Integer tick values spanning [min_v, max_v] -- every integer if that's
    few enough to fit, else an evenly spaced, deduplicated subset."""
    if max_v - min_v + 1 <= max_count:
        return list(range(min_v, max_v + 1))
    return sorted({round(v) for v in _linear_ticks(min_v, max_v, max_count)})


def _tick_decimal_places(step: float) -> int:
    """Decimal places so ticks `step` apart render as visually distinct
    values -- 2 in the normal case, more only when a run's identity values
    are converging in a very tight band (e.g. 99.988%-100.000%), where 2
    decimals would round two different ticks to the same displayed text."""
    places = 2
    while places < 6 and step < 10 ** (-places):
        places += 1
    return places


def _chart_svg(iterations: list[dict]) -> str:
    points = [
        (it["iteration"], it["identity_to_previous"])
        for it in iterations
        if it["identity_to_previous"] is not None
    ]
    if len(points) < 2:
        return ""

    width, height = 640, 220
    # pad_left needs room for both the y tick labels (up to 10 chars in the
    # rare case of a very tightly converging run, e.g. "100.000000") and the
    # rotated axis title to their left, or the two overlap.
    pad_left, pad_right, pad_top, pad_bottom = 96, 16, 16, 44
    plot_w = width - pad_left - pad_right
    plot_h = height - pad_top - pad_bottom
    plot_bottom = pad_top + plot_h

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_min, x_max = min(xs), max(xs)
    data_y_min, data_y_max = min(ys), max(ys)
    if data_y_max == data_y_min:
        # A flat line still needs a non-zero domain to plot against, but the
        # tick shown should be the actual (single) value, not a padded one.
        y_min, y_max = data_y_min - 1, data_y_max + 1
        y_ticks = [data_y_min]
        y_places = 2
    else:
        y_min, y_max = data_y_min, data_y_max
        y_ticks = _linear_ticks(y_min, y_max, 4)
        y_places = _tick_decimal_places((y_max - y_min) / 3)

    def sx(x: float) -> float:
        if x_max == x_min:
            return pad_left
        return pad_left + (x - x_min) / (x_max - x_min) * plot_w

    def sy(y: float) -> float:
        return pad_top + plot_h - (y - y_min) / (y_max - y_min) * plot_h

    path_d = " ".join(
        f"{'M' if i == 0 else 'L'}{sx(x):.1f},{sy(y):.1f}" for i, (x, y) in enumerate(points)
    )
    circles = "".join(
        f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.5" class="chart-point">'
        f"<title>iter {x}: {y:.2f}%</title></circle>"
        for x, y in points
    )

    axes = (
        f'<line x1="{pad_left:.1f}" y1="{pad_top:.1f}" x2="{pad_left:.1f}" y2="{plot_bottom:.1f}" '
        f'class="chart-axis"/>'
        f'<line x1="{pad_left:.1f}" y1="{plot_bottom:.1f}" x2="{pad_left + plot_w:.1f}" '
        f'y2="{plot_bottom:.1f}" class="chart-axis"/>'
    )

    x_ticks = _integer_ticks(int(x_min), int(x_max), max_count=8)
    x_tick_marks = "".join(
        f'<line x1="{sx(t):.1f}" y1="{plot_bottom:.1f}" x2="{sx(t):.1f}" y2="{plot_bottom + 5:.1f}" '
        f'class="chart-tick"/>'
        f'<text x="{sx(t):.1f}" y="{plot_bottom + 18:.1f}" class="chart-tick-label" '
        f'text-anchor="middle">{t}</text>'
        for t in x_ticks
    )
    y_tick_marks = "".join(
        f'<line x1="{pad_left - 5:.1f}" y1="{sy(t):.1f}" x2="{pad_left:.1f}" y2="{sy(t):.1f}" '
        f'class="chart-tick"/>'
        f'<text x="{pad_left - 9:.1f}" y="{sy(t):.1f}" class="chart-tick-label" '
        f'text-anchor="end" dominant-baseline="middle">{t:.{y_places}f}</text>'
        for t in y_ticks
    )

    x_title = (
        f'<text x="{pad_left + plot_w / 2:.1f}" y="{height - 6:.1f}" class="chart-axis-title" '
        f'text-anchor="middle">Iteration</text>'
    )
    y_title_cy = pad_top + plot_h / 2
    y_title_x = 12
    y_title = (
        f'<text x="{y_title_x}" y="{y_title_cy:.1f}" class="chart-axis-title" text-anchor="middle" '
        f'transform="rotate(-90 {y_title_x} {y_title_cy:.1f})">Identity to previous (%)</text>'
    )

    return (
        f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" '
        f'aria-label="Identity to previous consensus, by iteration">'
        f"{axes}{x_tick_marks}{y_tick_marks}{x_title}{y_title}"
        f'<path d="{path_d}" class="chart-line" fill="none"/>{circles}</svg>'
    )


def format_ambiguous_count(count: int | None) -> str:
    return "—" if count is None else f"{count:,}"


def _ambiguous_count(composition: dict[str, int]) -> int | None:
    """Count of consensus characters that aren't plain A/C/G/T -- IUPAC
    ambiguity codes, N, gaps, anything else. None if there's no composition
    data at all (a resumed run's summary.json may predate this field)."""
    return _composition_ambiguous_count(composition) if composition else None


def _table_rows(iterations: list[dict]) -> str:
    rows = []
    for it in iterations:
        # .get, not [] -- a resumed run's summary.json may predate this field.
        md5 = it.get("consensus_md5") or "—"
        ambiguous = format_ambiguous_count(_ambiguous_count(it.get("composition", {})))
        rows.append(
            "<tr>"
            f"<td>{it['iteration']}</td>"
            f"<td>{it['reads_mapped']:,}</td>"
            f"<td>{it['consensus_length']:,}</td>"
            f"<td>{ambiguous}</td>"
            f"<td>{format_identity(it['identity_to_previous'])}</td>"
            f"<td>{format_elapsed(it['elapsed_seconds'])}</td>"
            f"<td><code>{md5}</code></td>"
            "</tr>"
        )
    return "\n".join(rows)


_UNAMBIGUOUS_BASES = ("A", "C", "G", "T")
# Canonical display order: the four unambiguous bases, then IUPAC ambiguity
# codes in their usual reference-table order, then anything else (e.g. a gap
# character) alphabetically, last.
_BASE_ORDER = [*_UNAMBIGUOUS_BASES, "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]


def _sorted_composition_items(composition: dict[str, int]) -> list[tuple[str, int]]:
    def sort_key(item: tuple[str, int]) -> tuple[int, str]:
        base = item[0]
        return (_BASE_ORDER.index(base), "") if base in _BASE_ORDER else (len(_BASE_ORDER), base)

    return sorted(composition.items(), key=sort_key)


def _composition_section(composition: dict[str, int]) -> str:
    """Nucleotide breakdown of the final consensus -- every character seen
    (including ambiguity codes/gaps), what fraction is unambiguous (plain
    A/C/G/T), and GC content computed over just those unambiguous bases
    (excluding N and other ambiguity codes, which aren't G or C by
    definition and would just dilute the percentage if counted)."""
    total = sum(composition.values())
    if not total:
        return ""

    unambiguous = sum(composition.get(b, 0) for b in _UNAMBIGUOUS_BASES)
    gc = composition.get("G", 0) + composition.get("C", 0)
    pct_unambiguous = 100.0 * unambiguous / total
    pct_gc = 100.0 * gc / unambiguous if unambiguous else 0.0

    rows = "".join(
        f"<tr><td>{html.escape(base)}</td><td>{count:,}</td>"
        f"<td>{100.0 * count / total:.2f}%</td></tr>"
        for base, count in _sorted_composition_items(composition)
    )
    return f"""<section class="card">
  <h2>Consensus composition</h2>
  <div class="mini-stats">
    <div class="mini-stat"><span class="mini-stat-value">{pct_unambiguous:.2f}%</span><span class="mini-stat-label">Unambiguous (A/C/G/T)</span></div>
    <div class="mini-stat"><span class="mini-stat-value">{pct_gc:.2f}%</span><span class="mini-stat-label">GC content</span></div>
  </div>
  <table>
    <thead><tr><th>Base</th><th>Count</th><th>% of consensus</th></tr></thead>
    <tbody>
{rows}
    </tbody>
  </table>
</section>"""


def _collect_tool_versions(iterations: list[dict]) -> dict[str, list[tuple[int, str]]]:
    """tool name -> [(iteration, version), ...], for every iteration that
    reported that tool -- a tool absent from earlier iterations (e.g. a
    BAM-start run's mapper-less iter_000) just doesn't contribute an entry."""
    by_tool: dict[str, list[tuple[int, str]]] = {}
    for it in iterations:
        for name, version in it.get("tool_versions", {}).items():
            by_tool.setdefault(name, []).append((it["iteration"], version))
    return by_tool


def _tool_versions_section(iterations: list[dict]) -> str:
    """Each tool's version, shown once -- unless it changed partway through
    the run, in which case every distinct value seen is flagged prominently
    instead of silently picking one."""
    by_tool = _collect_tool_versions(iterations)
    if not by_tool:
        return ""

    rows = []
    for name in sorted(by_tool):
        occurrences = by_tool[name]
        distinct = {version for _, version in occurrences}
        safe_name = html.escape(name)
        if len(distinct) == 1:
            version = html.escape(occurrences[0][1])
            rows.append(
                f'<div class="tool-version">'
                f'<span class="tool-name">{safe_name}</span>'
                f'<pre class="tool-version-value">{version}</pre>'
                "</div>"
            )
        else:
            entries = "".join(
                f"<li>iter {iteration}: <pre>{html.escape(version)}</pre></li>"
                for iteration, version in occurrences
            )
            rows.append(
                f'<div class="tool-version tool-version-changed">'
                f'<div class="version-warning">Version changed during this run</div>'
                f'<span class="tool-name">{safe_name}</span>'
                f"<ul>{entries}</ul>"
                "</div>"
            )
    return f'<div class="tool-versions"><h3>Tool versions</h3>{"".join(rows)}</div>'


def _command_log_html(cmd: dict) -> str:
    log = cmd.get("log", "")
    # The log's first line is always "$ <command>" (see commands.run_command)
    # -- always shown, so the command itself stays visible even when it
    # produced nothing; only the rest counts as actual output.
    header, _, output = log.partition("\n")
    body = (
        f"<pre>{html.escape(log)}</pre>"
        if output.strip()
        else f'<pre>{html.escape(header)}</pre><p class="no-output">(no output)</p>'
    )
    return (
        '<details class="cmd-log"><summary>'
        f'{html.escape(cmd["name"])} — {format_elapsed(cmd["elapsed_seconds"])}'
        f"</summary>{body}</details>"
    )


def _command_group_html(title: str, commands: list[dict]) -> str:
    if not commands:
        return ""
    count = len(commands)
    bodies = "".join(_command_log_html(c) for c in commands)
    return (
        f'<details class="iter-log"><summary>{html.escape(title)} '
        f'({count} command{"s" if count != 1 else ""})</summary>{bodies}</details>'
    )


def _logs_section(iterations: list[dict], output_commands: list[dict]) -> str:
    body = _tool_versions_section(iterations)
    body += "".join(
        _command_group_html(f"Iteration {it['iteration']}", it.get("commands", []))
        for it in iterations
    )
    body += _command_group_html("Final output", output_commands)
    if not body:
        return ""
    return f'<section class="card"><h2>Logs</h2>{body}</section>'


_CSS = """
:root {
  --bg: #f7f8fa; --card-bg: #ffffff; --text: #1a1d23; --muted: #6b7280;
  --border: #e5e7eb; --accent: #2563eb;
  --converged: #16a34a; --converged-bg: rgba(22,163,74,0.14);
  --stopped: #b45309; --stopped-bg: rgba(217,119,6,0.16);
  --danger: #dc2626; --danger-bg: rgba(220,38,38,0.12);
  --cycle: #7c3aed; --cycle-bg: rgba(124,58,237,0.12);
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #14161a; --card-bg: #1e2127; --text: #e5e7eb; --muted: #9ca3af;
    --border: #2d323b; --accent: #60a5fa;
    --converged: #4ade80; --converged-bg: rgba(74,222,128,0.16);
    --stopped: #fbbf24; --stopped-bg: rgba(251,191,36,0.16);
    --danger: #f87171; --danger-bg: rgba(248,113,113,0.18);
    --cycle: #a78bfa; --cycle-bg: rgba(167,139,250,0.18);
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 2rem 1.5rem 4rem; background: var(--bg); color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  line-height: 1.5;
}
header { max-width: 960px; margin: 0 auto 2rem; display: flex; align-items: center; gap: 1rem; flex-wrap: wrap; }
h1 { font-size: 1.5rem; margin: 0; }
.badge { padding: 0.25rem 0.75rem; border-radius: 999px; font-size: 0.85rem; font-weight: 600; }
.badge.converged { background: var(--converged-bg); color: var(--converged); }
.badge.stopped { background: var(--stopped-bg); color: var(--stopped); }
.badge.cycle { background: var(--cycle-bg); color: var(--cycle); }
.cycle-notice { border: 2px solid var(--cycle); background: var(--cycle-bg); }
.cycle-notice h2 { color: var(--cycle); }
.cycle-notice p { margin: 0; }
.stats {
  max-width: 960px; margin: 0 auto 1.5rem;
  display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 1rem;
}
.stat { background: var(--card-bg); border: 1px solid var(--border); border-radius: 12px; padding: 1rem 1.25rem; }
.stat-value { display: block; font-size: 1.5rem; font-weight: 700; }
.stat-label { display: block; font-size: 0.85rem; color: var(--muted); margin-top: 0.25rem; }
.card {
  max-width: 960px; margin: 0 auto 1.5rem; background: var(--card-bg); border: 1px solid var(--border);
  border-radius: 12px; padding: 1.25rem 1.5rem; overflow-x: auto;
}
.card h2 { font-size: 1.05rem; margin: 0 0 1rem; }
.mini-stats { display: flex; flex-wrap: wrap; gap: 1.5rem; margin: 0 0 1rem; }
.mini-stat-value { display: block; font-size: 1.25rem; font-weight: 700; }
.mini-stat-label { display: block; font-size: 0.8rem; color: var(--muted); margin-top: 0.15rem; }
.chart { width: 100%; height: auto; }
.chart-line { stroke: var(--accent); stroke-width: 2; }
.chart-point { fill: var(--accent); }
.chart-axis { stroke: var(--muted); stroke-width: 1; }
.chart-tick { stroke: var(--border); stroke-width: 1; }
.chart-tick-label { fill: var(--muted); font-size: 10px; }
.chart-axis-title { fill: var(--muted); font-size: 11px; }
table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
th, td { text-align: left; padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); font-size: 0.9rem; }
th { color: var(--muted); font-weight: 600; font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.03em; }
tr:last-child td { border-bottom: none; }
code, pre { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }
.tool-versions { margin-bottom: 1rem; }
.tool-versions h3 { font-size: 0.95rem; margin: 0 0 0.75rem; }
.tool-version { margin-bottom: 0.75rem; }
.tool-name { font-weight: 600; font-size: 0.9rem; }
.tool-version-value {
  margin: 0.35rem 0 0; padding: 0.5rem 0.75rem; background: var(--bg); border: 1px solid var(--border);
  border-radius: 6px; font-size: 0.8rem; white-space: pre-wrap; word-break: break-word;
}
.tool-version-changed {
  border: 2px solid var(--danger); background: var(--danger-bg); border-radius: 8px; padding: 0.75rem 1rem;
}
.tool-version-changed ul { margin: 0.5rem 0 0; padding-left: 1.25rem; }
.tool-version-changed li { margin-bottom: 0.4rem; }
.version-warning {
  color: var(--danger); font-weight: 700; text-transform: uppercase; font-size: 0.8rem; letter-spacing: 0.03em;
  margin-bottom: 0.35rem;
}
.iter-log {
  border: 1px solid var(--border); border-radius: 8px; margin-bottom: 0.6rem; padding: 0 0.9rem 0.1rem;
}
.iter-log > summary {
  cursor: pointer; font-weight: 600; padding: 0.6rem 0; font-size: 0.9rem; list-style: none;
}
.iter-log > summary::-webkit-details-marker { display: none; }
.iter-log > summary::before { content: "▸ "; }
.iter-log[open] > summary::before { content: "▾ "; }
.cmd-log { margin: 0 0 0.6rem 1rem; border-left: 2px solid var(--border); padding-left: 0.85rem; }
.cmd-log > summary {
  cursor: pointer; font-size: 0.85rem; color: var(--text); list-style: none; padding: 0.2rem 0;
}
.cmd-log > summary::-webkit-details-marker { display: none; }
.cmd-log pre {
  margin: 0.4rem 0 0.6rem; padding: 0.6rem 0.75rem; background: var(--bg); border: 1px solid var(--border);
  border-radius: 6px; font-size: 0.78rem; overflow-x: auto; white-space: pre-wrap; word-break: break-word;
}
.no-output { color: var(--muted); font-style: italic; font-size: 0.8rem; margin: 0.3rem 0 0.6rem; }
footer { max-width: 960px; margin: 2rem auto 0; color: var(--muted); font-size: 0.8rem; text-align: center; }
"""


def _cycle_notice(cycle: dict) -> str:
    first = cycle["first_iteration"]
    repeat = cycle["repeat_iteration"]
    period = repeat - first
    md5 = html.escape(cycle["consensus_md5"])
    message = (
        f"Iteration {repeat}'s consensus is identical to iteration {first}'s "
        f"(MD5 <code>{md5}</code>) -- the run is oscillating among a set of {period} "
        "sequence(s) rather than settling on one, so it was stopped rather than run to "
        f"max_iterations chasing a fixed point that doesn't exist. Iteration {first}'s "
        "consensus was used as the final result below; see the per-iteration MD5 column "
        "for the full repeating pattern."
    )
    return f'<section class="card cycle-notice"><h2>Cycle detected</h2><p>{message}</p></section>'


def render_report_html(summary: dict) -> str:
    iterations: list[dict] = summary["iterations"]
    converged: bool = summary["converged"]
    total_elapsed: float = summary["total_elapsed_seconds"]
    cycle: dict | None = summary.get("cycle")
    # .get, not [] -- a resumed run's summary.json may predate final_iteration.
    final_iteration_num = summary.get("final_iteration", iterations[-1]["iteration"] if iterations else None)
    final = next((it for it in iterations if it["iteration"] == final_iteration_num), None)

    if cycle is not None:
        period = cycle["repeat_iteration"] - cycle["first_iteration"]
        status_label = f"Cycle detected (period {period})"
        status_class = "cycle"
    elif converged:
        status_label = "Converged"
        status_class = "converged"
    else:
        status_label = "Stopped (max iterations reached)"
        status_class = "stopped"
    final_length = final["consensus_length"] if final else 0
    final_reads_mapped = final["reads_mapped"] if final else 0

    cycle_notice = _cycle_notice(cycle) if cycle is not None else ""
    chart_svg = _chart_svg(iterations)
    chart_section = (
        f"<section class=\"card\"><h2>Identity to previous consensus</h2>{chart_svg}</section>"
        if chart_svg
        else ""
    )
    composition_section = _composition_section(final.get("composition", {}) if final else {})
    logs_section = _logs_section(iterations, summary.get("output_commands", []))

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>iterated-consensus run report</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>{_CSS}</style>
</head>
<body>
<header>
  <h1>iterated-consensus run report</h1>
  <span class="badge {status_class}">{status_label}</span>
</header>

{cycle_notice}

<section class="stats">
  <div class="stat"><span class="stat-value">{len(iterations)}</span><span class="stat-label">Iterations run</span></div>
  <div class="stat"><span class="stat-value">{final_length:,}</span><span class="stat-label">Final consensus length</span></div>
  <div class="stat"><span class="stat-value">{final_reads_mapped:,}</span><span class="stat-label">Final mapped reads</span></div>
  <div class="stat"><span class="stat-value">{format_elapsed(total_elapsed)}</span><span class="stat-label">Total time</span></div>
</section>

{chart_section}

{composition_section}

<section class="card">
  <h2>Per-iteration detail</h2>
  <table>
    <thead><tr><th>Iteration</th><th>Reads mapped</th><th>Consensus length</th><th>Ambiguous</th><th>Identity to previous</th><th>Elapsed</th><th>Consensus MD5</th></tr></thead>
    <tbody>
{_table_rows(iterations)}
    </tbody>
  </table>
</section>

{logs_section}

<footer>Generated by iterated-consensus.</footer>
</body>
</html>
"""


def write_report(out_dir: Path) -> Path:
    """Read `out_dir/summary.json` (plus each iteration's own stats.json, for
    tool-versions/command-log/composition detail that isn't duplicated into
    summary.json) and write `out_dir/index.html`."""
    summary = json.loads((out_dir / "summary.json").read_text())
    for it in summary["iterations"]:
        stats_path = out_dir / f"iter_{it['iteration']:03d}" / "stats.json"
        if stats_path.exists():
            stats = json.loads(stats_path.read_text())
            it["tool_versions"] = stats.get("tool_versions", {})
            it["commands"] = stats.get("commands", [])
            it["composition"] = stats.get("composition", {})
    report_path = out_dir / "index.html"
    report_path.write_text(render_report_html(summary))
    return report_path
