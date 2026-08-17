"""Renders summary.json into a self-contained, human-readable index.html."""

from __future__ import annotations

import json
from pathlib import Path


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


def _chart_svg(iterations: list[dict]) -> str:
    points = [
        (it["iteration"], it["identity_to_previous"])
        for it in iterations
        if it["identity_to_previous"] is not None
    ]
    if len(points) < 2:
        return ""

    width, height, pad = 640, 160, 24
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if y_max == y_min:
        y_min, y_max = y_min - 1, y_max + 1

    def sx(x: float) -> float:
        if x_max == x_min:
            return pad
        return pad + (x - x_min) / (x_max - x_min) * (width - 2 * pad)

    def sy(y: float) -> float:
        return height - pad - (y - y_min) / (y_max - y_min) * (height - 2 * pad)

    path_d = " ".join(
        f"{'M' if i == 0 else 'L'}{sx(x):.1f},{sy(y):.1f}" for i, (x, y) in enumerate(points)
    )
    circles = "".join(
        f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.5" class="chart-point">'
        f"<title>iter {x}: {y:.2f}%</title></circle>"
        for x, y in points
    )
    return (
        f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" '
        f'aria-label="Identity to previous consensus, by iteration">'
        f'<path d="{path_d}" class="chart-line" fill="none"/>{circles}</svg>'
    )


def _table_rows(iterations: list[dict]) -> str:
    rows = []
    for it in iterations:
        rows.append(
            "<tr>"
            f"<td>{it['iteration']}</td>"
            f"<td>{it['reads_mapped']:,}</td>"
            f"<td>{it['consensus_length']:,}</td>"
            f"<td>{format_identity(it['identity_to_previous'])}</td>"
            f"<td>{format_elapsed(it['elapsed_seconds'])}</td>"
            "</tr>"
        )
    return "\n".join(rows)


_CSS = """
:root {
  --bg: #f7f8fa; --card-bg: #ffffff; --text: #1a1d23; --muted: #6b7280;
  --border: #e5e7eb; --accent: #2563eb;
  --converged: #16a34a; --converged-bg: rgba(22,163,74,0.14);
  --stopped: #b45309; --stopped-bg: rgba(217,119,6,0.16);
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #14161a; --card-bg: #1e2127; --text: #e5e7eb; --muted: #9ca3af;
    --border: #2d323b; --accent: #60a5fa;
    --converged: #4ade80; --converged-bg: rgba(74,222,128,0.16);
    --stopped: #fbbf24; --stopped-bg: rgba(251,191,36,0.16);
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
.chart { width: 100%; height: auto; }
.chart-line { stroke: var(--accent); stroke-width: 2; }
.chart-point { fill: var(--accent); }
table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
th, td { text-align: left; padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); font-size: 0.9rem; }
th { color: var(--muted); font-weight: 600; font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.03em; }
tr:last-child td { border-bottom: none; }
footer { max-width: 960px; margin: 2rem auto 0; color: var(--muted); font-size: 0.8rem; text-align: center; }
"""


def render_report_html(summary: dict) -> str:
    iterations: list[dict] = summary["iterations"]
    converged: bool = summary["converged"]
    total_elapsed: float = summary["total_elapsed_seconds"]
    final = iterations[-1] if iterations else None

    status_label = "Converged" if converged else "Stopped (max iterations reached)"
    status_class = "converged" if converged else "stopped"
    final_length = final["consensus_length"] if final else 0
    final_identity = format_identity(final["identity_to_previous"]) if final else "—"

    chart_svg = _chart_svg(iterations)
    chart_section = (
        f"<section class=\"card\"><h2>Identity to previous consensus</h2>{chart_svg}</section>"
        if chart_svg
        else ""
    )

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

<section class="stats">
  <div class="stat"><span class="stat-value">{len(iterations)}</span><span class="stat-label">Iterations run</span></div>
  <div class="stat"><span class="stat-value">{final_length:,}</span><span class="stat-label">Final consensus length</span></div>
  <div class="stat"><span class="stat-value">{final_identity}</span><span class="stat-label">Final identity to previous</span></div>
  <div class="stat"><span class="stat-value">{format_elapsed(total_elapsed)}</span><span class="stat-label">Total time</span></div>
</section>

{chart_section}

<section class="card">
  <h2>Per-iteration detail</h2>
  <table>
    <thead><tr><th>Iteration</th><th>Reads mapped</th><th>Consensus length</th><th>Identity to previous</th><th>Elapsed</th></tr></thead>
    <tbody>
{_table_rows(iterations)}
    </tbody>
  </table>
</section>

<footer>Generated by iterated-consensus.</footer>
</body>
</html>
"""


def write_report(out_dir: Path) -> Path:
    """Read `out_dir/summary.json` and write `out_dir/index.html` from it."""
    summary = json.loads((out_dir / "summary.json").read_text())
    report_path = out_dir / "index.html"
    report_path.write_text(render_report_html(summary))
    return report_path
