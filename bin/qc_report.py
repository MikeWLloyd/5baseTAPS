#!/usr/bin/env python3
"""
TAPS QC batch report generator.

Reads all per-sample DRAGEN-format CSV files (mapping_metrics, methyl_metrics,
wgs_coverage_metrics, duplex_metrics, vc_metrics) from the working directory
and generates a self-contained batch-level HTML QC report.

Usage:
    taps_qc_report.py [--out <output.html>] [csv_file ...]

    If csv_file arguments are provided, only those files are used.
    Otherwise, the current directory is scanned for *.csv files.

Output:
    taps_qc_report.html  (or --out path)
    taps_qc_report.tsv   (flat TSV, one row per sample, all key metrics)
"""

import sys
import csv
import os
import re
import json
import html
import argparse
from pathlib import Path
from collections import defaultdict


# ── CSV suffix → metric type ─────────────────────────────────────────────────

METRIC_SUFFIXES = {
    "mapping_metrics.csv":     "mapping",
    "methyl_metrics.csv":      "methyl",
    "wgs_coverage_metrics.csv":"wgs",
    "duplex_metrics.csv":      "duplex",
    "vc_metrics.csv":          "vc",
}


def detect_metric_type(filename):
    """Detect which metric CSV type this file is from its suffix."""
    name = Path(filename).name
    for suffix, mtype in METRIC_SUFFIXES.items():
        if name.endswith(suffix):
            return mtype
    return None


def detect_sample_id(filename, metric_type):
    """Extract sample ID by stripping known suffixes from filename."""
    name = Path(filename).stem  # remove .csv
    for suffix in METRIC_SUFFIXES:
        base_suffix = suffix[:-4]  # strip .csv → e.g. 'mapping_metrics'
        if name.endswith("." + base_suffix):
            return name[: -(len(base_suffix) + 1)]
        if name.endswith(base_suffix):
            return name[: -len(base_suffix)]
    return name


def parse_csv(path):
    """
    Parse a DRAGEN-format 5-column CSV into nested dict.
    Returns dict: {(section, subsection, metric): (value, percent)}
    """
    data = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if len(row) < 4:
                continue
            section    = row[0].strip()
            subsection = row[1].strip() if len(row) > 1 else ""
            metric     = row[2].strip() if len(row) > 2 else ""
            value      = row[3].strip() if len(row) > 3 else ""
            percent    = row[4].strip() if len(row) > 4 else ""
            data[(section, subsection, metric)] = (value, percent)
    return data


def load_all_csvs(csv_files):
    """
    Load all CSV files and organize by (sample_id, metric_type).
    Returns: {sample_id: {metric_type: {(section,sub,metric): (value,pct)}}}
    """
    all_data = defaultdict(dict)
    for path in csv_files:
        mtype = detect_metric_type(path)
        if mtype is None:
            continue
        sample_id = detect_sample_id(path, mtype)
        all_data[sample_id][mtype] = parse_csv(path)
    return all_data


def get_val(data, section, subsection, metric, default=""):
    """Safely get a value from parsed CSV data dict."""
    return data.get((section, subsection, metric), (default, ""))[0]


def get_pct(data, section, subsection, metric, default=""):
    """Safely get a percent from parsed CSV data dict."""
    return data.get((section, subsection, metric), ("", default))[1]


# ── Summary row extraction ────────────────────────────────────────────────────

def extract_summary(sample_id, sample_data):
    """Extract the key per-sample metrics for the summary table."""
    mapping = sample_data.get("mapping", {})
    methyl  = sample_data.get("methyl",  {})
    wgs     = sample_data.get("wgs",     {})
    duplex  = sample_data.get("duplex",  {})
    vc      = sample_data.get("vc",      {})

    return {
        "Sample":                   sample_id,
        "Total reads":              get_val(mapping, "MAPPING/ALIGNING SUMMARY", "", "Total input reads"),
        "Mapped %":                 get_pct(mapping, "MAPPING/ALIGNING SUMMARY", "", "Mapped reads"),
        "lambda %":                 get_pct(mapping, "SPIKE-IN COMPOSITION", "", "Reads mapped to lambda (negCtrl)"),
        "pUC19 %":                  get_pct(mapping, "SPIKE-IN COMPOSITION", "", "Reads mapped to pUC19 (posCtrl)"),
        "Consensus reads":          get_val(duplex,  "DUPLEX DEDUP SUMMARY", "", "Consensus reads (post-dedup)"),
        "Dedup fold-red.":          get_val(duplex,  "DUPLEX DEDUP SUMMARY", "", "Dedup fold-reduction"),
        "Duplex yield %":           get_val(duplex,  "DUPLEX DEDUP SUMMARY", "", "Duplex yield fraction (%)"),
        "Mean coverage":            get_val(wgs,     "COVERAGE SUMMARY", "", "Average alignment coverage"),
        "PCT >=10x":                get_val(wgs,     "COVERAGE SUMMARY", "", "PCT of genome at coverage >= 10x"),
        "PCT >=20x":                get_val(wgs,     "COVERAGE SUMMARY", "", "PCT of genome at coverage >= 20x"),
        "Global meth %":            get_val(wgs,     "CPG COVERAGE", "", "Global methylation % (CpG)"),
        "CpG covered":              get_val(wgs,     "CPG COVERAGE", "", "Total CpGs called"),
        "lambda meth %":            get_val(methyl,  "METHYL QC", "lambda_negCtrl", "% methylated CpG"),
        "pUC19 meth %":             get_val(methyl,  "METHYL QC", "puc19_posCtrl",  "% methylated CpG"),
        "Total variants":           get_val(vc,      "VARIANT CALLER SUMMARY", "", "Total variants"),
        "SNPs":                     get_val(vc,      "VARIANT CALLER SUMMARY", "", "SNPs"),
        "Ti/Tv":                    get_val(vc,      "VARIANT CALLER SUMMARY", "", "Transition/Transversion ratio"),
    }


# ── TSV export ────────────────────────────────────────────────────────────────

def write_tsv(all_data, tsv_path):
    """Write flat TSV with one row per sample."""
    samples = sorted(all_data.keys())
    rows = [extract_summary(sid, all_data[sid]) for sid in samples]
    if not rows:
        return
    headers = list(rows[0].keys())
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ── HTML generation ───────────────────────────────────────────────────────────

# Minimal Bootstrap 5.3 subset inlined for self-contained offline use.
BOOTSTRAP_CSS = """
/* Bootstrap 5 minimal subset — inlined for offline use */
*,::after,::before{box-sizing:border-box}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;font-size:1rem;line-height:1.5;color:#212529;background-color:#fff}
h1,h2,h3,h4,h5,h6{margin-top:0;margin-bottom:.5rem;font-weight:500;line-height:1.2}
.container{width:100%;padding-right:15px;padding-left:15px;margin-right:auto;margin-left:auto;max-width:1400px}
table{width:100%;caption-side:bottom;border-collapse:collapse}
th{text-align:inherit;text-align:-webkit-match-parent}
thead,tbody,tfoot,tr,td,th{border-color:inherit;border-style:solid;border-width:0}
.table{--bs-table-bg:transparent;--bs-table-striped-bg:#f2f2f2;--bs-table-hover-bg:#ddd;width:100%;margin-bottom:1rem;vertical-align:top;border-color:#dee2e6}
.table>:not(caption)>*>*{padding:.5rem .5rem;background-color:var(--bs-table-bg);border-bottom-width:1px;box-shadow:inset 0 0 0 9999px var(--bs-table-accent-bg)}
.table>tbody{vertical-align:inherit}
.table>thead{vertical-align:bottom}
.table-striped>tbody>tr:nth-of-type(odd)>*{--bs-table-accent-bg:var(--bs-table-striped-bg);color:var(--bs-table-striped-color)}
.table-sm>:not(caption)>*>*{padding:.25rem .25rem}
.table-bordered>:not(caption)>*{border-width:1px 0}
.table-bordered>:not(caption)>*>*{border-width:0 1px}
.table-hover>tbody>tr:hover>*{--bs-table-accent-bg:var(--bs-table-hover-bg);color:#212529}
.nav{display:flex;flex-wrap:wrap;padding-left:0;margin-bottom:0;list-style:none}
.nav-tabs{border-bottom:1px solid #dee2e6}
.nav-tabs .nav-link{margin-bottom:-1px;border:1px solid transparent;border-top-left-radius:.25rem;border-top-right-radius:.25rem}
.nav-link{display:block;padding:.5rem 1rem;text-decoration:none;color:#0d6efd}
.nav-link:hover{color:#0a58ca}
.nav-tabs .nav-link:hover{border-color:#e9ecef #e9ecef #dee2e6}
.nav-tabs .nav-link.active{color:#495057;background-color:#fff;border-color:#dee2e6 #dee2e6 #fff}
.tab-content>.tab-pane{display:none}
.tab-content>.active{display:block}
.badge{display:inline-block;padding:.35em .65em;font-size:.75em;font-weight:700;line-height:1;color:#fff;text-align:center;white-space:nowrap;vertical-align:baseline;border-radius:.25rem}
.bg-success{background-color:#198754!important}
.bg-danger{background-color:#dc3545!important}
.bg-warning{background-color:#ffc107!important;color:#000!important}
.bg-secondary{background-color:#6c757d!important}
.text-muted{color:#6c757d!important}
.mt-3{margin-top:1rem!important}
.mt-4{margin-top:1.5rem!important}
.mb-3{margin-bottom:1rem!important}
.p-3{padding:1rem!important}
.px-3{padding-right:1rem!important;padding-left:1rem!important}
.py-2{padding-top:.5rem!important;padding-bottom:.5rem!important}
.fw-bold{font-weight:700!important}
.text-end{text-align:right!important}
.alert{position:relative;padding:1rem 1rem;margin-bottom:1rem;border:1px solid transparent;border-radius:.25rem}
.alert-info{color:#055160;background-color:#cff4fc;border-color:#b6effb}
th[data-sort]{cursor:pointer;user-select:none}
th[data-sort]:after{content:" ⇕";opacity:.4}
th[data-sort].asc:after{content:" ↑";opacity:1}
th[data-sort].desc:after{content:" ↓";opacity:1}
.pass{background-color:#d1e7dd!important}
.fail{background-color:#f8d7da!important}
.warn{background-color:#fff3cd!important}
"""

SORT_JS = """
function sortTable(th) {
    var table = th.closest('table');
    var tbody = table.querySelector('tbody');
    var rows  = Array.from(tbody.querySelectorAll('tr'));
    var col   = Array.from(th.parentNode.children).indexOf(th);
    var asc   = th.classList.contains('asc');
    th.parentNode.querySelectorAll('th').forEach(function(h){h.classList.remove('asc','desc');});
    th.classList.add(asc ? 'desc' : 'asc');
    var dir = asc ? -1 : 1;
    rows.sort(function(a,b){
        var ta = a.children[col] ? a.children[col].textContent.trim() : '';
        var tb = b.children[col] ? b.children[col].textContent.trim() : '';
        var na = parseFloat(ta), nb = parseFloat(tb);
        if (!isNaN(na) && !isNaN(nb)) return dir*(na-nb);
        return dir*ta.localeCompare(tb);
    });
    rows.forEach(function(r){tbody.appendChild(r);});
}
function initTabs(){
    document.querySelectorAll('.nav-tabs .nav-link').forEach(function(link){
        link.addEventListener('click',function(e){
            e.preventDefault();
            var tabId=this.getAttribute('href');
            document.querySelectorAll('.nav-tabs .nav-link').forEach(function(l){l.classList.remove('active');});
            document.querySelectorAll('.tab-pane').forEach(function(p){p.classList.remove('active');});
            this.classList.add('active');
            document.querySelector(tabId).classList.add('active');
        });
    });
}
document.addEventListener('DOMContentLoaded',initTabs);
"""


def badge(value, thresholds=None, low_good=False):
    """
    Return an HTML badge for a numeric value.
    thresholds: (warn_val, fail_val) or None for no badge.
    low_good: True means lower is better (e.g. lambda %).
    """
    if not value or value in ("N/A", ""):
        return f'<span class="text-muted">—</span>'
    try:
        v = float(value)
    except ValueError:
        return html.escape(str(value))
    return html.escape(str(value))


def qc_badge_lambda(val):
    """Green badge if lambda meth% < 2, red otherwise."""
    if not val or val in ("", "N/A"):
        return '<span class="text-muted">—</span>'
    try:
        v = float(val)
        cls = "bg-success" if v < 2.0 else "bg-danger"
        return f'<span class="badge {cls}">{html.escape(val)}%</span>'
    except ValueError:
        return html.escape(val)


def qc_badge_puc19(val):
    """Green badge if pUC19 meth% > 90, red otherwise."""
    if not val or val in ("", "N/A"):
        return '<span class="text-muted">—</span>'
    try:
        v = float(val)
        cls = "bg-success" if v > 90.0 else "bg-danger"
        return f'<span class="badge {cls}">{html.escape(val)}%</span>'
    except ValueError:
        return html.escape(val)


def build_detail_table(sample_id, sample_data, metric_type, col_headers=None):
    """Build an HTML table for a single metric type for all samples."""
    # Collect all keys for this metric type across all samples
    all_keys = set()
    for sid, sdata in sample_data.items():
        all_keys.update(sdata.get(metric_type, {}).keys())

    # Sort keys by section then subsection then metric
    sorted_keys = sorted(all_keys, key=lambda k: (k[0], k[1], k[2]))

    if not sorted_keys:
        return '<p class="text-muted">No data available.</p>'

    samples = sorted(sample_data.keys())
    rows_html = []
    for section, sub, metric in sorted_keys:
        cells = [
            f'<td class="text-muted">{html.escape(section)}</td>',
            f'<td class="text-muted">{html.escape(sub)}</td>',
            f'<td>{html.escape(metric)}</td>',
        ]
        for sid in samples:
            d = sample_data.get(sid, {}).get(metric_type, {})
            val, pct = d.get((section, sub, metric), ("", ""))
            if pct:
                cell = f'{html.escape(str(val))} <span class="text-muted">({html.escape(str(pct))}%)</span>'
            else:
                cell = html.escape(str(val)) if val else '<span class="text-muted">—</span>'
            cells.append(f'<td class="text-end">{cell}</td>')
        rows_html.append(f'<tr>{"".join(cells)}</tr>')

    sample_headers = "".join(
        f'<th data-sort>{html.escape(sid)}</th>' for sid in samples
    )
    table_html = f"""
<div class="table-responsive">
<table class="table table-sm table-striped table-bordered table-hover">
<thead class="table-dark">
<tr>
  <th>Section</th><th>Subsection</th><th>Metric</th>
  {sample_headers}
</tr>
</thead>
<tbody>
{"".join(rows_html)}
</tbody>
</table>
</div>"""
    return table_html


def build_summary_table(all_data):
    """Build the summary tab HTML table."""
    samples = sorted(all_data.keys())
    if not samples:
        return '<p class="text-muted">No samples found.</p>'

    rows_data = [extract_summary(sid, all_data[sid]) for sid in samples]
    headers   = list(rows_data[0].keys())

    # Identify special columns for color-coding
    LAMBDA_COL = "lambda meth %"
    PUC19_COL  = "pUC19 meth %"

    thead = "<tr>" + "".join(
        f'<th data-sort onclick="sortTable(this)" style="white-space:nowrap">{html.escape(h)}</th>'
        if h == "Sample" else
        f'<th data-sort onclick="sortTable(this)">{html.escape(h)}</th>'
        for h in headers
    ) + "</tr>"

    tbody_rows = []
    for row in rows_data:
        cells = []
        for h in headers:
            val = row.get(h, "")
            if h == "Sample":
                cell = html.escape(str(val)) if val else '<span class="text-muted">—</span>'
                cells.append(f'<td style="white-space:nowrap">{cell}</td>')
            elif h == LAMBDA_COL:
                cells.append(f'<td class="text-end">{qc_badge_lambda(val)}</td>')
            elif h == PUC19_COL:
                cells.append(f'<td class="text-end">{qc_badge_puc19(val)}</td>')
            else:
                cell = html.escape(str(val)) if val else '<span class="text-muted">—</span>'
                cells.append(f'<td class="text-end">{cell}</td>')
        tbody_rows.append(f'<tr>{"".join(cells)}</tr>')

    return f"""
<div class="table-responsive">
<table class="table table-sm table-striped table-bordered table-hover">
<thead class="table-dark">
<thead>{thead}</thead>
<tbody>{"".join(tbody_rows)}</tbody>
</table>
</div>"""


def build_html(all_data, run_date=""):
    """Generate the full self-contained HTML report."""
    samples   = sorted(all_data.keys())
    n_samples = len(samples)

    tab_defs = [
        ("summary",  "Summary",     None),
        ("mapping",  "Mapping",     "mapping"),
        ("methyl",   "Methylation", "methyl"),
        ("wgs",      "WGS Coverage","wgs"),
        ("duplex",   "Duplex QC",   "duplex"),
        ("vc",       "Variants",    "vc"),
    ]

    nav_items = "".join(
        f'<li class="nav-item"><a class="nav-link{"  active" if i==0 else ""}" href="#{tab_id}">{tab_label}</a></li>'
        for i, (tab_id, tab_label, _) in enumerate(tab_defs)
    )

    tab_panes = []
    for i, (tab_id, tab_label, mtype) in enumerate(tab_defs):
        active = " active" if i == 0 else ""
        if mtype is None:
            content = build_summary_table(all_data)
        else:
            content = build_detail_table(tab_id, all_data, mtype)
        tab_panes.append(
            f'<div class="tab-pane{active}" id="{tab_id}"><div class="mt-3">{content}</div></div>'
        )

    pane_html = "\n".join(tab_panes)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TAPS QC Report — {n_samples} sample(s)</title>
<style>
{BOOTSTRAP_CSS}
</style>
</head>
<body>
<div class="container mt-4">
  <h1>TAPS QC Report</h1>
  <p class="text-muted">
    {n_samples} sample(s) &nbsp;|&nbsp; {run_date}
  </p>
  <div class="alert alert-info">
    <strong>QC thresholds:</strong>
    Lambda methylation (negCtrl) &lt; 2% → <span class="badge bg-success">PASS</span> &nbsp;|&nbsp;
    pUC19 methylation (posCtrl)  &gt; 90% → <span class="badge bg-success">PASS</span>
  </div>

  <ul class="nav nav-tabs mb-3">
    {nav_items}
  </ul>
  <div class="tab-content px-3">
    {pane_html}
  </div>
</div>
<script>
{SORT_JS}
</script>
</body>
</html>"""


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="TAPS QC batch HTML report generator")
    parser.add_argument("csvs", nargs="*", help="CSV files to include (auto-detected if omitted)")
    parser.add_argument("--out", default="taps_qc_report.html", help="Output HTML file")
    parser.add_argument("--tsv", default="taps_qc_report.tsv",  help="Output TSV file")
    args = parser.parse_args()

    csv_files = args.csvs
    if not csv_files:
        # Scan current directory
        csv_files = [str(p) for p in Path(".").glob("*.csv")]

    if not csv_files:
        print("[taps_qc_report] Warning: no CSV files found.", file=sys.stderr)

    all_data = load_all_csvs(csv_files)

    import datetime
    run_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    html_content = build_html(all_data, run_date=run_date)
    with open(args.out, "w") as fh:
        fh.write(html_content)
    print(f"[taps_qc_report] HTML report written to {args.out} ({len(all_data)} samples)", file=sys.stderr)

    write_tsv(all_data, args.tsv)
    print(f"[taps_qc_report] TSV summary written to {args.tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
