#!/usr/bin/env python3
"""
generate_report.py - Render a cohort outlier TSV (run_cohort.py output) as a
single self-contained, sortable/searchable/filterable HTML report, with each
gene tagged by its best PanelApp confidence level (Green/Amber/Red).

The report has three tabs: Results, Run configuration (cohort YAML + run info),
and Column descriptions.

Usage:
    python generate_report.py --input cohort_results.tsv --output report.html --config cohort.yaml
    python generate_report.py --input cohort_results.tsv --output report.html --config cohort.yaml --outliers-only
"""
import argparse
import platform
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from jinja2 import Template

from panelapp import annotate_genes, CONFIDENCE_RANK

TAG_CLASS = {"Green": "tag-green", "Amber": "tag-amber", "Red": "tag-red", "Unknown": "tag-unknown"}

# Displayed with reduced-precision formatting + a narrow right-aligned cell, full
# precision kept in the cell's title= tooltip. pval/padj use significant-figure
# formatting (.2g) rather than fixed decimals (.2f) -- a p-value of 3.16e-22
# rounded to 2 decimal *places* is just "0.00", which defeats the purpose.
NUMERIC_DISPLAY_COLS = {
    "coverage": "{:.2f}", "methylation": "{:.2f}", "group_mean": "{:.2f}",
    "delta": "{:.2f}", "zscore": "{:.2f}", "mod_score": "{:.2f}",
    "pval": "{:.2g}", "padj": "{:.2g}",
    "hap1_mod_score": "{:.2f}", "hap2_mod_score": "{:.2f}", "asm_pval": "{:.2g}",
}

# Sample/location/gene/region first (the "what and where" of a hit), everything
# else -- stats, coverage counts, cohort-wide summary columns -- follows in
# whatever order they already occur in. transcript is dropped entirely (gene_id
# covers it); disorder is folded into the dmr_name cell, never its own column.
COLUMN_PRIORITY = ["sample", "chrom", "start", "end", "gene", "cpg_island", "dmr_name", "gene_id"]
COLUMNS_DROPPED = {"transcript", "disorder"}
# Haplotype working columns stay in the TSV but out of the report: only
# hap1_mod_score, hap2_mod_score, asm_flag and asm_pval are displayed.
COLUMNS_DROPPED |= {"hap1_mod_count", "hap1_unmod_count", "hap2_mod_count", "hap2_unmod_count",
                    "n_hap_sites", "asm_delta", "asm_padj"}

# Descriptions shown in the "Column descriptions" tab. Edit freely -- these are
# first-draft wordings; any displayed column missing here is listed with a
# placeholder and reported on stderr so nothing silently goes undocumented.
COLUMN_DESCRIPTIONS = {
    "sample": "Sample ID.",
    "chrom": "Chromosome of the tested region.",
    "start": "Start coordinate of the tested region.",
    "end": "End coordinate of the tested region.",
    "gene": "Gene symbol(s) overlapping the region. Tags show the best PanelApp confidence "
            "(Green/Amber/Red) for each panel the gene is on; no tag = not in PanelApp.",
    "cpg_island": "CpG island (tested region) identifier / coordinates.",
    "dmr_name": "Known DMR overlapped by the region (purple tag), with its associated disorder (pink tag).",
    "gene_id": "Ensembl gene ID(s) overlapping the region, one per line.",
    "coverage": "Read / CpG-call coverage supporting the region for this sample.",
    "methylation": "Methylation level of the region in this sample (fraction, 0-1).",
    "group_mean": "Mean methylation of the comparison group (tissue- and sex-matched) for this region.",
    "delta": "methylation minus group_mean. Positive = hypermethylated (red row shading), "
             "negative = hypomethylated (blue row shading).",
    "zscore": "Z-score of the sample's methylation relative to the comparison group.",
    "mod_score": "Region-level modification score for this sample.",
    "pval": "p-value for the sample vs comparison group test.",
    "padj": "Multiple-testing-adjusted p-value.",
    "outlier": "TRUE if the region passed the outlier thresholds for this sample.",
    "any_outlier": "TRUE if the region was an outlier by any of the applied criteria.",
    "hap1_mod_score": "Modification score for haplotype 1 across the region.",
    "hap2_mod_score": "Modification score for haplotype 2 across the region.",
    "asm_flag": "TRUE if the region shows allele-specific methylation (significant difference between "
                "haplotypes, minimum coverage of 5 per haplotype).",
    "asm_pval": "p-value for the between-haplotype (chi-squared) test.",
}
MISSING_DESCRIPTION = "No description yet (add to COLUMN_DESCRIPTIONS in generate_report.py)."

TEMPLATE = Template(r"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css">
<style>
  body { font-family: -apple-system, "Segoe UI", Roboto, sans-serif; margin: 24px; color: #1a1a1a; }
  h1 { font-size: 20px; margin-bottom: 4px; }
  .subtitle { color: #666; margin-bottom: 12px; font-size: 13px; }
  .legend { margin-bottom: 14px; font-size: 12px; color: #444; }
  table.dataTable { font-size: 12.5px; }
  td, th { white-space: nowrap; }

  /* Tabs */
  .tab-bar { display: flex; gap: 4px; border-bottom: 1px solid #ccc; margin-bottom: 16px; }
  .tab-btn { font-size: 13px; padding: 7px 16px; border: 1px solid transparent; border-bottom: none;
             background: none; cursor: pointer; color: #555; border-radius: 6px 6px 0 0; }
  .tab-btn:hover { background: #f2f2f2; }
  .tab-btn.active { background: #fff; border-color: #ccc; color: #1a1a1a; font-weight: 600;
                    margin-bottom: -1px; border-bottom: 1px solid #fff; }
  .tab-pane { display: none; }
  .tab-pane.active { display: block; }

  /* Row shading by methylation state (outlier rows; sign of delta) */
  tr.row-hyper > td { background-color: #fdecee; }
  tr.row-hypo  > td { background-color: #e9f1fc; }
  tr.outlier-row > td { background-color: #fff7f0; }  /* fallback when no delta column */
  .swatch { display: inline-block; width: 12px; height: 12px; border: 1px solid #bbb;
            vertical-align: -2px; margin: 0 3px 0 0; }
  .swatch-hyper { background: #fdecee; }
  .swatch-hypo  { background: #e9f1fc; }

  .tag { display: inline-block; padding: 1px 7px; border-radius: 10px; font-size: 11px;
         font-weight: 600; color: #fff; margin: 1px 3px 1px 0; cursor: default; }
  .tag-green   { background: #2e7d32; }
  .tag-amber   { background: #e08e00; }
  .tag-red     { background: #c62828; }
  .tag-unknown { background: #9e9e9e; }
  .tag-dmr      { background: #5e35b1; }
  .tag-disorder { background: #ad1457; }

  /* One gene per line -- each gene + its panel tags is its own block. */
  .gene-cell   { max-width: 320px; white-space: normal; }
  .dmr-cell    { max-width: 260px; white-space: normal; }

  /* No max-width here -- scrollX (below) gives every column its natural width
     instead of columns fighting for space inside a fixed 100%-wide table. */
  /* gene_id: one ID per line, each a nowrap block. */
  .list-cell   { white-space: normal; }
  .list-cell .id-line { display: block; white-space: nowrap; }

  /* Reduced-precision numeric columns: narrow, right-aligned, full value on hover. */
  .num-cell { max-width: 64px; text-align: right; font-variant-numeric: tabular-nums; cursor: default; }

  tr.filter-row th { padding: 2px 4px; font-weight: normal; }
  tr.filter-row input, tr.filter-row select {
    width: 100%; box-sizing: border-box; font-size: 11px; padding: 2px 3px;
  }
  .clear-btn {
    font-size: 11px; padding: 2px 9px; border-radius: 10px; border: 1px solid #ccc;
    background: #fff; cursor: pointer; color: #444;
  }
  .clear-btn:hover { background: #f2f2f2; }

  .panel-filter-row { margin-bottom: 12px; font-size: 12.5px; display: flex; align-items: center; gap: 8px; }
  .panel-filter-row select { font-size: 12.5px; padding: 3px 6px; }

  /* Config + column description tabs */
  .info-table { border-collapse: collapse; font-size: 13px; margin-bottom: 20px; }
  .info-table th, .info-table td { text-align: left; vertical-align: top; padding: 5px 14px 5px 0;
                                   border-bottom: 1px solid #eee; white-space: normal; }
  .info-table th { font-weight: 600; color: #444; white-space: nowrap; }
  .info-table code, .config-pre { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12px; }
  .config-pre { background: #f6f8fa; border: 1px solid #e1e4e8; border-radius: 6px; padding: 12px 14px;
                overflow-x: auto; white-space: pre; max-width: 1100px; }
  .section-title { font-size: 14px; font-weight: 600; margin: 4px 0 8px; }
  .muted { color: #888; }
  .col-name { font-family: ui-monospace, Menlo, Consolas, monospace; font-size: 12px; white-space: nowrap; }
</style>
</head>
<body>
<h1>{{ title }}</h1>
<div class="subtitle">{{ n_rows }} rows &middot; generated {{ generated }}</div>

<div class="tab-bar">
  <button type="button" class="tab-btn active" data-tab="tab-results">Results</button>
  <button type="button" class="tab-btn" data-tab="tab-config">Run configuration</button>
  <button type="button" class="tab-btn" data-tab="tab-columns">Column descriptions</button>
</div>

<!-- ============================ TAB 1: RESULTS ============================ -->
<div id="tab-results" class="tab-pane active">

<div class="legend">
  PanelApp confidence:
  <span class="tag tag-green">Green</span> established &nbsp;
  <span class="tag tag-amber">Amber</span> moderate evidence &nbsp;
  <span class="tag tag-red">Red</span> limited evidence &nbsp;
  <span style="color:#888">no tag = not in PanelApp</span>
  &nbsp;|&nbsp;
  <span class="tag tag-dmr">DMR name</span>
  <span class="tag tag-disorder">associated disorder</span>
  (known DMRs, Table 1)
  &nbsp;|&nbsp;
  Row shading: <span class="swatch swatch-hypo"></span>hypomethylated
  <span class="swatch swatch-hyper" style="margin-left:8px"></span>hypermethylated
  {% if not shade_all %}<span class="muted">(outlier rows only)</span>{% endif %}
  &nbsp;|&nbsp;
  <button type="button" id="clear-filters" class="clear-btn">Clear filters</button>
</div>

{% if panel_counts %}
<div class="panel-filter-row">
  <label for="panel-select"><strong>Filter by gene panel:</strong></label>
  <select id="panel-select">
    <option value="">(All panels) &ndash; {{ n_rows }} rows</option>
    {% for panel, count in panel_counts %}<option value="{{ panel }}">{{ panel }} ({{ count }})</option>{% endfor %}
  </select>
</div>
{% endif %}

<table id="report" class="display" style="width:100%">
  <thead>
    <tr>{% for col in columns %}<th>{{ col }}</th>{% endfor %}</tr>
    <tr class="filter-row">
      {% for col in columns %}
      <th>
        {% if column_filters[col].type == 'select' %}
        <select data-col="{{ loop.index0 }}" class="col-filter">
          <option value="">(All)</option>
          {% for opt in column_filters[col].options %}<option value="{{ opt }}">{{ opt }}</option>{% endfor %}
        </select>
        {% elif column_filters[col].type == 'numeric' %}
        <input type="text" data-col="{{ loop.index0 }}" class="col-filter col-filter-numeric" placeholder="e.g. abs&gt;=.3">
        {% else %}
        <input type="text" data-col="{{ loop.index0 }}" class="col-filter" placeholder="Filter...">
        {% endif %}
      </th>
      {% endfor %}
    </tr>
  </thead>
  <tbody>
    {% for row in rows %}
    <tr class="{{ 'row-hyper' if row['_state'] == 'hyper' else ('row-hypo' if row['_state'] == 'hypo' else ('outlier-row' if row['_outlier'] and not has_delta else '')) }}" data-panels="{{ row['_panels'] }}">
      {% for col in columns %}
        {% if col == 'gene' %}<td class="gene-cell">{{ row['_gene_html'] }}</td>
        {% elif col == 'dmr_name' %}<td class="dmr-cell">{{ row['_dmr_html'] }}</td>
        {% elif col == 'gene_id' %}<td class="list-cell">{{ row['_id_html'] }}</td>
        {% elif col in numeric_display_cols %}<td class="num-cell" title="{{ row[col] }}">{{ row['_fmt_' + col] }}</td>
        {% else %}<td>{{ row[col] }}</td>{% endif %}
      {% endfor %}
    </tr>
    {% endfor %}
  </tbody>
</table>
</div>

<!-- ======================= TAB 2: RUN CONFIGURATION ======================= -->
<div id="tab-config" class="tab-pane">
  <div class="section-title">Report run</div>
  <table class="info-table">
    {% for k, v in run_info %}<tr><th>{{ k }}</th><td><code>{{ v | e }}</code></td></tr>
    {% endfor %}
  </table>

  <div class="section-title">Cohort configuration{% if config_path %} <span class="muted">({{ config_path | e }})</span>{% endif %}</div>
  {% if config_text %}
  <pre class="config-pre">{{ config_text | e }}</pre>
  {% else %}
  <div class="muted">No cohort YAML supplied &ndash; re-run generate_report.py with <code>--config cohort.yaml</code>.</div>
  {% endif %}
</div>

<!-- ======================= TAB 3: COLUMN DESCRIPTIONS ===================== -->
<div id="tab-columns" class="tab-pane">
  <table class="info-table">
    <tr><th>Column</th><th>Description</th></tr>
    {% for name, desc in column_docs %}
    <tr><td class="col-name">{{ name }}</td><td>{{ desc | e }}</td></tr>
    {% endfor %}
  </table>
</div>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
<script>
$(document).ready(function () {
  // scrollX/scrollY (DataTables' native scroll handling, header pinned via its
  // own internal split header/body tables) instead of a fixed 100%-wide table.
  var table = $('#report').DataTable({
    pageLength: 25, order: [], orderCellsTop: true,
    scrollX: true, scrollY: '600px', scrollCollapse: true
  });

  // Tabs. DataTables mis-sizes columns if it was laid out while hidden, so
  // re-adjust whenever the Results tab is shown again.
  $('.tab-btn').on('click', function () {
    var target = $(this).data('tab');
    $('.tab-btn').removeClass('active');
    $(this).addClass('active');
    $('.tab-pane').removeClass('active');
    $('#' + target).addClass('active');
    if (target === 'tab-results') { table.columns.adjust(); }
  });

  // Numeric column filters (coverage, delta, zscore, padj, etc.) support
  // >=, <=, >, <, = comparisons -- a bare number with no operator means exact match.
  // DataTables' built-in column().search() only does substring/regex text matching,
  // so numeric comparisons need a custom search plugin instead; it runs alongside
  // (ANDs with) the normal per-column searches used for select/text filters below.
  var numericFilters = {};

  function parseNumericFilter(raw) {
    // Optional leading "abs" makes the comparison apply to |value| -- useful
    // for delta/zscore, where a hit can be hyper- or hypomethylated.
    var m = String(raw).trim().match(/^(abs)?\s*(>=|<=|>|<|=)?\s*(-?\d+\.?\d*)$/i);
    if (!m) return null;
    return { abs: !!m[1], op: m[2] || '=', val: parseFloat(m[3]) };
  }

  $.fn.dataTable.ext.search.push(function (settings, rowData) {
    for (var colIdx in numericFilters) {
      var f = numericFilters[colIdx];
      var cell = parseFloat(rowData[colIdx]);
      if (isNaN(cell)) return false;
      if (f.abs) cell = Math.abs(cell);
      if (f.op === '>=' && !(cell >= f.val)) return false;
      if (f.op === '<=' && !(cell <= f.val)) return false;
      if (f.op === '>'  && !(cell >  f.val)) return false;
      if (f.op === '<'  && !(cell <  f.val)) return false;
      if (f.op === '='  && !(cell === f.val)) return false;
    }
    return true;
  });

  $('.col-filter').on('change keyup', function () {
    var colIdx = $(this).data('col');
    var val = $(this).val();
    if ($(this).is('select')) {
      // exact match on the selected value, empty selection clears the filter
      var pattern = val ? '^' + $.fn.dataTable.util.escapeRegex(val) + '$' : '';
      table.column(colIdx).search(pattern, true, false).draw();
    } else if ($(this).hasClass('col-filter-numeric')) {
      if (!val) {
        delete numericFilters[colIdx];
      } else {
        var parsed = parseNumericFilter(val);
        // invalid/incomplete input (e.g. still typing "-") drops any active
        // constraint for this column rather than filtering everything out
        if (parsed) { numericFilters[colIdx] = parsed; } else { delete numericFilters[colIdx]; }
      }
      table.draw();
    } else {
      table.column(colIdx).search(val).draw();
    }
  });

  // Gene-panel filter (single alphabetized dropdown).
  var activePanel = '';
  $.fn.dataTable.ext.search.push(function (settings, data, dataIndex) {
    if (!activePanel) return true;
    var panels = table.row(dataIndex).node().getAttribute('data-panels') || '';
    return panels.split(';').indexOf(activePanel) !== -1;
  });

  $('#panel-select').on('change', function () {
    activePanel = $(this).val();
    table.draw();
  });

  $('#clear-filters').on('click', function () {
    $('.col-filter').val('');
    numericFilters = {};
    activePanel = '';
    $('#panel-select').val('');
    table.columns().search('').draw();
  });
});
</script>
</body>
</html>
""")


def render_id_list(field: str) -> str:
    """One ID per line, each in its own nowrap span so an ID like
    ENSG00000239857.8 can never break mid-string."""
    ids = [i for i in (field or "").split(";") if i]
    return "".join(f'<span class="id-line">{i}</span>' for i in ids)


def split_symbols(gene_field: str) -> list[str]:
    """Gene column is ';'-joined and sometimes falls back to an Ensembl ID
    when annotate.py had no clean symbol -- skip those, PanelApp won't match them."""
    return [g for g in (gene_field or "").split(";") if g and not g.startswith("ENSG")]


def render_gene_cell(gene_field: str, panelapp_cache: dict) -> str:
    symbols = split_symbols(gene_field)
    if not symbols:
        return gene_field or ""
    parts = []
    for sym in symbols:
        hits = panelapp_cache.get(sym, [])
        if not hits:
            parts.append(sym)
            continue
        tags = "".join(
            f'<span class="tag {TAG_CLASS.get(h["confidence"], "tag-unknown")}" '
            f'title="{h["confidence"]} evidence">{h["panel"]}</span>'
            for h in hits
        )
        parts.append(f'{sym} {tags}')
    return "<br>".join(parts)


def row_panels(gene_field: str, panelapp_cache: dict) -> str:
    """';'-joined set of every distinct panel any gene in this row is on --
    drives both the panel filter and the row's data-panels attribute."""
    panels = set()
    for sym in split_symbols(gene_field):
        for hit in panelapp_cache.get(sym, []):
            panels.add(hit["panel"])
    return ";".join(sorted(panels))


def render_dmr_cell(dmr_field: str, disorder_field: str) -> str:
    """dmr_name and disorder are positionally paired ';'-joined lists (see
    annotate_dmrs()) -- zip them back together rather than rendering
    independently, so a disorder tag never ends up next to the wrong DMR."""
    names = [n for n in (dmr_field or "").split(";") if n]
    if not names:
        return ""
    disorders = (disorder_field or "").split(";")
    disorders += [""] * (len(names) - len(disorders))  # defensive padding
    parts = []
    for name, disorder in zip(names, disorders):
        parts.append(f'<span class="tag tag-dmr">{name}</span>')
        if disorder:
            parts.append(f'<span class="tag tag-disorder">{disorder}</span>')
    return " ".join(parts)


def methylation_state(delta, is_outlier: bool, shade_all: bool) -> str:
    """'hyper' / 'hypo' / '' from the sign of delta (methylation - group_mean).
    By default only outlier rows are shaded -- otherwise nearly every row in a
    full cohort TSV would carry a tint and the shading would stop meaning anything."""
    if pd.isna(delta) or not (is_outlier or shade_all):
        return ""
    if delta > 0:
        return "hyper"
    if delta < 0:
        return "hypo"
    return ""


def compute_column_filters(df: pd.DataFrame, columns: list[str], max_options: int = 20) -> dict:
    """For each displayed column, decide numeric (>=, <=, >, <, = comparisons --
    coverage_*, delta_*, zscore_*, padj_*, etc.), dropdown (low-cardinality
    non-numeric -- chrom, tissue, outlier, etc.), or free-text (everything
    else, including gene/dmr_name -- those render as HTML tags, not their raw
    value, so a dropdown of raw values wouldn't match what's on screen;
    free-text search still works there since DataTables searches rendered
    cell text, which usefully includes the tag labels).
    """
    filters = {}
    for col in columns:
        if col in ("gene", "dmr_name"):
            filters[col] = {"type": "text", "options": []}
            continue
        if col not in df.columns:
            filters[col] = {"type": "text", "options": []}
            continue
        # bool is numeric-dtype in pandas too (outlier flags) -- keep those as a
        # dropdown, not a >=/<= comparison, since TRUE/FALSE isn't ordinal here.
        if pd.api.types.is_numeric_dtype(df[col]) and not pd.api.types.is_bool_dtype(df[col]):
            filters[col] = {"type": "numeric", "options": []}
            continue
        series = df[col].fillna("").astype(str)
        uniques = sorted(v for v in series.unique() if v != "")
        if 0 < len(uniques) <= max_options:
            filters[col] = {"type": "select", "options": uniques}
        else:
            filters[col] = {"type": "text", "options": []}
    return filters


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Cohort outlier TSV (run_cohort.py output)")
    p.add_argument("--output", required=True, help="Output HTML report path")
    p.add_argument("--config", default=None,
                    help="Cohort YAML used for the run; embedded verbatim in the 'Run configuration' tab")
    p.add_argument("--panelapp-cache", default="panelapp_cache.json",
                    help="Local JSON cache of PanelApp lookups (reused/extended across runs)")
    p.add_argument("--title", default="promethyl cohort outlier report")
    p.add_argument("--outliers-only", action="store_true",
                    help="Only include rows flagged as an outlier")
    p.add_argument("--shade-all", action="store_true",
                    help="Shade hyper/hypomethylated state on every row, not just outlier rows")
    args = p.parse_args()

    # dtype=str on the text-y columns stops pandas from inferring a numeric/
    # mixed dtype when a handful of rows have a blank dmr_name/disorder/etc
    # (that's what the "Columns (0,1) have mixed types" warning was about).
    # Any of these that are still NaN after this (fully-empty column, or a
    # column not present in this input) get fillna'd below before templating.
    df = pd.read_csv(
        args.input, sep="\t",
        dtype={"gene": str, "gene_id": str, "dmr_name": str, "disorder": str},
    )
    n_input_rows = len(df)

    outlier_col = "outlier" if "outlier" in df.columns else ("any_outlier" if "any_outlier" in df.columns else None)
    if outlier_col:
        df["_outlier"] = df[outlier_col].astype(str).str.upper() == "TRUE"
        if args.outliers_only:
            df = df[df["_outlier"]]
    else:
        df["_outlier"] = False

    all_symbols = set()
    for field in df.get("gene", pd.Series(dtype=str)).fillna(""):
        all_symbols.update(split_symbols(field))

    panelapp_cache = annotate_genes(all_symbols, Path(args.panelapp_cache))

    # Every text-y column that can legitimately be blank needs fillna("") here --
    # an empty string written to TSV comes back as NaN on this fresh read, and
    # dtype-detection-based fillna misses pandas' StringDtype. Just fillna all
    # of them unconditionally; it's dtype-agnostic and version-proof.
    for col in ("gene", "gene_id", "dmr_name", "disorder"):
        if col in df.columns:
            df[col] = df[col].fillna("")

    df["_gene_html"] = df.get("gene", "").apply(lambda g: render_gene_cell(g, panelapp_cache))
    df["_panels"] = df.get("gene", "").apply(lambda g: row_panels(g, panelapp_cache))
    df["_id_html"] = df.get("gene_id", "").apply(render_id_list)

    if "dmr_name" in df.columns:
        if df.empty:
            df["_dmr_html"] = ""
        else:
            df["_dmr_html"] = df.apply(
                lambda r: render_dmr_cell(r["dmr_name"], r.get("disorder", "")), axis=1
            )

    # Row shading by methylation state (sign of delta).
    has_delta = "delta" in df.columns
    if has_delta:
        df["_state"] = [
            methylation_state(d, o, args.shade_all) for d, o in zip(df["delta"], df["_outlier"])
        ]
    else:
        df["_state"] = ""

    # Panel dropdown: alphabetized (case-insensitive) with per-panel row counts.
    panel_counter = Counter()
    for panels_str in df["_panels"]:
        for panel in panels_str.split(";"):
            if panel:
                panel_counter[panel] += 1
    panel_counts = sorted(panel_counter.items(), key=lambda kv: kv[0].lower())

    for col, fmt in NUMERIC_DISPLAY_COLS.items():
        if col in df.columns:
            df[f"_fmt_{col}"] = df[col].apply(lambda v, fmt=fmt: "" if pd.isna(v) else fmt.format(v))

    columns = [c for c in df.columns if not c.startswith("_") and c not in COLUMNS_DROPPED]
    columns = [c for c in COLUMN_PRIORITY if c in columns] + [c for c in columns if c not in COLUMN_PRIORITY]

    # Column descriptions tab: displayed columns, in displayed order.
    column_docs = [(c, COLUMN_DESCRIPTIONS.get(c, MISSING_DESCRIPTION)) for c in columns]
    undocumented = [c for c in columns if c not in COLUMN_DESCRIPTIONS]
    if undocumented:
        print(f"! No description for column(s): {', '.join(undocumented)}", file=sys.stderr)

    # Run configuration tab.
    config_text, config_path = "", None
    if args.config:
        config_path = str(Path(args.config).resolve())
        try:
            config_text = Path(args.config).read_text()
        except OSError as e:
            print(f"! Could not read --config {args.config}: {e}", file=sys.stderr)
    run_info = [
        ("Command", " ".join([Path(sys.argv[0]).name] + sys.argv[1:])),
        ("Input TSV", str(Path(args.input).resolve())),
        ("Rows in input", str(n_input_rows)),
        ("Rows displayed", f"{len(df)}" + (" (outliers only)" if args.outliers_only else "")),
        ("Row shading", "all rows" if args.shade_all else "outlier rows only"),
        ("PanelApp cache", str(Path(args.panelapp_cache).resolve())),
        ("Generated", pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")),
        ("Host", platform.node()),
        ("Python / pandas", f"{platform.python_version()} / {pd.__version__}"),
    ]

    column_filters = compute_column_filters(df, columns)
    html = TEMPLATE.render(
        title=args.title,
        n_rows=len(df),
        generated=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
        columns=columns,
        column_filters=column_filters,
        numeric_display_cols=set(NUMERIC_DISPLAY_COLS),
        panel_counts=panel_counts,
        rows=df.to_dict(orient="records"),
        has_delta=has_delta,
        shade_all=args.shade_all,
        run_info=run_info,
        config_text=config_text,
        config_path=config_path,
        column_docs=column_docs,
    )
    Path(args.output).write_text(html)
    print(f"✓ Report written -> {args.output}")


if __name__ == "__main__":
    main()