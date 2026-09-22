#!/usr/bin/env python3
"""
generate_report.py - Render a cohort outlier TSV (run_cohort.py output) as a
single self-contained, sortable/searchable/filterable HTML report, with each
gene tagged by its best PanelApp confidence level (Green/Amber/Red).

Usage:
    python generate_report.py --input cohort_results.tsv --output report.html
    python generate_report.py --input cohort_results.tsv --output report.html --outliers-only
"""
import argparse
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
    "delta": "{:.2f}", "zscore": "{:.2f}",
    "pval": "{:.2g}", "padj": "{:.2g}",
}

# Sample/location/gene/region first (the "what and where" of a hit), everything
# else -- stats, coverage counts, cohort-wide summary columns -- follows in
# whatever order they already occur in. transcript is dropped entirely (gene_id
# covers it); disorder is folded into the dmr_name cell, never its own column.
COLUMN_PRIORITY = ["sample", "chrom", "start", "end", "gene", "cpg_island", "dmr_name", "gene_id"]
COLUMNS_DROPPED = {"transcript", "disorder"}

TEMPLATE = Template(r"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css">
<style>
  body { font-family: -apple-system, "Segoe UI", Roboto, sans-serif; margin: 24px; color: #1a1a1a; }
  h1 { font-size: 20px; margin-bottom: 4px; }
  .subtitle { color: #666; margin-bottom: 4px; font-size: 13px; }
  .legend { margin-bottom: 14px; font-size: 12px; color: #444; }
  table.dataTable { font-size: 12.5px; }
  td, th { white-space: nowrap; }
  tr.outlier-row { background: #fff7f0; }
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
     instead of columns fighting for space inside a fixed 100%-wide table, which
     was forcing gene_id to wrap even when this had a generous max-width. */
  .list-cell   { white-space: normal; word-break: break-word; }

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

  .panel-summary { margin-bottom: 16px; }
  .panel-summary-title { font-size: 12.5px; font-weight: 600; color: #333; margin-bottom: 6px; }
  .panel-pill {
    display: inline-flex; align-items: center; gap: 6px; background: #f0f2f5; border: 1px solid #dfe3e8;
    border-radius: 14px; padding: 4px 12px; margin: 0 6px 6px 0; font-size: 12px; cursor: pointer; color: #333;
  }
  .panel-pill:hover { background: #e4e8ed; }
  .panel-pill.active { background: #1a1a1a; color: #fff; border-color: #1a1a1a; }
  .panel-pill .count { font-weight: 700; }

  .panel-filter-row { margin-bottom: 12px; font-size: 12.5px; display: flex; align-items: center; gap: 8px; }
  .panel-filter-row select { font-size: 12.5px; padding: 3px 6px; }
</style>
</head>
<body>
<h1>{{ title }}</h1>
<div class="subtitle">{{ n_rows }} rows &middot; generated {{ generated }}</div>

{% if panel_counts %}
<div class="panel-summary">
  <div class="panel-summary-title">Panels with hits</div>
  <span class="panel-pill active" data-panel="">All panels <span class="count">({{ n_rows }})</span></span>
  {% for panel, count in panel_counts %}
  <span class="panel-pill" data-panel="{{ panel }}">{{ panel }} <span class="count">({{ count }})</span></span>
  {% endfor %}
</div>
{% endif %}

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
  <button type="button" id="clear-filters" class="clear-btn">Clear filters</button>
</div>

{% if panel_counts %}
<div class="panel-filter-row">
  <label for="panel-select"><strong>Filter by gene panel:</strong></label>
  <select id="panel-select">
    <option value="">(All panels)</option>
    {% for panel, count in panel_counts %}<option value="{{ panel }}">{{ panel }}</option>{% endfor %}
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
        <input type="text" data-col="{{ loop.index0 }}" class="col-filter col-filter-numeric" placeholder="e.g. &gt;=15">
        {% else %}
        <input type="text" data-col="{{ loop.index0 }}" class="col-filter" placeholder="Filter...">
        {% endif %}
      </th>
      {% endfor %}
    </tr>
  </thead>
  <tbody>
    {% for row in rows %}
    <tr class="{{ 'outlier-row' if row['_outlier'] else '' }}" data-panels="{{ row['_panels'] }}">
      {% for col in columns %}
        {% if col == 'gene' %}<td class="gene-cell">{{ row['_gene_html'] }}</td>
        {% elif col == 'dmr_name' %}<td class="dmr-cell">{{ row['_dmr_html'] }}</td>
        {% elif col == 'gene_id' %}<td class="list-cell">{{ row[col].replace(';', ';<br>') | safe }}</td>
        {% elif col in numeric_display_cols %}<td class="num-cell" title="{{ row[col] }}">{{ row['_fmt_' + col] }}</td>
        {% else %}<td>{{ row[col] }}</td>{% endif %}
      {% endfor %}
    </tr>
    {% endfor %}
  </tbody>
</table>
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
<script>
$(document).ready(function () {
  // scrollX/scrollY (DataTables' native scroll handling, header pinned via its
  // own internal split header/body tables) instead of a fixed 100%-wide table --
  // this is also what stops narrower text columns like gene_id from being
  // squeezed by the other ~15 columns competing for a fixed viewport width.
  var table = $('#report').DataTable({
    pageLength: 25, order: [], orderCellsTop: true,
    scrollX: true, scrollY: '600px', scrollCollapse: true
  });

  // Numeric column filters (coverage, delta, zscore, padj, etc.) support
  // >=, <=, >, <, = comparisons -- a bare number with no operator means exact match.
  // DataTables' built-in column().search() only does substring/regex text matching,
  // so numeric comparisons need a custom search plugin instead; it runs alongside
  // (ANDs with) the normal per-column searches used for select/text filters below.
  var numericFilters = {};

  function parseNumericFilter(raw) {
    var m = String(raw).trim().match(/^(>=|<=|>|<|=)?\s*(-?\d+\.?\d*)$/);
    if (!m) return null;
    return { op: m[1] || '=', val: parseFloat(m[2]) };
  }

  $.fn.dataTable.ext.search.push(function (settings, rowData) {
    for (var colIdx in numericFilters) {
      var f = numericFilters[colIdx];
      var cell = parseFloat(rowData[colIdx]);
      if (isNaN(cell)) return false;
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

  // Gene-panel filter -- separate from the per-column filter row above, driven
  // by either the dropdown or clicking a summary pill; both set the same filter.
  var activePanel = '';
  $.fn.dataTable.ext.search.push(function (settings, data, dataIndex) {
    if (!activePanel) return true;
    var panels = table.row(dataIndex).node().getAttribute('data-panels') || '';
    return panels.split(';').indexOf(activePanel) !== -1;
  });

  $('#panel-select').on('change', function () {
    activePanel = $(this).val();
    $('.panel-pill').removeClass('active').filter(function () {
      return $(this).data('panel') === activePanel;
    }).addClass('active');
    table.draw();
  });

  $('.panel-pill').on('click', function () {
    activePanel = $(this).data('panel') || '';
    $('#panel-select').val(activePanel);
    $('.panel-pill').removeClass('active');
    $(this).addClass('active');
    table.draw();
  });

  $('#clear-filters').on('click', function () {
    $('.col-filter').val('');
    numericFilters = {};
    activePanel = '';
    $('#panel-select').val('');
    $('.panel-pill').removeClass('active');
    $('.panel-pill[data-panel=""]').addClass('active');
    table.columns().search('').draw();
  });
});
</script>
</body>
</html>
""")


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
    p.add_argument("--panelapp-cache", default="panelapp_cache.json",
                    help="Local JSON cache of PanelApp lookups (reused/extended across runs)")
    p.add_argument("--title", default="promethyl cohort outlier report")
    p.add_argument("--outliers-only", action="store_true",
                    help="Only include rows flagged as an outlier")
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

    # Every text-y column that can legitimately be blank (most rows have no gene,
    # dmr_name, disorder or gene_id overlap) needs fillna("") here -- the writer
    # already does this in memory, but an empty string written to TSV comes back
    # as NaN on this fresh read, regardless of the dtype= hint above: dtype=str
    # produces pandas' StringDtype here (pandas 3.x), and
    # pd.api.types.is_object_dtype() -- used by an earlier version of this dtype
    # check -- returns False for StringDtype, so a dtype-detection-based fillna
    # silently misses exactly the columns dtype=str was applied to. Just fillna
    # all of them unconditionally instead; it's dtype-agnostic and version-proof.
    for col in ("gene", "gene_id", "dmr_name", "disorder"):
        if col in df.columns:
            df[col] = df[col].fillna("")

    df["_gene_html"] = df.get("gene", "").apply(lambda g: render_gene_cell(g, panelapp_cache))
    df["_panels"] = df.get("gene", "").apply(lambda g: row_panels(g, panelapp_cache))

    if "dmr_name" in df.columns:
        df["_dmr_html"] = df.apply(
            lambda r: render_dmr_cell(r["dmr_name"], r.get("disorder", "")), axis=1
        )

    panel_counter = Counter()
    for panels_str in df["_panels"]:
        for panel in panels_str.split(";"):
            if panel:
                panel_counter[panel] += 1
    panel_counts = sorted(panel_counter.items(), key=lambda kv: (-kv[1], kv[0]))

    for col, fmt in NUMERIC_DISPLAY_COLS.items():
        if col in df.columns:
            df[f"_fmt_{col}"] = df[col].apply(lambda v, fmt=fmt: "" if pd.isna(v) else fmt.format(v))

    columns = [c for c in df.columns if not c.startswith("_") and c not in COLUMNS_DROPPED]
    columns = [c for c in COLUMN_PRIORITY if c in columns] + [c for c in columns if c not in COLUMN_PRIORITY]

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
    )
    Path(args.output).write_text(html)
    print(f"✓ Report written -> {args.output}")

if __name__ == "__main__":
    main()
