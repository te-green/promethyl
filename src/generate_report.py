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
from pathlib import Path

import pandas as pd
from jinja2 import Template

from panelapp import annotate_genes, CONFIDENCE_RANK

TAG_CLASS = {"Green": "tag-green", "Amber": "tag-amber", "Red": "tag-red", "Unknown": "tag-unknown"}

TEMPLATE = Template(r"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css">
<style>
  body { font-family: -apple-system, "Segoe UI", Roboto, sans-serif; margin: 24px; color: #1a1a1a; }
  h1 { font-size: 20px; margin-bottom: 4px; }
  .subtitle { color: #666; margin-bottom: 18px; font-size: 13px; }
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
  .gene-cell   { max-width: 280px; white-space: normal; }
  .dmr-cell    { max-width: 260px; white-space: normal; }
  tr.filter-row th { padding: 2px 4px; font-weight: normal; }
  tr.filter-row input, tr.filter-row select {
    width: 100%; box-sizing: border-box; font-size: 11px; padding: 2px 3px;
  }
  .clear-btn {
    font-size: 11px; padding: 2px 9px; border-radius: 10px; border: 1px solid #ccc;
    background: #fff; cursor: pointer; color: #444;
  }
  .clear-btn:hover { background: #f2f2f2; }
</style>
</head>
<body>
<h1>{{ title }}</h1>
<div class="subtitle">{{ n_rows }} rows &middot; generated {{ generated }}</div>
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
        {% else %}
        <input type="text" data-col="{{ loop.index0 }}" class="col-filter" placeholder="Filter...">
        {% endif %}
      </th>
      {% endfor %}
    </tr>
  </thead>
  <tbody>
    {% for row in rows %}
    <tr class="{{ 'outlier-row' if row['_outlier'] else '' }}">
      {% for col in columns %}
        {% if col == 'gene' %}<td class="gene-cell">{{ row['_gene_html'] }}</td>
        {% elif col == 'dmr_name' %}<td class="dmr-cell">{{ row['_dmr_html'] }}</td>
        {% elif col == 'disorder' %}{# folded into the dmr_name cell above #}
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
  var table = $('#report').DataTable({ pageLength: 25, order: [], orderCellsTop: true });

  $('.col-filter').on('change keyup', function () {
    var colIdx = $(this).data('col');
    var val = $(this).val();
    if ($(this).is('select')) {
      // exact match on the selected value, empty selection clears the filter
      var pattern = val ? '^' + $.fn.dataTable.util.escapeRegex(val) + '$' : '';
      table.column(colIdx).search(pattern, true, false).draw();
    } else {
      table.column(colIdx).search(val).draw();
    }
  });

  $('#clear-filters').on('click', function () {
    $('.col-filter').val('');
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
        best = hits[0]  # pre-sorted best-evidence-first
        title = f'{best["confidence"]} evidence'
        if len(hits) > 1:
            title += f" — +{len(hits) - 1} more panel(s)"
        parts.append(
            f'{sym} <span class="tag {TAG_CLASS.get(best["confidence"], "tag-unknown")}" '
            f'title="{title}">{best["panel"]}</span>'
        )
    return " ".join(parts)


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
    """For each displayed column, decide dropdown (low-cardinality -- chrom,
    tissue, outlier, etc.) vs free-text (everything else, including numeric
    columns and gene/dmr_name -- those render as HTML tags, not their raw
    value, so a dropdown of raw values wouldn't match what's on screen;
    free-text search still works there since DataTables searches rendered
    cell text, which usefully includes the tag labels).
    """
    filters = {}
    for col in columns:
        if col in ("gene", "dmr_name"):
            filters[col] = {"type": "text", "options": []}
            continue
        series = df[col].fillna("").astype(str) if col in df.columns else pd.Series(dtype=str)
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

    df = pd.read_csv(args.input, sep="\t")

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

    df["_gene_html"] = df.get("gene", "").fillna("").apply(lambda g: render_gene_cell(g, panelapp_cache))

    if "dmr_name" in df.columns:
        df["dmr_name"] = df["dmr_name"].fillna("")
        df["disorder"] = df["disorder"].fillna("") if "disorder" in df.columns else ""
        df["_dmr_html"] = df.apply(
            lambda r: render_dmr_cell(r["dmr_name"], r["disorder"]), axis=1
        )

    columns = [c for c in df.columns if not c.startswith("_") and c != "disorder"]
    column_filters = compute_column_filters(df, columns)
    html = TEMPLATE.render(
        title=args.title,
        n_rows=len(df),
        generated=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
        columns=columns,
        column_filters=column_filters,
        rows=df.to_dict(orient="records"),
    )
    Path(args.output).write_text(html)
    print(f"✓ Report written -> {args.output}")

if __name__ == "__main__":
    main()
