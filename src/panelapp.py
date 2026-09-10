"""
panelapp.py - Annotate gene symbols with PanelApp Australia panel membership
and confidence level (Green/Amber/Red), with local JSON caching so a cohort
run doesn't re-query the API for every sample x island x gene.

API: https://panelapp-aus.org/api/v1/genes/?entity_name=<symbol>
NOTE: I could not confirm from outside the HPC network whether the
`entity_name` server-side filter reliably narrows results on every
deployment/version -- as a safety net, results are always re-checked
client-side against gene_symbol/entity_name/alias before being accepted,
so a non-filtering server just costs an extra local comparison rather than
silently contaminating a gene's tags with unrelated genes.

Verify connectivity from wherever this actually runs before relying on it:
HPC compute nodes often have no outbound internet access even when the
login node does.
"""
import json
import time
from pathlib import Path

import requests

PANELAPP_BASE = "https://panelapp-aus.org/api/v1/genes/"

# PanelApp confidence_level: "3"=Green, "2"=Amber, "1"=Red, "0"=other/removed.
CONFIDENCE_LABELS = {"3": "Green", "2": "Amber", "1": "Red", "0": "Red"}
CONFIDENCE_RANK = {"Green": 3, "Amber": 2, "Red": 1, "Unknown": 0}


def _matches_symbol(record: dict, symbol: str) -> bool:
    gd = record.get("gene_data") or {}
    candidates = {record.get("entity_name"), gd.get("gene_symbol"), gd.get("hgnc_symbol")}
    candidates |= set(gd.get("alias") or [])
    return symbol.upper() in {c.upper() for c in candidates if c}


def _query_gene(symbol: str, session: requests.Session, timeout: int = 15) -> list[dict]:
    """Query PanelApp for one gene symbol; returns raw per-panel hit records."""
    resp = session.get(
        PANELAPP_BASE,
        params={"entity_name": symbol, "format": "json"},
        timeout=timeout,
    )
    resp.raise_for_status()
    results = resp.json().get("results", [])
    return [r for r in results if _matches_symbol(r, symbol)]


def annotate_gene(symbol: str, session: requests.Session) -> list[dict]:
    """Return this gene's panel hits as [{panel, confidence, moi, phenotypes}, ...],
    sorted best-evidence-first."""
    hits = []
    for h in _query_gene(symbol, session):
        conf = CONFIDENCE_LABELS.get(str(h.get("confidence_level")), "Unknown")
        hits.append({
            "panel": (h.get("panel") or {}).get("name", "?"),
            "confidence": conf,
            "moi": h.get("mode_of_inheritance", "") or "",
            "phenotypes": h.get("phenotypes", []) or [],
        })
    hits.sort(key=lambda x: -CONFIDENCE_RANK.get(x["confidence"], 0))
    return hits


def bulk_download_all_genes(session: requests.Session, page_size: int = 1000, sleep: float = 0.34) -> dict:
    """Walk every page of /api/v1/genes/ once and build a complete
    symbol -> [panel hits] index locally, instead of one filtered request
    per gene symbol.

    This is a better fit than annotate_genes()'s per-gene queries once
    your gene list is more than a couple dozen: request COUNT here is
    bounded by the size of PanelApp's whole database (fixed, one-time),
    not by how many distinct genes show up in your cohort's outliers
    (which only grows). It also sidesteps a real uncertainty in
    annotate_gene(): I was not able to confirm from outside your network
    that the `entity_name=` server-side filter reliably narrows results
    -- this walks the unfiltered listing and matches every record
    client-side by gene_symbol, so there's nothing to trust there.

    Safe to re-run occasionally to pick up PanelApp updates (new genes
    added to panels, confidence changes) -- not needed on every pipeline
    run once you've populated the cache once.
    """
    index: dict[str, list[dict]] = {}
    url = PANELAPP_BASE
    params = {"format": "json", "page_size": page_size}
    page_n = 0
    total = None
    while url:
        resp = session.get(url, params=params if page_n == 0 else None, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        total = total if total is not None else data.get("count")

        for rec in data.get("results", []):
            gd = rec.get("gene_data") or {}
            symbol = gd.get("gene_symbol") or rec.get("entity_name")
            if not symbol:
                continue
            conf = CONFIDENCE_LABELS.get(str(rec.get("confidence_level")), "Unknown")
            index.setdefault(symbol, []).append({
                "panel": (rec.get("panel") or {}).get("name", "?"),
                "confidence": conf,
                "moi": rec.get("mode_of_inheritance", "") or "",
                "phenotypes": rec.get("phenotypes", []) or [],
            })

        page_n += 1
        seen = min(page_n * page_size, total) if total else page_n * page_size
        print(f"  bulk PanelApp download: page {page_n} ({seen:,}/{total or '?'} records, "
              f"{len(index):,} symbols so far)")

        url = data.get("next")
        params = None  # page params are already folded into `next`
        time.sleep(sleep)

    for sym in index:
        index[sym].sort(key=lambda x: -CONFIDENCE_RANK.get(x["confidence"], 0))
    print(f"  PanelApp bulk download complete: {len(index):,} unique gene symbols indexed")
    return index


def load_cache(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def save_cache(path: Path, cache: dict) -> None:
    path.write_text(json.dumps(cache, indent=2, sort_keys=True))


def annotate_genes(symbols: set[str], cache_path: Path, sleep: float = 0.34) -> dict:
    """Annotate a set of gene symbols against PanelApp, reusing/updating a
    local JSON cache keyed by symbol. `sleep` throttles requests (~3 req/s
    default -- polite for a shared public service, no published hard limit)."""
    cache = load_cache(cache_path)
    session = requests.Session()
    n_new = 0
    for sym in sorted(symbols):
        if sym in cache:
            continue
        try:
            cache[sym] = annotate_gene(sym, session)
        except requests.RequestException as e:
            print(f"  WARNING: PanelApp lookup failed for {sym}: {e}")
            cache[sym] = []
        n_new += 1
        if n_new % 20 == 0:
            save_cache(cache_path, cache)
        time.sleep(sleep)
    save_cache(cache_path, cache)
    print(f"  PanelApp: {n_new} new lookups this run, {len(cache)} genes cached total")
    return cache

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="One-time bulk download of the entire PanelApp gene database into a "
                     "local cache, so generate_report.py never needs a live per-gene lookup "
                     "for genes already covered by it.")
    p.add_argument("--cache", default="panelapp_cache.json",
                    help="Cache file to write (merged with any existing entries)")
    p.add_argument("--page-size", type=int, default=1000)
    args = p.parse_args()

    cache_path = Path(args.cache)
    existing = load_cache(cache_path)
    bulk = bulk_download_all_genes(requests.Session(), page_size=args.page_size)
    existing.update(bulk)  # bulk data wins on overlap -- it's the freshest, complete pull
    save_cache(cache_path, existing)
    print(f"  Cache written -> {cache_path} ({len(existing):,} total symbols)")
