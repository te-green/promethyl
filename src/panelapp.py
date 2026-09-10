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
