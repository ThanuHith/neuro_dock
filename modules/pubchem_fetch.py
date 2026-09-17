# =============================================================================
# modules/pubchem_fetch.py
# =============================================================================
# Automatically fetches SMILES strings from PubChem by compound name.
#
# Why this matters for publication:
#   Using PubChem's PUG REST API means every compound has a traceable
#   PubChem CID (Compound ID) that reviewers can verify independently.
#   This is far more credible than hand-typed SMILES with no source.
#
# API used: PubChem PUG REST
#   https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<name>/property/...
#
# Rate limit: PubChem allows ~5 requests/second — we add a small delay
# to stay well within limits and avoid being blocked.
# =============================================================================

import logging
import time
import json
from pathlib import Path
from typing import List, Dict, Optional

log = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_DELAY = 0.3   # seconds between requests — stays under PubChem's rate limit


def fetch_compound_by_name(name: str) -> Optional[Dict]:
    """
    Fetch CID, canonical SMILES, and molecular weight for a compound by name.

    Args:
        name : compound name (e.g. "Quercetin", "Rasagiline")

    Returns:
        Dict with keys: name, cid, smiles, mw  — or None if not found
    """
    import requests
    from urllib.parse import quote

    encoded_name = quote(name)
    url = (
        f"{PUBCHEM_BASE}/compound/name/{encoded_name}/property/"
        f"IsomericSMILES,ConnectivitySMILES,MolecularWeight,IUPACName/JSON"
    )

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/120.0.0.0 Safari/537.36"
        )
    }

    try:
        resp = requests.get(url, headers=headers, timeout=15)
    except Exception as e:
        log.error(f"  [pubchem] Request failed for '{name}': {e}")
        return None

    if resp.status_code != 200:
        log.warning(f"  [pubchem] '{name}' not found on PubChem (HTTP {resp.status_code})")
        return None

    try:
        data = resp.json()
        props = data["PropertyTable"]["Properties"][0]
        # Prefer IsomericSMILES (includes stereochemistry) when available,
        # fall back to ConnectivitySMILES (PubChem's newer plain-SMILES field)
        smiles = props.get("IsomericSMILES") or props.get("ConnectivitySMILES")
        if not smiles:
            log.warning(f"  [pubchem] No SMILES field found for '{name}'")
            return None
        return {
            "name"   : name,
            "cid"    : props["CID"],
            "smiles" : smiles,
            "mw"     : props.get("MolecularWeight", "N/A"),
            "iupac"  : props.get("IUPACName", ""),
        }
    except (KeyError, IndexError, json.JSONDecodeError) as e:
        log.warning(f"  [pubchem] Could not parse response for '{name}': {e}")
        return None


def fetch_compound_library(
    names      : List[str],
    cache_path : Path = None,
) -> List[Dict]:
    """
    Fetch SMILES for a list of compound names from PubChem, with caching.

    Args:
        names      : list of compound names to fetch
        cache_path : optional JSON file to cache results
                     (avoids re-querying PubChem on every pipeline run)

    Returns:
        List of dicts: {"name", "cid", "smiles", "mw", "note"}
        Compounds that fail to fetch are skipped with a warning.
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  PUBCHEM COMPOUND LIBRARY FETCH ({len(names)} compounds)")
    log.info(f"{'═'*50}")

    # Load cache if it exists
    cache = {}
    if cache_path and cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
            log.info(f"  [pubchem] Loaded {len(cache)} cached compounds from {cache_path.name}")
        except Exception:
            cache = {}

    results = []
    fetched_count = 0

    for name in names:
        if name in cache:
            results.append(cache[name])
            log.info(f"  [pubchem] {name} — using cached SMILES")
            continue

        log.info(f"  [pubchem] Fetching: {name} ...")
        compound = fetch_compound_by_name(name)

        if compound is None:
            log.warning(f"  [pubchem] ✗ Skipping {name} — not found")
            continue

        entry = {
            "name"   : compound["name"],
            "smiles" : compound["smiles"],
            "note"   : f"PubChem CID {compound['cid']} | MW {compound['mw']}",
        }
        results.append(entry)
        cache[name] = entry
        fetched_count += 1

        log.info(f"  [pubchem] ✓ {name} — CID {compound['cid']}, "
                 f"MW {compound['mw']}")

        time.sleep(REQUEST_DELAY)   # Respect PubChem rate limits

    # Save updated cache
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=2))
        log.info(f"\n  [pubchem] Cache updated: {cache_path}")

    log.info(f"\n  [pubchem] Fetch complete: {len(results)}/{len(names)} compounds "
             f"({fetched_count} newly fetched, {len(results)-fetched_count} from cache)")

    return results


# =============================================================================
# Curated CNS / Neuroprotective Compound List
# =============================================================================
# 25 compounds with established literature relevance to neuroprotection,
# MAO-B inhibition, antioxidant activity, or CNS drug-likeness.
# Mix of approved drugs (positive controls) and natural polyphenols.

NEUROPROTECTIVE_LIBRARY = [
    # ── Approved MAO-B inhibitors (positive controls) ─────────────────────────
    "Selegiline",
    "Rasagiline",
    "Safinamide",

    # ── Approved AChE inhibitors (cross-reference) ────────────────────────────
    "Donepezil",
    "Galantamine",
    "Rivastigmine",

    # ── Natural polyphenols — neuroprotective / antioxidant ───────────────────
    "Quercetin",
    "Resveratrol",
    "Curcumin",
    "Kaempferol",
    "Myricetin",
    "Apigenin",
    "Luteolin",
    "Genistein",
    "Catechin",
    "Epicatechin",
    "Epigallocatechin gallate",
    "Fisetin",
    "Naringenin",
    "Hesperetin",
    "Chrysin",
    "Baicalein",
    "Rutin",

    # ── Neurotrophic / TrkB-related compounds ─────────────────────────────────
    "7,8-Dihydroxyflavone",
    "Hesperidin",

    # ── Other CNS-relevant natural compounds ──────────────────────────────────
    "Ferulic acid",
    "Caffeic acid",
]