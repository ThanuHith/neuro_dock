# =============================================================================
# modules/scorer.py
# =============================================================================
# Filters and ranks docking results using pre-specified thresholds.
#
# From your paper (Section 3.7 — Decision Criteria):
#   "Binding affinity cutoff of ≤ −7.0 kcal/mol"
#
# This module:
#   1. Takes all raw docking results (every mode for every ligand-target pair)
#   2. Keeps only mode 1 (the best pose) for each pair
#   3. Applies the binding energy threshold
#   4. Ranks survivors by affinity (best first)
#   5. Adds a "PASS" / "FAIL" label
# =============================================================================

import logging
from typing import List, Dict

log = logging.getLogger(__name__)


def filter_best_poses(all_results: List[Dict]) -> List[Dict]:
    """
    From all docking modes, keep only mode 1 (the best pose) per ligand-target pair.

    AutoDock Vina returns up to num_modes poses ranked by energy.
    Mode 1 is always the lowest (best) energy — the one we care about.

    Args:
        all_results : flat list of all docking result dicts from docking.py

    Returns:
        List of dicts — one per unique ligand-target pair (mode 1 only)
    """
    best = {}

    for row in all_results:
        if row["mode"] == 1:
            key = (row["ligand"], row["target"])
            best[key] = row

    best_list = list(best.values())
    log.info(f"  [scorer] {len(best_list)} unique ligand-target pairs (mode 1 selected)")
    return best_list


def apply_threshold(best_poses: List[Dict], cutoff: float) -> List[Dict]:
    """
    Apply the binding energy cutoff and label each result PASS or FAIL.

    Args:
        best_poses : list of best-pose dicts from filter_best_poses()
        cutoff     : binding energy threshold (e.g. -7.0 kcal/mol)

    Returns:
        All poses with an added "decision" key ("PASS" or "FAIL")
        Sorted from best (most negative) to worst affinity.
    """
    for row in best_poses:
        # More negative = better binding
        row["decision"] = "PASS" if row["affinity_kcal_mol"] <= cutoff else "FAIL"

    # Sort: best affinity first (most negative number first)
    scored = sorted(best_poses, key=lambda x: x["affinity_kcal_mol"])

    n_pass = sum(1 for r in scored if r["decision"] == "PASS")
    n_fail = len(scored) - n_pass

    log.info(f"  [scorer] Threshold: ≤ {cutoff} kcal/mol")
    log.info(f"  [scorer] Results: {n_pass} PASS, {n_fail} FAIL")

    return scored


def score_results(all_results: List[Dict], cutoff: float) -> List[Dict]:
    """
    Full scoring pipeline: filter best poses → apply threshold → rank.

    This is the main function called by run_pipeline.py.

    Args:
        all_results : raw results list from docking.dock_all()
        cutoff      : binding energy threshold from config/settings.py

    Returns:
        Ranked, labelled list of results ready for reporter.py
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  SCORING & FILTERING")
    log.info(f"{'═'*50}")

    if not all_results:
        log.warning("  [scorer] No docking results to score.")
        return []

    best  = filter_best_poses(all_results)
    final = apply_threshold(best, cutoff)
    return final