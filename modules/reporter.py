# =============================================================================
# modules/reporter.py
# =============================================================================
# Saves docking results to CSV and prints a formatted summary table.
#
# Outputs:
#   data/results/docking_results.csv   — full ranked results table
#   data/results/hits_only.csv         — only PASS candidates
#   data/results/pipeline.log          — full pipeline log
# =============================================================================

import csv
import logging
from datetime import datetime
from pathlib import Path
from typing import List, Dict

log = logging.getLogger(__name__)


# Columns we want in the CSV (in order)
CSV_COLUMNS = [
    "rank",
    "ligand",
    "target",
    "affinity_kcal_mol",
    "rmsd_lb",
    "rmsd_ub",
    "decision",
    # ADMET columns (present only if admet module was run)
    "mw",
    "logp",
    "hbd",
    "hba",
    "tpsa",
    "rot_bonds",
    "logbb",
    "sa_score",
    "lipinski_violations",
    "pains_count",
    "admet_decision",
    "admet_failures",
    "final_decision",
]


def save_csv(results: List[Dict], out_path: Path, columns: List[str] = CSV_COLUMNS) -> None:
    """
    Save results to a CSV file.

    Args:
        results  : list of result dicts
        out_path : path to write the CSV
        columns  : which keys to include (and their order)
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    log.info(f"  [reporter] CSV saved: {out_path}")


def add_ranks(results: List[Dict]) -> List[Dict]:
    """Add a rank number (1 = best) to each result row."""
    for i, row in enumerate(results, start=1):
        row["rank"] = i
    return results


def print_summary(results: List[Dict], cutoff: float) -> None:
    """
    Print a formatted summary table to the terminal.

    Uses tabulate if available, falls back to plain text.
    """
    try:
        from tabulate import tabulate
        use_tabulate = True
    except ImportError:
        use_tabulate = False

    # Check if ADMET was run
    admet_ran = "admet_decision" in (results[0] if results else {})

    hits  = [r for r in results if r["decision"] == "PASS"]
    fails = [r for r in results if r["decision"] == "FAIL"]
    leads = [r for r in results if r.get("final_decision") == "✓ LEAD CANDIDATE"]

    print("\n" + "═" * 72)
    print(f"  NEURO DOCK — RESULTS SUMMARY")
    print(f"  Run at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Docking threshold: ≤ {cutoff} kcal/mol")
    if admet_ran:
        print(f"  ADMET: CNS filters applied (TPSA, logBB, PAINS, Lipinski)")
    print("═" * 72)

    if not results:
        print("  No results to display.")
        return

    # ── Lead candidates table (passed BOTH docking + ADMET) ─────────────────
    if admet_ran:
        print(f"\n  🏆 LEAD CANDIDATES — passed docking AND all CNS filters ({len(leads)})\n")
        if leads:
            rows = [
                [
                    r["rank"],
                    r["ligand"],
                    r["target"][:30] + ("…" if len(r["target"]) > 30 else ""),
                    f"{r['affinity_kcal_mol']:.1f}",
                    f"{r.get('logbb', 'N/A')}",
                    f"{r.get('tpsa', 'N/A')}",
                    f"{r.get('sa_score', 'N/A')}",
                    r.get("admet_decision", "N/A"),
                ]
                for r in leads
            ]
            headers = ["Rank", "Ligand", "Target", "ΔG", "logBB", "TPSA", "SA", "ADMET"]
            if use_tabulate:
                print(tabulate(rows, headers=headers, tablefmt="rounded_outline"))
            else:
                for row in rows:
                    print(f"  {row}")
        else:
            print("  No compounds passed all filters.")
            print("  Consider relaxing thresholds in config/settings.py or adding more ligands.")

    # ── Docking hits table ───────────────────────────────────────────────────
    print(f"\n  ✓ DOCKING HITS ({len(hits)} compounds ≤ {cutoff} kcal/mol)\n")
    if hits:
        if admet_ran:
            rows = [
                [
                    r["rank"],
                    r["ligand"],
                    r["target"][:28] + ("…" if len(r["target"]) > 28 else ""),
                    f"{r['affinity_kcal_mol']:.1f}",
                    f"{r.get('logbb', 'N/A')}",
                    f"{r.get('tpsa', 'N/A')} Å²",
                    r.get("admet_decision", "N/A"),
                    r.get("final_decision", "N/A"),
                ]
                for r in hits
            ]
            headers = ["Rank", "Ligand", "Target", "ΔG", "logBB", "TPSA", "ADMET", "Final"]
        else:
            rows = [
                [
                    r["rank"],
                    r["ligand"],
                    r["target"][:35] + ("…" if len(r["target"]) > 35 else ""),
                    f"{r['affinity_kcal_mol']:.1f}",
                    f"{r['rmsd_lb']:.2f}",
                ]
                for r in hits
            ]
            headers = ["Rank", "Ligand", "Target", "ΔG (kcal/mol)", "RMSD l.b."]

        if use_tabulate:
            print(tabulate(rows, headers=headers, tablefmt="rounded_outline"))
        else:
            for row in rows:
                print(f"  {row}")

    # ── ADMET failure details ────────────────────────────────────────────────
    if admet_ran and hits:
        admet_fails = [r for r in hits if r.get("admet_decision") == "FAIL"]
        if admet_fails:
            print(f"\n  ✗ ADMET failure reasons for docking hits:\n")
            for r in admet_fails:
                print(f"    {r['ligand']} vs {r['target'][:30]}:")
                for reason in r.get("admet_failures", "").split("; "):
                    if reason and reason != "None":
                        print(f"      • {reason}")

    # ── Did not pass docking ─────────────────────────────────────────────────
    print(f"\n  ✗ DID NOT PASS DOCKING THRESHOLD ({len(fails)} compounds)\n")
    if fails:
        rows = [
            [r["rank"], r["ligand"], r["target"][:35], f"{r['affinity_kcal_mol']:.1f}"]
            for r in fails
        ]
        headers = ["Rank", "Ligand", "Target", "ΔG (kcal/mol)"]
        if use_tabulate:
            print(tabulate(rows, headers=headers, tablefmt="simple"))
        else:
            for row in rows:
                print(f"  {row[0]:<4} {row[1]:<20} {row[2]:<36} {row[3]}")

    print("\n" + "═" * 72)
    if admet_ran:
        print(f"  Lead candidates saved to: data/results/leads_only.csv")
    print(f"  Full results: data/results/docking_results.csv")
    print("═" * 72 + "\n")


def generate_report(results: List[Dict], results_dir: Path, cutoff: float) -> None:
    """
    Full reporting pipeline: rank → save CSVs → print summary.

    This is the main function called by run_pipeline.py.

    Args:
        results     : scored results list from scorer.score_results()
        results_dir : folder to write output files
        cutoff      : binding energy cutoff (for display)
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  REPORTING")
    log.info(f"{'═'*50}")

    if not results:
        log.warning("  [reporter] No results to report.")
        print("\n  ⚠ No docking results were generated. Check the pipeline log.")
        return

    # Add rank numbers
    results = add_ranks(results)

    # Save full results
    all_csv = results_dir / "docking_results.csv"
    save_csv(results, all_csv)

    # Save hits only
    hits = [r for r in results if r["decision"] == "PASS"]
    if hits:
        hits_csv = results_dir / "hits_only.csv"
        save_csv(hits, hits_csv)
        log.info(f"  [reporter] Hits-only CSV saved: {hits_csv}")

    # Save lead candidates (passed docking + ADMET)
    leads = [r for r in results if r.get("final_decision") == "✓ LEAD CANDIDATE"]
    if leads:
        leads_csv = results_dir / "leads_only.csv"
        save_csv(leads, leads_csv)
        log.info(f"  [reporter] Lead candidates CSV saved: {leads_csv}")

    # Print to terminal
    print_summary(results, cutoff)