# =============================================================================
# modules/figures.py
# =============================================================================
# Generates publication-quality figures from NeuroDock pipeline results.
#
# Reads directly from data/results/docking_results.csv (already produced
# by reporter.py) so this module can be run standalone after a pipeline run,
# without needing to redo docking.
#
# Output: data/results/figures/
#   fig1_binding_energy_bar.png       — ranked ΔG bar chart, all compounds
#   fig2_admet_radar.png              — radar plot, ADMET properties of leads
#   fig3_logp_vs_affinity_scatter.png — SAR-style scatter plot
#   fig4_admet_heatmap.png            — pass/fail heatmap across filters
#   fig5_pipeline_flowchart.png       — pipeline overview diagram
#
# Design choices (journal-figure conventions, not web-UI conventions):
#   - White background, black axes, serif-adjacent sans font (Arial/Helvetica)
#   - Colorblind-safe palette (Okabe-Ito inspired)
#   - 300 DPI export — standard minimum for journal submission
#   - No unnecessary gridlines, no 3D effects, no gradients
# =============================================================================

import logging
import csv
from pathlib import Path
from typing import List, Dict

log = logging.getLogger(__name__)


# ── Modern refined palette (colorblind-safe, higher saturation/contrast) ────
COLORS = {
    "pass"      : "#2563EB",   # vivid indigo-blue
    "fail"      : "#EF4444",   # clean red
    "lead"      : "#059669",   # emerald green
    "neutral"   : "#94A3B8",   # cool slate grey
    "accent"    : "#F59E0B",   # amber
    "highlight" : "#DB2777",   # magenta-pink
    "bg_panel"  : "#F8FAFC",   # very light cool grey for panel backgrounds
    "grid"      : "#E2E8F0",   # soft gridline grey
    "text_dark" : "#0F172A",   # near-black slate for text
}

# Radar plot multi-series palette (distinct, modern, still colorblind-considerate)
RADAR_PALETTE = ["#059669", "#2563EB", "#F59E0B", "#DB2777", "#7C3AED", "#0891B2"]


def _setup_style():
    """Configure matplotlib for a modern, polished journal-acceptable look."""
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend, safe for scripts
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family"       : "sans-serif",
        "font.sans-serif"   : ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size"         : 10.5,
        "text.color"        : COLORS["text_dark"],
        "axes.linewidth"    : 1.0,
        "axes.edgecolor"    : "#334155",
        "axes.labelweight"  : "bold",
        "axes.labelcolor"   : COLORS["text_dark"],
        "axes.titleweight"  : "bold",
        "axes.titlesize"    : 13.5,
        "axes.titlecolor"   : COLORS["text_dark"],
        "xtick.color"       : "#334155",
        "ytick.color"       : "#334155",
        "xtick.labelsize"   : 9.5,
        "ytick.labelsize"   : 9.5,
        "figure.facecolor"  : "white",
        "savefig.facecolor" : "white",
        "savefig.dpi"       : 300,
        "legend.frameon"    : False,
    })
    return plt


def load_results(csv_path: Path) -> List[Dict]:
    """Load docking_results.csv into a list of dicts with correct types."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Results CSV not found: {csv_path}")

    rows = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert numeric fields
            for key in ["affinity_kcal_mol", "rmsd_lb", "rmsd_ub", "mw", "logp",
                        "tpsa", "logbb", "sa_score"]:
                if row.get(key) not in (None, "", "N/A"):
                    try:
                        row[key] = float(row[key])
                    except ValueError:
                        pass
            for key in ["hbd", "hba", "rot_bonds", "lipinski_violations",
                        "pains_count", "rank"]:
                if row.get(key) not in (None, "", "N/A"):
                    try:
                        row[key] = int(float(row[key]))
                    except ValueError:
                        pass
            rows.append(row)

    # Defensive filter: drop any row where the core numeric field
    # (binding affinity) didn't parse to a float — guards against
    # malformed CSV rows (e.g. unquoted commas in compound names)
    # silently corrupting downstream sorting/plotting.
    clean_rows = []
    dropped = 0
    for row in rows:
        if isinstance(row.get("affinity_kcal_mol"), (int, float)):
            clean_rows.append(row)
        else:
            dropped += 1
    if dropped:
        log.warning(f"  [figures] Dropped {dropped} malformed row(s) from CSV "
                   f"(non-numeric affinity_kcal_mol — check for unquoted "
                   f"commas in compound names)")

    return clean_rows


# =============================================================================
# Figure 1 — Binding Energy Bar Chart (all compounds, ranked)
# =============================================================================

def fig1_binding_energy_bar(results: List[Dict], out_dir: Path, cutoff: float) -> Path:
    plt = _setup_style()
    import matplotlib.pyplot as plt

    sorted_results = sorted(results, key=lambda x: x["affinity_kcal_mol"])

    names      = [r["ligand"] for r in sorted_results]
    energies   = [r["affinity_kcal_mol"] for r in sorted_results]
    is_lead    = [r.get("final_decision") == "✓ LEAD CANDIDATE" for r in sorted_results]
    is_hit     = [r["decision"] == "PASS" for r in sorted_results]

    bar_colors = []
    for lead, hit in zip(is_lead, is_hit):
        if lead:
            bar_colors.append(COLORS["lead"])
        elif hit:
            bar_colors.append(COLORS["pass"])
        else:
            bar_colors.append(COLORS["neutral"])

    fig_height = max(5, len(names) * 0.34)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    fig.patch.set_facecolor("white")
    ax.set_facecolor(COLORS["bg_panel"])

    y_pos = range(len(names))
    ax.barh(y_pos, energies, color=bar_colors, edgecolor="white",
           linewidth=0.6, height=0.68, zorder=3)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8.8, color=COLORS["text_dark"])
    ax.invert_yaxis()

    ax.axvline(x=cutoff, color=COLORS["text_dark"], linestyle=(0, (5, 3)),
              linewidth=1.3, zorder=4, alpha=0.75)
    ax.text(cutoff, -1.3, f" Threshold: {cutoff} kcal/mol", fontsize=8.5,
           color=COLORS["text_dark"], fontweight="bold", ha="left", va="bottom")

    ax.set_xlabel("Binding Affinity  ΔG  (kcal/mol)", fontsize=11.5, labelpad=10)
    ax.set_title("Molecular Docking Results — MAO-B Virtual Screen",
                pad=16, fontsize=14)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS["lead"], edgecolor="white",
              label="Lead candidate (docking + ADMET pass)"),
        Patch(facecolor=COLORS["pass"], edgecolor="white",
              label="Docking hit (ADMET fail)"),
        Patch(facecolor=COLORS["neutral"], edgecolor="white",
              label="Below threshold"),
    ]
    leg = ax.legend(handles=legend_elements, loc="upper center",
                    bbox_to_anchor=(0.5, -0.075), ncol=1, fontsize=8.5,
                    frameon=True, edgecolor=COLORS["grid"], facecolor="white")
    leg.get_frame().set_linewidth(0.8)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#334155")

    ax.grid(axis="x", linestyle="-", linewidth=0.6, color=COLORS["grid"], alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=3, color="#334155")

    plt.tight_layout()
    out_path = out_dir / "fig1_binding_energy_bar.png"
    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info(f"  [figures] Saved: {out_path.name}")
    return out_path


# =============================================================================
# Figure 2 — ADMET Radar Plot (lead candidates only)
# =============================================================================

def fig2_admet_radar(results: List[Dict], out_dir: Path) -> Path:
    import matplotlib.pyplot as plt
    import numpy as np

    leads = [r for r in results if r.get("final_decision") == "✓ LEAD CANDIDATE"]
    if not leads:
        log.warning("  [figures] No lead candidates — skipping radar plot")
        return None

    def normalise(value, lo, hi, invert=False):
        try:
            v = float(value)
        except (TypeError, ValueError):
            return 0.5
        v = max(lo, min(hi, v))
        score = (v - lo) / (hi - lo)
        return 1 - score if invert else score

    categories = ["BBB Permeability\n(logBB)", "Low TPSA\n(CNS access)",
                  "Synthetic\nAccessibility", "Lipinski\nCompliance",
                  "Low\nMolecular\nWeight", "PAINS-free"]
    n = len(categories)
    angles = [i / n * 2 * np.pi for i in range(n)]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(8.2, 7), subplot_kw=dict(polar=True))
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.16, right=0.74, top=0.88, bottom=0.08)

    # Soft background rings instead of plain grey gridlines
    ax.set_facecolor("#FAFBFC")

    for i, r in enumerate(leads):
        values = [
            normalise(r.get("logbb", 0), -2, 1),
            normalise(r.get("tpsa", 90), 0, 150, invert=True),
            normalise(r.get("sa_score", 5), 1, 8, invert=True),
            normalise(r.get("lipinski_violations", 0), 0, 3, invert=True),
            normalise(r.get("mw", 300), 100, 600, invert=True),
            1.0 if int(r.get("pains_count", 0) or 0) == 0 else 0.0,
        ]
        values += values[:1]

        color = RADAR_PALETTE[i % len(RADAR_PALETTE)]
        ax.plot(angles, values, linewidth=2.4, label=r["ligand"], color=color,
               solid_capstyle="round", zorder=3)
        ax.scatter(angles[:-1], values[:-1], s=28, color=color,
                  edgecolor="white", linewidth=1.0, zorder=4)
        ax.fill(angles, values, alpha=0.06, color=color, zorder=1)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=9.5, color=COLORS["text_dark"],
                       fontweight="medium")
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=8, color="#94A3B8")
    ax.spines["polar"].set_color("#CBD5E1")
    ax.spines["polar"].set_linewidth(1.0)
    ax.grid(color="#E2E8F0", linewidth=0.9)

    ax.set_title("CNS Drug-Likeness Profile — Lead Candidates",
                pad=28, fontsize=14.5)

    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.18, 1.08), fontsize=9.5,
                    frameon=True, edgecolor=COLORS["grid"], facecolor="white",
                    title="Lead Candidates", title_fontsize=10)
    leg.get_frame().set_linewidth(0.8)

    out_path = out_dir / "fig2_admet_radar.png"
    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info(f"  [figures] Saved: {out_path.name}")
    return out_path


# =============================================================================
# Figure 3 — logP vs Binding Affinity Scatter (SAR-style)
# =============================================================================

def fig3_logp_vs_affinity_scatter(results: List[Dict], out_dir: Path, cutoff: float) -> Path:
    import matplotlib.pyplot as plt

    valid = [r for r in results if r.get("logp") not in (None, "", "N/A")
             and r.get("affinity_kcal_mol") is not None]

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor(COLORS["bg_panel"])

    ax.axvspan(-0.5, 5.0, color=COLORS["accent"], alpha=0.07, zorder=0,
              label="CNS-favourable logP range")

    for r in valid:
        is_lead = r.get("final_decision") == "✓ LEAD CANDIDATE"
        is_hit  = r["decision"] == "PASS"

        if is_lead:
            color, marker, size, z = COLORS["lead"], "D", 130, 4
        elif is_hit:
            color, marker, size, z = COLORS["pass"], "o", 75, 3
        else:
            color, marker, size, z = COLORS["neutral"], "o", 50, 2

        ax.scatter(r["logp"], r["affinity_kcal_mol"], c=color, marker=marker,
                  s=size, edgecolor="white", linewidth=1.0, zorder=z, alpha=0.92)

        if is_lead:
            ax.annotate(r["ligand"], (r["logp"], r["affinity_kcal_mol"]),
                       fontsize=8.5, xytext=(7, 5), textcoords="offset points",
                       fontweight="bold", color=COLORS["text_dark"],
                       bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                                edgecolor=COLORS["lead"], linewidth=0.8, alpha=0.92))

    ax.axhline(y=cutoff, color=COLORS["text_dark"], linestyle=(0, (5, 3)),
              linewidth=1.3, alpha=0.75, zorder=1,
              label=f"Docking threshold ({cutoff} kcal/mol)")

    ax.set_xlabel("Calculated logP (lipophilicity)", fontsize=11.5, labelpad=8)
    ax.set_ylabel("Binding Affinity  ΔG  (kcal/mol)", fontsize=11.5, labelpad=8)
    ax.set_title("Lipophilicity vs. Binding Affinity — MAO-B Screen", pad=16, fontsize=14)
    ax.invert_yaxis()

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COLORS["lead"],
              markeredgecolor="white", markersize=11, label="Lead candidate"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS["pass"],
              markeredgecolor="white", markersize=10, label="Docking hit"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS["neutral"],
              markeredgecolor="white", markersize=9, label="Below threshold"),
    ]
    leg = ax.legend(handles=legend_elements, loc="lower left", fontsize=9,
                    frameon=True, edgecolor=COLORS["grid"], facecolor="white")
    leg.get_frame().set_linewidth(0.8)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#334155")
    ax.spines["bottom"].set_color("#334155")
    ax.grid(linestyle="-", linewidth=0.6, color=COLORS["grid"], alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, color="#334155")

    plt.tight_layout()
    out_path = out_dir / "fig3_logp_vs_affinity_scatter.png"
    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info(f"  [figures] Saved: {out_path.name}")
    return out_path


# =============================================================================
# Figure 4 — ADMET Pass/Fail Heatmap
# =============================================================================

def fig4_admet_heatmap(results: List[Dict], out_dir: Path) -> Path:
    import matplotlib.pyplot as plt
    import numpy as np

    seen = {}
    for r in results:
        if r["ligand"] not in seen:
            seen[r["ligand"]] = r
    unique = list(seen.values())
    unique = sorted(unique, key=lambda x: x["affinity_kcal_mol"])

    def safe_num(value, default):
        if value is None or value == "":
            return default
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    criteria = [
        ("MW \u2264450",        lambda r: safe_num(r.get("mw"), 9999) <= 450),
        ("logP -0.5\u20135.0",  lambda r: -0.5 <= safe_num(r.get("logp"), -99) <= 5.0),
        ("HBD \u22643",         lambda r: safe_num(r.get("hbd"), 99) <= 3),
        ("HBA \u22647",         lambda r: safe_num(r.get("hba"), 99) <= 7),
        ("TPSA <90",           lambda r: safe_num(r.get("tpsa"), 999) < 90),
        ("logBB >-1.0",        lambda r: safe_num(r.get("logbb"), -99) > -1.0),
        ("Lipinski \u22641",    lambda r: safe_num(r.get("lipinski_violations"), 9) <= 1),
        ("PAINS-free",         lambda r: safe_num(r.get("pains_count"), 9) == 0),
        ("Docking \u2264-7.0",  lambda r: r["decision"] == "PASS"),
    ]

    matrix = []
    for r in unique:
        row = []
        for _, fn in criteria:
            try:
                row.append(1 if fn(r) else 0)
            except (TypeError, ValueError):
                row.append(0)
        matrix.append(row)

    matrix = np.array(matrix)
    names = [r["ligand"] for r in unique]
    labels = [c[0] for c in criteria]

    fig_height = max(5.5, len(names) * 0.32)
    fig, ax = plt.subplots(figsize=(8, fig_height))
    fig.patch.set_facecolor("white")

    cmap_colors = [COLORS["fail"], COLORS["lead"]]
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(cmap_colors)

    ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=0, vmax=1, alpha=0.88)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8.8, rotation=38, ha="right",
                       color=COLORS["text_dark"])
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8.8, color=COLORS["text_dark"])

    ax.set_xticks([x - 0.5 for x in range(1, len(labels))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(names))], minor=True)
    ax.grid(which="minor", color="white", linewidth=2.2)
    ax.tick_params(which="minor", size=0)
    ax.tick_params(which="major", length=0)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title("ADMET & Docking Criteria — Pass/Fail Matrix", pad=16, fontsize=14)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS["lead"], edgecolor="none", label="Pass", alpha=0.88),
        Patch(facecolor=COLORS["fail"], edgecolor="none", label="Fail", alpha=0.88),
    ]
    leg = ax.legend(handles=legend_elements, loc="upper center",
                    bbox_to_anchor=(0.5, -0.13 - 0.008 * len(names)), ncol=2,
                    fontsize=9.5, frameon=True, edgecolor=COLORS["grid"],
                    facecolor="white")
    leg.get_frame().set_linewidth(0.8)

    plt.tight_layout()
    out_path = out_dir / "fig4_admet_heatmap.png"
    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info(f"  [figures] Saved: {out_path.name}")
    return out_path


# =============================================================================
# Figure 5 — Pipeline Flowchart
# =============================================================================

def fig5_pipeline_flowchart(out_dir: Path) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    fig, ax = plt.subplots(figsize=(4.6, 10.5))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 24)
    ax.axis("off")

    steps = [
        ("Protein Structures",         "(RCSB PDB)",                    COLORS["neutral"]),
        ("Receptor Preparation",       "(Open Babel \u2192 PDBQT)",     COLORS["pass"]),
        ("Compound Library",           "(PubChem, n=29)",               COLORS["neutral"]),
        ("Ligand Preparation",         "(RDKit ETKDG / MMFF94)",        COLORS["pass"]),
        ("Molecular Docking",          "(AutoDock Vina 1.2.5)",         COLORS["accent"]),
        ("Binding Energy Filter",      "(\u2264 \u22127.0 kcal/mol)",   COLORS["accent"]),
        ("ADMET & CNS Filtering",      "(RDKit descriptors)",           COLORS["highlight"]),
        ("Interaction Analysis",       "(H-bonds, hydrophobic, \u03c0-stacking)", COLORS["highlight"]),
        ("Lead Candidates",            "(n = 6)",                       COLORS["lead"]),
    ]

    n = len(steps)
    box_h = 2.05
    gap = 0.55
    total_h = n * box_h + (n - 1) * gap
    y_start = (24 - total_h) / 2 + total_h - box_h

    positions = []
    y = y_start
    for title, subtitle, color in steps:
        positions.append((y, title, subtitle, color))
        y -= (box_h + gap)

    for i, (y, title, subtitle, color) in enumerate(positions):
        # Soft shadow behind each box for depth
        shadow = FancyBboxPatch((1.07, y - 0.07), 8, box_h,
                                boxstyle="round,pad=0.05,rounding_size=0.22",
                                linewidth=0, facecolor="#00000018", zorder=1)
        ax.add_patch(shadow)

        box = FancyBboxPatch((1, y), 8, box_h,
                             boxstyle="round,pad=0.05,rounding_size=0.22",
                             linewidth=0, facecolor=color, alpha=0.95, zorder=2)
        ax.add_patch(box)

        # Step number badge
        ax.add_patch(plt.Circle((1.65, y + box_h - 0.42), 0.32,
                                facecolor="white", edgecolor="none", zorder=3, alpha=0.92))
        ax.text(1.65, y + box_h - 0.42, str(i + 1), ha="center", va="center",
               fontsize=10, fontweight="bold", color=color, zorder=4)

        ax.text(5.3, y + box_h / 2 + 0.18, title, ha="center", va="center",
               fontsize=10.5, fontweight="bold", color="white", zorder=4)
        ax.text(5.3, y + box_h / 2 - 0.35, subtitle, ha="center", va="center",
               fontsize=8, color="white", alpha=0.92, zorder=4, style="italic")

        if i < len(positions) - 1:
            next_y = positions[i + 1][0]
            ax.annotate("", xy=(5, next_y + box_h + 0.06), xytext=(5, y - 0.06),
                       arrowprops=dict(arrowstyle="-|>", linewidth=1.6,
                                       color="#475569", shrinkA=0, shrinkB=0),
                       zorder=1)

    ax.set_title("NeuroDock Pipeline Overview", fontsize=14.5,
                fontweight="bold", pad=18, color=COLORS["text_dark"])

    plt.tight_layout()
    out_path = out_dir / "fig5_pipeline_flowchart.png"
    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info(f"  [figures] Saved: {out_path.name}")
    return out_path


# =============================================================================
# Main entry point
# =============================================================================

def generate_all_figures(results_dir: Path, cutoff: float) -> List[Path]:
    """
    Generate all publication figures from the pipeline's CSV output.

    Args:
        results_dir : data/results directory (contains docking_results.csv)
        cutoff      : binding energy threshold used in the screen

    Returns:
        List of paths to generated figure files
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  FIGURE GENERATION")
    log.info(f"{'═'*50}")

    csv_path = results_dir / "docking_results.csv"
    fig_dir = results_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    try:
        results = load_results(csv_path)
    except FileNotFoundError as e:
        log.error(f"  [figures] {e}")
        return []

    if not results:
        log.warning("  [figures] No results found in CSV — skipping figures")
        return []

    generated = []

    try:
        generated.append(fig1_binding_energy_bar(results, fig_dir, cutoff))
    except Exception as e:
        log.error(f"  [figures] Figure 1 failed: {e}")

    try:
        p = fig2_admet_radar(results, fig_dir)
        if p:
            generated.append(p)
    except Exception as e:
        log.error(f"  [figures] Figure 2 failed: {e}")

    try:
        generated.append(fig3_logp_vs_affinity_scatter(results, fig_dir, cutoff))
    except Exception as e:
        log.error(f"  [figures] Figure 3 failed: {e}")

    try:
        generated.append(fig4_admet_heatmap(results, fig_dir))
    except Exception as e:
        log.error(f"  [figures] Figure 4 failed: {e}")

    try:
        generated.append(fig5_pipeline_flowchart(fig_dir))
    except Exception as e:
        log.error(f"  [figures] Figure 5 failed: {e}")

    log.info(f"\n  [figures] Generated {len(generated)} figures in {fig_dir}")
    return generated