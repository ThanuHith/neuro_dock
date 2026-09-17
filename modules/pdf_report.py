# =============================================================================
# modules/pdf_report.py
# =============================================================================
# Generates a professional PDF report of the full NeuroDock pipeline results.
#
# Output: data/results/NeuroDock_Report.pdf
#
# Report structure:
#   Page 1 — Title page
#   Page 2 — Executive summary + lead candidate highlight
#   Page 3 — Full docking results table
#   Page 4 — ADMET analysis table
#   Page 5 — Methods & thresholds used
# =============================================================================

import logging
from datetime import datetime
from pathlib import Path
from typing import List, Dict

log = logging.getLogger(__name__)


# ── Colour palette ────────────────────────────────────────────────────────────
DARK_BLUE   = (0.08, 0.18, 0.35)    # Title / headers
MID_BLUE    = (0.18, 0.38, 0.62)    # Subheadings
ACCENT      = (0.04, 0.62, 0.55)    # Pass / highlight green-teal
RED         = (0.75, 0.15, 0.15)    # Fail
LIGHT_GREY  = (0.94, 0.95, 0.97)    # Table row shading
MID_GREY    = (0.55, 0.55, 0.60)    # Secondary text
WHITE       = (1.0,  1.0,  1.0)
BLACK       = (0.0,  0.0,  0.0)


def _rgb(triplet):
    """Convert 0-1 float triplet to reportlab Color."""
    from reportlab.lib.colors import Color
    return Color(*triplet)


def build_pdf(
    results      : List[Dict],
    results_dir  : Path,
    cutoff       : float,
    targets      : List[Dict],
    ligands      : List[Dict],
    interactions : list = None,
) -> Path:
    """
    Build the full NeuroDock PDF report.

    Args:
        results     : scored + ADMET-analysed results from the pipeline
        results_dir : folder to save the PDF
        cutoff      : binding energy threshold used
        targets     : target list from settings.py
        ligands     : ligand list from settings.py

    Returns:
        Path to the generated PDF file
    """
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import cm
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
            PageBreak, HRFlowable
        )
        from reportlab.lib.colors import Color, HexColor
    except ImportError:
        raise ImportError(
            "reportlab not installed. Run: pip install reportlab"
        )

    results_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = results_dir / "NeuroDock_Report.pdf"

    # ── Document setup ────────────────────────────────────────────────────────
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize    = A4,
        leftMargin  = 2.2 * cm,
        rightMargin = 2.2 * cm,
        topMargin   = 2.0 * cm,
        bottomMargin= 2.0 * cm,
        title       = "NeuroDock — Drug Discovery Pipeline Report",
        author      = "NeuroDock Pipeline",
    )

    W = A4[0] - 4.4 * cm   # usable page width

    # ── Styles ────────────────────────────────────────────────────────────────
    styles = getSampleStyleSheet()

    def style(name, **kw):
        return ParagraphStyle(name, **kw)

    S = {
        "title"  : style("T",  fontSize=26, leading=32, textColor=_rgb(DARK_BLUE),
                          alignment=TA_CENTER, fontName="Helvetica-Bold", spaceAfter=6),
        "subtitle": style("ST", fontSize=13, leading=18, textColor=_rgb(MID_BLUE),
                          alignment=TA_CENTER, fontName="Helvetica", spaceAfter=4),
        "date"   : style("D",  fontSize=10, textColor=_rgb(MID_GREY),
                          alignment=TA_CENTER, fontName="Helvetica"),
        "h1"     : style("H1", fontSize=14, leading=18, textColor=_rgb(DARK_BLUE),
                          fontName="Helvetica-Bold", spaceBefore=14, spaceAfter=6),
        "h2"     : style("H2", fontSize=11, leading=14, textColor=_rgb(MID_BLUE),
                          fontName="Helvetica-Bold", spaceBefore=10, spaceAfter=4),
        "body"   : style("B",  fontSize=9,  leading=13, textColor=_rgb(BLACK),
                          fontName="Helvetica", spaceAfter=4),
        "small"  : style("SM", fontSize=8,  leading=11, textColor=_rgb(MID_GREY),
                          fontName="Helvetica"),
        "lead"   : style("LD", fontSize=11, leading=15, textColor=_rgb(ACCENT),
                          fontName="Helvetica-Bold", spaceAfter=4),
        "th"     : style("TH", fontSize=8,  leading=10, textColor=_rgb(WHITE),
                          fontName="Helvetica-Bold", alignment=TA_CENTER),
        "td"     : style("TD", fontSize=8,  leading=10, textColor=_rgb(BLACK),
                          fontName="Helvetica", alignment=TA_CENTER),
        "td_l"   : style("TDL",fontSize=8,  leading=10, textColor=_rgb(BLACK),
                          fontName="Helvetica", alignment=TA_LEFT),
    }

    story = []

    def hr(color=DARK_BLUE, thickness=1):
        return HRFlowable(width="100%", thickness=thickness,
                          color=_rgb(color), spaceAfter=8, spaceBefore=4)

    def sp(h=0.3):
        return Spacer(1, h * cm)

    # =========================================================================
    # PAGE 1 — TITLE PAGE
    # =========================================================================
    story += [
        sp(3),
        Paragraph("NeuroDock", S["title"]),
        Paragraph("AI-Driven Molecular Docking Pipeline", S["subtitle"]),
        Paragraph("MDD &amp; Alzheimer's Disease Drug Discovery", S["subtitle"]),
        sp(0.5),
        hr(ACCENT, 2),
        sp(0.5),
        Paragraph(f"Report generated: {datetime.now().strftime('%d %B %Y, %H:%M')}", S["date"]),
        sp(0.3),
        Paragraph(f"Pipeline: AutoDock Vina v1.2.5 &nbsp;|&nbsp; RDKit 2026 &nbsp;|&nbsp; Open Babel 3.1", S["date"]),
        sp(3),
    ]

    # Summary box
    hits  = [r for r in results if r.get("decision") == "PASS"]
    leads = [r for r in results if r.get("final_decision") == "✓ LEAD CANDIDATE"]

    summary_data = [
        ["Targets screened", str(len(targets))],
        ["Compounds screened", str(len(ligands))],
        ["Docking runs", str(len(results))],
        ["Docking hits (threshold)", str(len(hits))],
        ["Lead candidates (docking + ADMET)", str(len(leads))],
        ["Binding energy threshold", f"≤ {cutoff} kcal/mol"],
    ]

    summary_table = Table(summary_data, colWidths=[10*cm, 5*cm])
    summary_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), _rgb(DARK_BLUE)),
        ("ROWBACKGROUNDS", (0,0), (-1,-1), [_rgb(LIGHT_GREY), _rgb(WHITE)]),
        ("TEXTCOLOR", (0,0), (-1,-1), _rgb(BLACK)),
        ("FONTNAME", (0,0), (0,-1), "Helvetica-Bold"),
        ("FONTNAME", (1,0), (1,-1), "Helvetica"),
        ("FONTSIZE", (0,0), (-1,-1), 10),
        ("ALIGN", (1,0), (1,-1), "CENTER"),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("ROWHEIGHT", (0,0), (-1,-1), 22),
        ("BOX", (0,0), (-1,-1), 0.5, _rgb(MID_BLUE)),
        ("INNERGRID", (0,0), (-1,-1), 0.3, _rgb(MID_GREY)),
        ("LEFTPADDING", (0,0), (-1,-1), 10),
        ("RIGHTPADDING", (0,0), (-1,-1), 10),
    ]))
    story += [summary_table, sp(1), PageBreak()]

    # =========================================================================
    # PAGE 2 — LEAD CANDIDATES
    # =========================================================================
    story += [
        Paragraph("Lead Candidates", S["h1"]),
        hr(),
        Paragraph(
            "The following compounds passed both the molecular docking energy threshold "
            f"(≤ {cutoff} kcal/mol) and all CNS ADMET filters including BBB permeability, "
            "TPSA, Lipinski Rule of Five, and PAINS screening.",
            S["body"]
        ),
        sp(0.4),
    ]

    if leads:
        for r in leads:
            story += [
                Paragraph(f"🏆  {r['ligand']}  vs  {r['target']}", S["lead"]),
                sp(0.2),
            ]
            props_data = [
                ["Property", "Value", "Threshold", "Status"],
                ["Binding Affinity (ΔG)", f"{r['affinity_kcal_mol']:.2f} kcal/mol",
                 f"≤ {cutoff}", "PASS"],
                ["logBB (BBB permeability)", str(r.get("logbb", "N/A")), "> -1.0", "PASS"],
                ["TPSA", f"{r.get('tpsa', 'N/A')} Å²", "< 90 Å²", "PASS"],
                ["MW", f"{r.get('mw', 'N/A')} Da", "≤ 450 Da", "PASS"],
                ["logP", str(r.get("logp", "N/A")), "-0.5 to 5.0", "PASS"],
                ["H-bond donors", str(r.get("hbd", "N/A")), "≤ 3", "PASS"],
                ["H-bond acceptors", str(r.get("hba", "N/A")), "≤ 7", "PASS"],
                ["Rotatable bonds", str(r.get("rot_bonds", "N/A")), "≤ 8", "PASS"],
                ["SA Score", str(r.get("sa_score", "N/A")), "≤ 6.0", "PASS"],
                ["Lipinski violations", str(r.get("lipinski_violations", "N/A")), "≤ 1", "PASS"],
                ["PAINS alerts", str(r.get("pains_count", "N/A")), "0", "PASS"],
            ]
            t = Table(props_data, colWidths=[6*cm, 3.5*cm, 3.5*cm, 2.5*cm])
            t.setStyle(TableStyle([
                ("BACKGROUND", (0,0), (-1,0), _rgb(DARK_BLUE)),
                ("TEXTCOLOR", (0,0), (-1,0), _rgb(WHITE)),
                ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                ("FONTSIZE", (0,0), (-1,-1), 8),
                ("ROWBACKGROUNDS", (0,1), (-1,-1), [_rgb(LIGHT_GREY), _rgb(WHITE)]),
                ("TEXTCOLOR", (-1,1), (-1,-1), _rgb(ACCENT)),
                ("FONTNAME", (-1,1), (-1,-1), "Helvetica-Bold"),
                ("ALIGN", (1,0), (-1,-1), "CENTER"),
                ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
                ("ROWHEIGHT", (0,0), (-1,-1), 16),
                ("BOX", (0,0), (-1,-1), 0.5, _rgb(MID_BLUE)),
                ("INNERGRID", (0,0), (-1,-1), 0.3, _rgb(MID_GREY)),
                ("LEFTPADDING", (0,0), (-1,-1), 6),
            ]))
            story += [t, sp(0.6)]
    else:
        story.append(Paragraph(
            "No compounds passed all filters in this screening run. "
            "Consider expanding the compound library or adjusting thresholds in config/settings.py.",
            S["body"]
        ))

    story.append(PageBreak())

    # =========================================================================
    # PAGE 3 — FULL DOCKING RESULTS TABLE
    # =========================================================================
    story += [
        Paragraph("Docking Results — All Compounds", S["h1"]),
        hr(),
        Paragraph(
            f"All {len(results)} ligand-target pairs ranked by binding affinity. "
            f"Compounds with ΔG ≤ {cutoff} kcal/mol are highlighted.",
            S["body"]
        ),
        sp(0.3),
    ]

    dock_header = ["Rank", "Ligand", "Target", "ΔG (kcal/mol)", "RMSD l.b.", "Docking"]
    dock_rows   = [dock_header]
    for r in results:
        decision_label = "PASS" if r["decision"] == "PASS" else "FAIL"
        dock_rows.append([
            str(r.get("rank", "")),
            r["ligand"],
            r["target"][:32] + ("…" if len(r["target"]) > 32 else ""),
            f"{r['affinity_kcal_mol']:.3f}",
            f"{r['rmsd_lb']:.2f}",
            decision_label,
        ])

    col_w = [1.2*cm, 3.5*cm, 6.5*cm, 2.8*cm, 2.2*cm, 1.8*cm]
    dock_table = Table(dock_rows, colWidths=col_w, repeatRows=1)

    ts = [
        ("BACKGROUND",  (0,0), (-1,0), _rgb(DARK_BLUE)),
        ("TEXTCOLOR",   (0,0), (-1,0), _rgb(WHITE)),
        ("FONTNAME",    (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE",    (0,0), (-1,-1), 8),
        ("ALIGN",       (0,0), (-1,-1), "CENTER"),
        ("ALIGN",       (1,0), (2,-1), "LEFT"),
        ("VALIGN",      (0,0), (-1,-1), "MIDDLE"),
        ("ROWHEIGHT",   (0,0), (-1,-1), 14),
        ("BOX",         (0,0), (-1,-1), 0.5, _rgb(MID_BLUE)),
        ("INNERGRID",   (0,0), (-1,-1), 0.25, _rgb(MID_GREY)),
        ("LEFTPADDING", (0,0), (-1,-1), 5),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [_rgb(LIGHT_GREY), _rgb(WHITE)]),
    ]
    # Highlight PASS rows in green tint
    for i, r in enumerate(results, start=1):
        if r["decision"] == "PASS":
            ts += [
                ("TEXTCOLOR", (-1, i), (-1, i), _rgb(ACCENT)),
                ("FONTNAME",  (-1, i), (-1, i), "Helvetica-Bold"),
            ]
        else:
            ts += [("TEXTCOLOR", (-1, i), (-1, i), _rgb(RED))]

    dock_table.setStyle(TableStyle(ts))
    story += [dock_table, PageBreak()]

    # =========================================================================
    # PAGE 4 — ADMET TABLE
    # =========================================================================
    story += [
        Paragraph("ADMET &amp; CNS Drug-Likeness Analysis", S["h1"]),
        hr(),
        Paragraph(
            "ADMET properties calculated per compound using RDKit. "
            "logBB estimated via Clark (1999) empirical model: "
            "logBB = -0.0148 × TPSA + 0.152 × logP + 0.139. "
            "PAINS alerts detected using the Baell &amp; Holloway (2010) filter set.",
            S["body"]
        ),
        sp(0.3),
    ]

    # One row per unique compound
    seen = set()
    unique = []
    for r in results:
        if r["ligand"] not in seen and "admet_decision" in r:
            seen.add(r["ligand"])
            unique.append(r)

    admet_header = ["Compound", "MW", "logP", "HBD", "HBA",
                    "TPSA", "logBB", "SA", "PAINS", "ADMET"]
    admet_rows = [admet_header]
    for r in unique:
        admet_rows.append([
            r["ligand"],
            str(r.get("mw", "N/A")),
            str(r.get("logp", "N/A")),
            str(r.get("hbd", "N/A")),
            str(r.get("hba", "N/A")),
            f"{r.get('tpsa', 'N/A')} A²",
            str(r.get("logbb", "N/A")),
            str(r.get("sa_score", "N/A")),
            str(r.get("pains_count", "N/A")),
            r.get("admet_decision", "N/A"),
        ])

    col_w2 = [3.2*cm, 1.5*cm, 1.5*cm, 1.2*cm, 1.2*cm,
              2.0*cm, 1.8*cm, 1.5*cm, 1.5*cm, 1.6*cm]
    admet_table = Table(admet_rows, colWidths=col_w2, repeatRows=1)

    ts2 = [
        ("BACKGROUND",  (0,0), (-1,0), _rgb(DARK_BLUE)),
        ("TEXTCOLOR",   (0,0), (-1,0), _rgb(WHITE)),
        ("FONTNAME",    (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE",    (0,0), (-1,-1), 8),
        ("ALIGN",       (0,0), (-1,-1), "CENTER"),
        ("ALIGN",       (0,0), (0,-1), "LEFT"),
        ("VALIGN",      (0,0), (-1,-1), "MIDDLE"),
        ("ROWHEIGHT",   (0,0), (-1,-1), 15),
        ("BOX",         (0,0), (-1,-1), 0.5, _rgb(MID_BLUE)),
        ("INNERGRID",   (0,0), (-1,-1), 0.25, _rgb(MID_GREY)),
        ("LEFTPADDING", (0,0), (-1,-1), 5),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [_rgb(LIGHT_GREY), _rgb(WHITE)]),
    ]
    for i, r in enumerate(unique, start=1):
        if r.get("admet_decision") == "PASS":
            ts2 += [
                ("TEXTCOLOR", (-1,i), (-1,i), _rgb(ACCENT)),
                ("FONTNAME",  (-1,i), (-1,i), "Helvetica-Bold"),
            ]
        else:
            ts2 += [("TEXTCOLOR", (-1,i), (-1,i), _rgb(RED))]

    admet_table.setStyle(TableStyle(ts2))
    story += [admet_table, sp(0.5)]

    # ADMET failure details
    story.append(Paragraph("Failure Details", S["h2"]))
    for r in unique:
        if r.get("admet_decision") == "FAIL" and r.get("admet_failures"):
            failures = r["admet_failures"]
            if failures and failures != "None":
                story.append(Paragraph(
                    f"<b>{r['ligand']}:</b> {failures}",
                    S["body"]
                ))

    story.append(PageBreak())

    # =========================================================================
    # PAGE 5 — METHODS
    # =========================================================================
    story += [
        Paragraph("Methods &amp; Pipeline Overview", S["h1"]),
        hr(),
        Paragraph("Protein Preparation", S["h2"]),
        Paragraph(
            "Protein structures were downloaded from the RCSB Protein Data Bank. "
            "Crystallographic water molecules and co-crystallised ligands were removed. "
            "Structures were converted to PDBQT format using Open Babel 3.1 with "
            "Gasteiger partial charge assignment.",
            S["body"]
        ),
        Paragraph("Ligand Preparation", S["h2"]),
        Paragraph(
            "Ligands were prepared from SMILES notation using RDKit 2026. "
            "3D conformers were generated using the ETKDGv3 algorithm with "
            "MMFF94 force field geometry optimisation. "
            "Conversion to PDBQT format was performed with Open Babel.",
            S["body"]
        ),
        Paragraph("Molecular Docking", S["h2"]),
        Paragraph(
            f"Molecular docking was performed using AutoDock Vina v1.2.5. "
            f"Exhaustiveness was set to {16} with {9} binding modes per run. "
            f"Compounds with binding affinity ≤ {cutoff} kcal/mol were considered hits.",
            S["body"]
        ),
        Paragraph("ADMET Filtering", S["h2"]),
        Paragraph(
            "ADMET properties were calculated using RDKit descriptors. "
            "Blood-brain barrier permeability was estimated using the Clark (1999) "
            "logBB model. PAINS (Pan-Assay Interference Compounds) screening used "
            "the Baell &amp; Holloway (2010) filter catalogue. "
            "CNS-specific thresholds were applied: MW ≤ 450 Da, TPSA &lt; 90 A², "
            "logBB &gt; -1.0, HBD ≤ 3, rotatable bonds ≤ 8.",
            S["body"]
        ),
        Paragraph("Targets", S["h2"]),
    ]

    target_data = [["PDB ID", "Target", "Disease Relevance"]]
    for t in targets:
        target_data.append([t["pdb_id"], t["name"], t["note"]])

    t_table = Table(target_data, colWidths=[2*cm, 5*cm, 9.5*cm])
    t_table.setStyle(TableStyle([
        ("BACKGROUND",  (0,0), (-1,0), _rgb(DARK_BLUE)),
        ("TEXTCOLOR",   (0,0), (-1,0), _rgb(WHITE)),
        ("FONTNAME",    (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE",    (0,0), (-1,-1), 8),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [_rgb(LIGHT_GREY), _rgb(WHITE)]),
        ("VALIGN",      (0,0), (-1,-1), "TOP"),
        ("ROWHEIGHT",   (0,0), (-1,-1), 22),
        ("BOX",         (0,0), (-1,-1), 0.5, _rgb(MID_BLUE)),
        ("INNERGRID",   (0,0), (-1,-1), 0.3, _rgb(MID_GREY)),
        ("LEFTPADDING", (0,0), (-1,-1), 6),
        ("WORDWRAP",    (0,0), (-1,-1), True),
    ]))
    story += [t_table, sp(0.5)]

    story += [
        Paragraph("Software &amp; Tools", S["h2"]),
        Paragraph(
            "AutoDock Vina 1.2.5 &nbsp;|&nbsp; "
            "Open Babel 3.1 &nbsp;|&nbsp; "
            "RDKit 2026 &nbsp;|&nbsp; "
            "Python 3.12 &nbsp;|&nbsp; "
            "ReportLab (PDF generation)",
            S["body"]
        ),
        sp(1),
        hr(MID_GREY, 0.5),
        Paragraph(
            f"Generated by NeuroDock Pipeline &nbsp;|&nbsp; "
            f"{datetime.now().strftime('%d %B %Y')}",
            S["small"]
        ),
    ]

    # ── Build PDF ─────────────────────────────────────────────────────────────
    doc.build(story)
    log.info(f"  [pdf_report] PDF saved: {pdf_path}")
    return pdf_path


def generate_pdf_report(
    results      : List[Dict],
    results_dir  : Path,
    cutoff       : float,
    targets      : List[Dict],
    ligands      : List[Dict],
    interactions : list = None,
) -> None:
    """
    Main function called by run_pipeline.py.
    """
    log.info(f"\n{'═'*50}")
    log.info(f"  PDF REPORT GENERATION")
    log.info(f"{'═'*50}")

    if not results:
        log.warning("  [pdf_report] No results — skipping PDF generation")
        return

    try:
        pdf_path = build_pdf(results, results_dir, cutoff, targets, ligands, interactions)
        print(f"\n  📄 PDF report saved: {pdf_path}\n")
    except Exception as e:
        log.error(f"  [pdf_report] Failed: {e}")
        raise