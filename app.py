import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import sys

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="NeuroDock — MAO-B Virtual Screening",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main-title {
        font-size: 2.4rem; font-weight: 800; color: #1F3A5F;
        text-align: center; margin-bottom: 0.2rem;
    }
    .subtitle {
        font-size: 1.1rem; color: #555; text-align: center; margin-bottom: 2rem;
    }
    .metric-card {
        background: #F0F4FF; border-radius: 12px; padding: 1rem;
        border-left: 5px solid #2563EB; margin-bottom: 1rem;
    }
    .lead-card {
        background: #DCEFE4; border-radius: 12px; padding: 1rem;
        border-left: 5px solid #059669; margin-bottom: 0.5rem;
    }
    .high-risk { color: #DC2626; font-weight: bold; }
    .mod-high-risk { color: #D97706; font-weight: bold; }
    .mod-risk { color: #CA8A04; font-weight: bold; }
    .low-risk { color: #059669; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.image("https://img.shields.io/badge/Python-3.10+-blue.svg")
    st.image("https://img.shields.io/badge/License-MIT-green.svg")
    st.markdown("---")
    st.markdown("### Navigation")
    page = st.radio("", [
        "Overview",
        "Docking Results",
        "ADMET Analysis",
        "Nitrosamine Risk",
        "Figures",
        "Run Pipeline",
    ])
    st.markdown("---")
    st.markdown("**Author:** Thanuhith J.")
    st.markdown("**Institution:** Dr. M.G.R. Educational and Research Institute")
    st.markdown("**Target:** MAO-B (PDB: 2V61)")
    st.markdown("[GitHub Repo](https://github.com/ThanuHith/neuro_dock)")

# ── Data ─────────────────────────────────────────────────────────────────────
leads_data = {
    "Compound":    ["Donepezil", "Chrysin", "Resveratrol", "Safinamide", "Galantamine", "Rivastigmine"],
    "Rank":        [4, 12, 16, 20, 22, 23],
    "dG (kcal/mol)": [-8.844, -8.484, -8.352, -8.134, -7.744, -7.124],
    "logBB":       [0.287, -0.470, -0.307, -0.453, -0.200, 0.073],
    "TPSA (A2)":   [38.77, 70.67, 60.69, 64.35, 41.93, 32.78],
    "MW (Da)":     [393.5, 254.2, 228.3, 302.4, 287.4, 250.3],
    "logP":        [4.75, 2.87, 2.97, 2.37, 1.85, 2.76],
    "SA Score":    [2.71, 2.10, 2.11, 2.27, 4.23, 2.67],
    "Approved Drug": ["No", "No", "No", "Yes", "Yes", "Yes"],
    "Nitrosamine Risk": ["MODERATE-HIGH", "LOW", "LOW", "MODERATE", "MODERATE", "HIGH"],
}
df_leads = pd.DataFrame(leads_data)

full_data = {
    "Rank": list(range(1, 30)),
    "Compound": [
        "Hesperetin","Hesperidin","Rutin","Donepezil","Luteolin",
        "Quercetin","Fisetin","Myricetin","Baicalein","Genistein",
        "Kaempferol","Chrysin","Apigenin","Catechin","Epicatechin",
        "Resveratrol","7,8-DHF","7,8-Dihydroxyflavone","Curcumin","Safinamide",
        "Epigallocatechin gallate","Galantamine","Rivastigmine","Naringenin","NP-alpha1",
        "Rasagiline","Caffeic acid","Ferulic acid","Selegiline"
    ],
    "dG (kcal/mol)": [
        -9.02,-8.975,-8.92,-8.844,-8.793,
        -8.692,-8.654,-8.606,-8.591,-8.541,
        -8.537,-8.484,-8.471,-8.418,-8.366,
        -8.352,-8.315,-8.281,-8.199,-8.134,
        -8.036,-7.744,-7.124,-6.905,-6.755,
        -6.618,-6.357,-6.263,-5.787
    ],
    "Final Decision": [
        "ADMET Fail","ADMET Fail","ADMET Fail","LEAD","ADMET Fail",
        "ADMET Fail","ADMET Fail","ADMET Fail","ADMET Fail","ADMET Fail",
        "ADMET Fail","LEAD","ADMET Fail","ADMET Fail","ADMET Fail",
        "LEAD","ADMET Fail","ADMET Fail","ADMET Fail","LEAD",
        "ADMET Fail","LEAD","LEAD","Below threshold","Below threshold",
        "Below threshold","Below threshold","Below threshold","Below threshold"
    ],
}
df_full = pd.DataFrame(full_data)

nitro_data = {
    "Compound":    ["Rivastigmine","Donepezil","Safinamide","Galantamine","Chrysin","Resveratrol"],
    "Risk Level":  ["HIGH","MODERATE-HIGH","MODERATE","MODERATE","LOW","LOW"],
    "Risk Score":  [4, 3, 2, 2, 0, 0],
    "Reactive N":  [1, 1, 1, 1, 0, 0],
    "Structural Basis": [
        "N,N-dimethyl tertiary amine (NDMA-type precursor)",
        "Piperidine ring N — benzylpiperidine class flagged EMA 2022",
        "Secondary amine (R2NH)",
        "7-membered ring tertiary amine",
        "No nitrogen atoms",
        "No nitrogen atoms",
    ],
    "ICH M7 Action": [
        "Mandatory confirmatory testing, TTC=18 ng/day",
        "Stepwise risk assessment + in vitro study",
        "Assess excipients, limit nitrite sources",
        "Assess excipients, limit nitrite sources",
        "No action required",
        "No action required",
    ],
}
df_nitro = pd.DataFrame(nitro_data)


# ── PAGES ─────────────────────────────────────────────────────────────────────

if page == "Overview":
    st.markdown('<div class="main-title">NeuroDock</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">AI-Driven MAO-B Virtual Screening Pipeline | Computational Drug Discovery</div>', unsafe_allow_html=True)

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Compounds Screened", "29")
    col2.metric("Docking Hits", "23", "79.3% hit rate")
    col3.metric("Lead Candidates", "6", "passed all ADMET")
    col4.metric("Best Binding Energy", "-8.84 kcal/mol", "Donepezil")

    st.markdown("---")
    col_a, col_b = st.columns(2)

    with col_a:
        st.markdown("### What is NeuroDock?")
        st.markdown("""
NeuroDock is a fully open-source, end-to-end virtual screening pipeline for
computational drug discovery targeting neurodegeneration.

**Target:** Monoamine Oxidase-B (MAO-B), a key enzyme elevated in Parkinson's
and Alzheimer's disease. MAO-B inhibition is an established therapeutic strategy.

**Approach:** 29 compounds (approved CNS drugs + dietary polyphenols) were screened
against the MAO-B crystal structure (PDB: 2V61) using molecular docking, followed
by strict CNS-focused ADMET filtering and nitrosamine safety assessment.
        """)

    with col_b:
        st.markdown("### Pipeline Steps")
        steps = [
            ("1", "Protein Preparation", "Download PDB, clean, convert to PDBQT"),
            ("2", "Ligand Library", "PubChem API retrieval + CID verification"),
            ("3", "3D Conformer Generation", "RDKit ETKDGv3 + MMFF94 minimization"),
            ("4", "Molecular Docking", "AutoDock Vina 1.2.5, exhaustiveness=16"),
            ("5", "Binding Energy Filter", "Threshold: ≤ -7.0 kcal/mol"),
            ("6", "ADMET Filtering", "CNS-specific: TPSA, logBB, PAINS, Lipinski"),
            ("7", "Interaction Analysis", "H-bonds, hydrophobic, pi-stacking"),
            ("8", "Nitrosamine Risk", "ICH M7(R2) / EMA-FDA 2023 framework"),
            ("9", "Results", "Figures, CSV, PDF report"),
        ]
        for num, title, desc in steps:
            st.markdown(f"**{num}.** {title} — *{desc}*")

    st.markdown("---")
    st.markdown("### Key Finding")
    st.info("Safinamide (an approved MAO-B inhibitor) and two AChE inhibitors (Galantamine, Rivastigmine) were correctly recovered from the blinded compound library — internally validating the pipeline's ability to identify known actives.")


elif page == "Docking Results":
    st.markdown("## Docking Results — All 29 Compounds")

    col1, col2, col3 = st.columns(3)
    col1.metric("Total Screened", "29")
    col2.metric("Docking Hits (≤-7.0)", "23")
    col3.metric("Lead Candidates", "6")

    st.markdown("### Full Ranked Results")
    def color_decision(val):
        if val == "LEAD":
            return "background-color: #DCEFE4; color: #059669; font-weight: bold"
        elif val == "Below threshold":
            return "color: #94A3B8"
        else:
            return "color: #2563EB"
    st.dataframe(
        df_full.style.map(color_decision, subset=["Final Decision"]),
        use_container_width=True, height=600
    )

    st.markdown("### Binding Affinity Distribution")
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ["#059669" if d=="LEAD" else "#2563EB" if d=="ADMET Fail" else "#94A3B8"
              for d in df_full["Final Decision"]]
    bars = ax.barh(df_full["Compound"], df_full["dG (kcal/mol)"],
                   color=colors, edgecolor="white", height=0.7)
    ax.axvline(x=-7.0, color="black", linestyle="--", linewidth=1.2, label="Threshold (-7.0 kcal/mol)")
    ax.set_xlabel("Binding Affinity dG (kcal/mol)", fontweight="bold")
    ax.set_title("Molecular Docking Results — MAO-B Virtual Screen", fontweight="bold", pad=12)
    ax.invert_yaxis()
    from matplotlib.patches import Patch
    legend = [Patch(color="#059669", label="Lead candidate"),
              Patch(color="#2563EB", label="Docking hit (ADMET fail)"),
              Patch(color="#94A3B8", label="Below threshold")]
    ax.legend(handles=legend, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("#F8FAFC")
    fig.patch.set_facecolor("white")
    st.pyplot(fig)


elif page == "ADMET Analysis":
    st.markdown("## ADMET Analysis — 6 Lead Candidates")
    st.markdown("All six lead candidates satisfied every CNS-specific ADMET criterion applied.")

    for _, row in df_leads.iterrows():
        with st.expander(f"**{row['Compound']}** — dG = {row['dG (kcal/mol)']} kcal/mol"):
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Binding Energy", f"{row['dG (kcal/mol)']} kcal/mol")
            c2.metric("logBB", f"{row['logBB']:.3f}")
            c3.metric("TPSA", f"{row['TPSA (A2)']} A2")
            c4.metric("logP", f"{row['logP']}")
            c1.metric("MW", f"{row['MW (Da)']} Da")
            c2.metric("SA Score", f"{row['SA Score']}")
            c3.metric("Approved Drug", row["Approved Drug"])
            c4.metric("Nitrosamine Risk", row["Nitrosamine Risk"])

    st.markdown("---")
    st.markdown("### Comparative Table")
    st.dataframe(df_leads, use_container_width=True)

    st.markdown("### logP vs Binding Affinity")
    fig, ax = plt.subplots(figsize=(8, 5))
    colors_lead = ["#059669"]*6
    ax.scatter(df_leads["logP"], df_leads["dG (kcal/mol)"],
               c=colors_lead, s=120, marker="D", edgecolor="white", zorder=3)
    for _, row in df_leads.iterrows():
        ax.annotate(row["Compound"],
                    (row["logP"], row["dG (kcal/mol)"]),
                    fontsize=8.5, xytext=(6, 4), textcoords="offset points",
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="#059669", linewidth=0.8))
    ax.axhline(y=-7.0, color="black", linestyle="--", alpha=0.5)
    ax.set_xlabel("Calculated logP", fontweight="bold")
    ax.set_ylabel("Binding Affinity dG (kcal/mol)", fontweight="bold")
    ax.set_title("Lipophilicity vs Binding Affinity — Lead Candidates", fontweight="bold")
    ax.set_facecolor("#F8FAFC")
    fig.patch.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    st.pyplot(fig)


elif page == "Nitrosamine Risk":
    st.markdown("## Nitrosamine Formation Risk Assessment")
    st.markdown("**Framework:** ICH M7(R2) / EMA-FDA Nitrosamine Impurities Guidance 2021/2023")
    st.markdown("**Method:** In silico structural alert screening via RDKit")

    st.warning("This is a preliminary structural screen only. Confirmatory in vitro assays are required for regulatory submissions under EMA/FDA 2023 guidance.")

    risk_colors = {
        "HIGH": "#FEE2E2",
        "MODERATE-HIGH": "#FEF3C7",
        "MODERATE": "#FEF9C3",
        "LOW": "#DCFCE7",
    }
    risk_text_colors = {
        "HIGH": "#DC2626",
        "MODERATE-HIGH": "#D97706",
        "MODERATE": "#CA8A04",
        "LOW": "#059669",
    }

    for _, row in df_nitro.iterrows():
        risk = row["Risk Level"]
        bg = risk_colors.get(risk, "#F8FAFC")
        tc = risk_text_colors.get(risk, "#000")
        st.markdown(f"""
        <div style="background:{bg}; border-radius:10px; padding:1rem; margin-bottom:0.8rem;
                    border-left:5px solid {tc};">
            <b style="font-size:1.1rem;">{row['Compound']}</b>
            <span style="float:right; color:{tc}; font-weight:bold; font-size:1rem;">{risk}</span><br>
            <small><b>Structural basis:</b> {row['Structural Basis']}</small><br>
            <small><b>Regulatory action:</b> {row['ICH M7 Action']}</small>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### Risk Summary Chart")
    fig, ax = plt.subplots(figsize=(8, 4))
    bar_colors = ["#DC2626","#D97706","#CA8A04","#CA8A04","#059669","#059669"]
    bars = ax.barh(df_nitro["Compound"], df_nitro["Risk Score"],
                   color=bar_colors, edgecolor="white", height=0.6)
    ax.set_xlabel("Risk Score (0=LOW, 4=HIGH)", fontweight="bold")
    ax.set_title("Nitrosamine Formation Risk — ICH M7(R2)", fontweight="bold")
    ax.set_xlim(0, 4.5)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(["LOW", "LOW-MOD", "MODERATE", "MOD-HIGH", "HIGH"])
    ax.set_facecolor("#F8FAFC")
    fig.patch.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    st.pyplot(fig)

    st.markdown("### Key Takeaway")
    st.success("Chrysin and Resveratrol — the two novel natural polyphenol leads — carry LOW nitrosamine risk with no reactive nitrogen atoms, making them the most regulatory-favourable candidates for further development.")


elif page == "Figures":
    st.markdown("## Publication Figures")
    st.markdown("All figures generated at 300 DPI using Matplotlib from real pipeline output data.")

    fig_dir = "data/results/figures"
    figures = {
        "fig1_binding_energy_bar.png": ("Figure 1", "Ranked binding affinities (dG) of all 29 screened compounds. Green = lead candidates, Blue = docking hits, Grey = below threshold."),
        "fig2_admet_radar.png": ("Figure 2", "CNS drug-likeness radar plot for all 6 lead candidates across 6 normalised ADMET dimensions."),
        "fig3_logp_vs_affinity_scatter.png": ("Figure 3", "Lipophilicity (logP) vs binding affinity scatter plot. Green diamonds = lead candidates."),
        "fig4_admet_heatmap.png": ("Figure 4", "Pass/fail matrix across all 9 screening criteria for all 29 compounds."),
        "fig5_pipeline_flowchart.png": ("Figure 5", "NeuroDock 9-stage computational pipeline overview."),
    }

    for fname, (title, caption) in figures.items():
        fpath = os.path.join(fig_dir, fname)
        if os.path.exists(fpath):
            st.markdown(f"### {title}")
            st.image(fpath, caption=caption, use_container_width=True)
        else:
            st.info(f"{title} — Run the pipeline first to generate: `{fpath}`")

    st.markdown("---")
    st.markdown("*Figures are generated automatically by the pipeline in `data/results/figures/`.*")


elif page == "Run Pipeline":
    st.markdown("## Run NeuroDock Pipeline")
    st.markdown("This app displays pre-computed results. To generate fresh results, run the pipeline locally.")

    st.markdown("### Quick Start")
    st.code("""# 1. Clone the repository
git clone https://github.com/ThanuHith/neuro_dock.git
cd neuro_dock

# 2. Create virtual environment
python3 -m venv venv
source venv/bin/activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Install AutoDock Vina + Open Babel
bash scripts/install_tools.sh

# 5. Run the pipeline
python3 run_pipeline.py""", language="bash")

    st.markdown("### System Requirements")
    req_data = {
        "Component": ["Python", "AutoDock Vina", "Open Babel", "RAM", "Storage", "OS"],
        "Requirement": ["3.10+", "1.2.5", "3.1+", "8 GB minimum", "5 GB free", "Ubuntu 20.04+ or WSL2"],
    }
    st.table(pd.DataFrame(req_data))

    st.markdown("### Expected Output")
    st.markdown("""
After running `python3 run_pipeline.py`, you will find in `data/results/`:
- `docking_results.csv` — Full ranked results for all 29 compounds
- `leads_only.csv` — The 6 lead candidates with all ADMET data
- `nitrosamine_risk_report.txt` — ICH M7 risk assessment for leads
- `NeuroDock_Report.pdf` — 6-page PDF pipeline report
- `figures/` — 5 publication-quality PNG figures at 300 DPI
    """)

    st.markdown("---")
    st.markdown("### Citation")
    st.code("""Thanuhith J. (2025). In Silico Identification of Donepezil and Natural
Polyphenols as Monoamine Oxidase-B (MAO-B) Inhibitors: A Structure-Based
Virtual Screening, ADMET Profiling, and Protein-Ligand Interaction Study.
Dr. M.G.R. Educational and Research Institute, Chennai, India.
GitHub: https://github.com/ThanuHith/neuro_dock""")
