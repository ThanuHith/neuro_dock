# modules/nitrosamine.py
# ICH M7(R2) / EMA-FDA Nitrosamine Guidance 2021/2023
import logging
from typing import Dict, List, Tuple
log = logging.getLogger(__name__)

REG_NOTES = {
    "HIGH": ("Mandatory confirmatory testing per EMA/FDA 2023. "
             "N,N-dimethyl amine = ICH M7 Cohort of Concern (NDMA-type). TTC=18 ng/day."),
    "MODERATE-HIGH": ("Piperidine/benzylpiperidine class flagged in EMA 2022. "
                      "Stepwise risk assessment + in vitro study recommended."),
    "MODERATE": ("Secondary or ring tertiary amine. Nitrosamine formation possible "
                 "under acidic pH + nitrite. Assess excipients per ICH M7(R2)."),
    "LOW-MODERATE": "Amine present, lower reactivity. Screen for nitrite excipients.",
    "LOW": "No reactive amine nitrogen. No nitrosamine formation pathway.",
}
SCORE_TO_LEVEL = {4:"HIGH", 3:"MODERATE-HIGH", 2:"MODERATE", 1:"LOW-MODERATE", 0:"LOW"}


def _result(compound, risk_level, risk_score, reactive_nitrogens, rationale, regulatory_note):
    return dict(compound=compound, risk_level=risk_level, risk_score=risk_score,
                reactive_nitrogens=reactive_nitrogens, rationale=rationale,
                regulatory_note=regulatory_note)


def assess_nitrosamine_risk(name: str, smiles: str) -> Dict:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return _result(name, "UNKNOWN", -1, [], "Cannot parse SMILES", "")

    nitrogens = [a for a in mol.GetAtoms() if a.GetSymbol() == "N"]
    if not nitrogens:
        return _result(name, "LOW", 0, [],
                       "No nitrogen atoms — no nitrosamine formation pathway.",
                       REG_NOTES["LOW"])

    reactive = []
    for atom in nitrogens:
        idx, h, in_ring = atom.GetIdx(), atom.GetTotalNumHs(), atom.IsInRing()
        ring_size = None
        if in_ring:
            for ring in mol.GetRingInfo().AtomRings():
                if idx in ring:
                    ring_size = len(ring)
                    break

        # Skip amide/carbamate N — deactivated by adjacent C=O
        adj_c = False
        for nb in atom.GetNeighbors():
            nba = mol.GetAtomWithIdx(nb.GetIdx())
            if nba.GetSymbol() == "C":
                for nb2 in nba.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(nba.GetIdx(), nb2.GetIdx())
                    if (mol.GetAtomWithIdx(nb2.GetIdx()).GetSymbol() == "O"
                            and bond.GetBondTypeAsDouble() == 2.0):
                        adj_c = True
                        break
        if adj_c:
            continue

        methyls = sum(
            1 for nb in atom.GetNeighbors()
            if (mol.GetAtomWithIdx(nb.GetIdx()).GetSymbol() == "C"
                and mol.GetAtomWithIdx(nb.GetIdx()).GetDegree() == 1)
        )

        if methyls >= 2 and h == 0:
            cls, risk, score = "N,N-dimethyl tertiary amine", "HIGH", 4
            note = "NDMA-type precursor. ICH M7 Cohort of Concern. TTC=18 ng/day."
        elif in_ring and h == 0 and ring_size == 6:
            cls, risk, score = "Piperidine/piperazine ring N", "MODERATE-HIGH", 3
            note = "6-membered ring amine. Benzylpiperidine class flagged EMA 2022."
        elif in_ring and h == 0:
            cls, risk, score = f"Ring tertiary amine ({ring_size}-membered)", "MODERATE", 2
            note = "Ring-constrained tertiary amine. Steric hindrance reduces but does not eliminate risk."
        elif h == 1:
            cls, risk, score = "Secondary amine (R2NH)", "MODERATE", 2
            note = "Secondary amines form nitrosamines under acidic pH + nitrite exposure."
        elif h == 0 and not in_ring:
            cls, risk, score = "Acyclic tertiary amine (R3N)", "LOW-MODERATE", 1
            note = "Lower direct risk; assess dealkylation potential under acidic conditions."
        else:
            cls, risk, score = "Primary amine (RNH2)", "LOW-MODERATE", 1
            note = "Primary nitrosamines weaker carcinogens; still ICH M7 flagged."

        reactive.append(dict(atom_idx=idx, cls=cls, risk=risk, score=score,
                             note=note, in_ring=in_ring, ring_size=ring_size,
                             h_count=h, methyl_count=methyls))

    if not reactive:
        return _result(name, "LOW", 0, [],
                       "All N atoms are amide/carbamate type — deactivated by carbonyl. Not precursors.",
                       "Low regulatory concern.")

    max_score = max(n["score"] for n in reactive)
    overall = SCORE_TO_LEVEL[max_score]
    rationale = (f"{len(reactive)} reactive N: " +
                 "; ".join(f"{n['cls']} at atom[{n['atom_idx']}] ({n['risk']})"
                           for n in reactive))
    return _result(name, overall, max_score, reactive, rationale,
                   REG_NOTES.get(overall, "Low regulatory concern."))


def run_nitrosamine_assessment(compounds: List[Tuple[str, str]]) -> List[Dict]:
    log.info("=" * 50)
    log.info(f"  NITROSAMINE RISK ASSESSMENT  ({len(compounds)} compounds)")
    log.info("  Framework : ICH M7(R2) / EMA-FDA 2023 Guidance")
    log.info("=" * 50)
    results = []
    for name, smiles in compounds:
        r = assess_nitrosamine_risk(name, smiles)
        results.append(r)
        log.info(f"  [{r['risk_level']:<14}] {name:<22} "
                 f"({len(r['reactive_nitrogens'])} reactive N)")
    results.sort(key=lambda x: -x["risk_score"])
    h  = sum(1 for r in results if r["risk_score"] >= 4)
    mh = sum(1 for r in results if r["risk_score"] == 3)
    m  = sum(1 for r in results if r["risk_score"] == 2)
    lo = sum(1 for r in results if r["risk_score"] <= 1)
    log.info(f"  Summary: HIGH={h}, MOD-HIGH={mh}, MOD={m}, LOW={lo}")
    log.info("  DISCLAIMER: Structural screen only. In vitro assays required for regulatory submissions.")
    return results


def format_report(results: List[Dict]) -> str:
    lines = ["=" * 70,
             "  NITROSAMINE FORMATION RISK ASSESSMENT REPORT",
             "  Reference: ICH M7(R2) / EMA-FDA Nitrosamine Guidance 2023",
             "=" * 70]
    for r in results:
        lines += [f"\nCompound  : {r['compound']}",
                  f"Risk Level: {r['risk_level']}",
                  f"Rationale : {r['rationale']}",
                  f"Reg. Note : {r['regulatory_note']}"]
        for n in r["reactive_nitrogens"]:
            lines.append(f"  atom[{n['atom_idx']}] {n['cls']}: {n['risk']} | {n['note']}")
    lines += ["", "-" * 70, "SUMMARY TABLE", "-" * 70,
              f"{'Compound':<22} {'Risk Level':<18} {'React.N':<9} Rationale",
              "-" * 70]
    for r in results:
        short = r["rationale"][:38] + "..." if len(r["rationale"]) > 38 else r["rationale"]
        lines.append(f"{r['compound']:<22} {r['risk_level']:<18} "
                     f"{len(r['reactive_nitrogens']):<9} {short}")
    lines += ["",
              "DISCLAIMER: In silico structural screening only. Does not replace",
              "confirmatory in vitro assays required per EMA/FDA 2023 guidance.",
              "=" * 70]
    return "\n".join(lines)
