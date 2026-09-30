"""Transparent, rule-based prioritisation inspired by ACMG/AMP evidence types.

IMPORTANT: this is an *educational prioritisation*, not a clinical
classification. It implements a handful of easy-to-audit rules; real
interpretation needs the full ACMG/AMP framework, gene-disease context,
segregation, functional and case data, and review by a qualified professional.

Notes on the evidence used:
* BA1 / BS1: allele frequency far too high for a rare Mendelian disease.
* PM2 (supporting): absent from / very rare in population data.
* LoF (PVS1-like): a loss-of-function consequence. Real PVS1 also requires that
  loss of function is a known disease mechanism for the gene, so we only flag it.
* ClinVar is shown as its own evidence line instead of PP5/BP6, because ClinGen's
  Sequence Variant Interpretation group recommended retiring PP5/BP6.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from .annotate import Annotation

LOF_TERMS = {
    "stop_gained",
    "nonsense",
    "frameshift_variant",
    "frameshift",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "initiator_codon_variant",
}

TIER_ORDER = {"high": 0, "review": 1, "vus": 2, "low": 3}
TIER_LABEL = {
    "high": "High priority",
    "review": "Needs manual review",
    "vus": "Uncertain / unannotated",
    "low": "Low priority",
}


@dataclass
class Thresholds:
    ba1: float = 0.05  # stand-alone benign frequency
    bs1: float = 0.01  # too common for a rare disorder
    pm2: float = 0.0001  # "very rare" cut-off


@dataclass
class Result:
    tier: str
    evidence: List[str] = field(default_factory=list)

    @property
    def label(self) -> str:
        return TIER_LABEL[self.tier]


def is_lof(consequences: List[str]) -> bool:
    return any(c in LOF_TERMS for c in consequences)


def prioritise(ann: Annotation, th: Thresholds = Thresholds()) -> Result:
    af = ann.allele_freq
    cv = ann.clinvar
    evidence: List[str] = []

    if af is None:
        evidence.append("Frequency: not found in population data")
    else:
        evidence.append(f"Frequency: {af:.6g}")

    rare = af is None or af <= th.pm2
    common_ba1 = af is not None and af >= th.ba1
    common_bs1 = af is not None and af >= th.bs1
    lof = is_lof(ann.consequences)

    if cv:
        stars = "*" * cv.stars if cv.stars else "no stars"
        evidence.append(f"ClinVar: {cv.significance.replace('_', ' ')} ({stars})")
    if lof:
        evidence.append("Loss-of-function consequence (PVS1 candidate; check gene mechanism)")
    if rare:
        evidence.append(f"PM2_supporting: absent or AF <= {th.pm2:g}")
    if common_ba1:
        evidence.append(f"BA1: AF >= {th.ba1:g} (stand-alone benign)")
    elif common_bs1:
        evidence.append(f"BS1: AF >= {th.bs1:g} (too common for a rare disorder)")

    # 1) Very common variants are benign for rare-disease purposes, whatever else is said.
    if common_ba1:
        return Result("low", evidence)

    # 2) ClinVar assertion drives the tier when present.
    if cv:
        if cv.category == "pathogenic":
            if common_bs1:
                evidence.append("Conflict: ClinVar pathogenic but frequency is high; review")
                return Result("review", evidence)
            return Result("high" if cv.stars >= 2 else "review", evidence)
        if cv.category == "conflicting":
            return Result("review", evidence)
        if cv.category == "benign":
            return Result("low" if cv.stars >= 1 else "vus", evidence)
        # ClinVar VUS / other falls through to frequency and consequence rules.

    # 3) Novel or unreviewed loss-of-function that is rare deserves a human look.
    if lof and rare:
        return Result("review", evidence)

    # 4) Too common for a rare disorder.
    if common_bs1:
        return Result("low", evidence)

    return Result("vus", evidence)
