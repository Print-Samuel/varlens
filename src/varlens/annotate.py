"""Attach ClinVar assertions, population frequency and consequence to variants."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union
from urllib.parse import unquote

from .vcf import Key, Variant, read_vcf


@dataclass
class ClinVarRecord:
    significance: str  # raw CLNSIG, e.g. "Pathogenic/Likely_pathogenic"
    category: str  # pathogenic | benign | vus | conflicting | other
    review_status: str  # raw CLNREVSTAT
    stars: int  # 0-4, ClinVar review-status stars
    gene: str
    disease: str
    consequences: List[str]


@dataclass
class Annotation:
    variant: Variant
    gene: str
    consequences: List[str]
    allele_freq: Optional[float]
    clinvar: Optional[ClinVarRecord]


def clinvar_category(significance: str) -> str:
    sig = significance.lower()
    if "conflicting" in sig:
        return "conflicting"
    if "benign" in sig and "pathogenic" not in sig:
        return "benign"
    if "pathogenic" in sig and "benign" not in sig:
        return "pathogenic"
    if "uncertain" in sig:
        return "vus"
    return "other"


def review_stars(review_status: str) -> int:
    """Map CLNREVSTAT to the 0-4 star scale ClinVar shows on its website."""
    rs = review_status.lower()
    if rs.startswith("no_assertion") or rs.startswith("no_classification"):
        return 0  # must come first: "no_assertion_criteria_provided" contains "criteria_provided"
    if "practice_guideline" in rs:
        return 4
    if "reviewed_by_expert_panel" in rs:
        return 3
    if "criteria_provided" in rs and "multiple_submitters" in rs and "no_conflicts" in rs:
        return 2
    if "criteria_provided" in rs:  # single submitter, or conflicting classifications
        return 1
    return 0


def _gene_from_geneinfo(geneinfo: str) -> str:
    # "BRCA1:672|BRCA1-AS1:123" -> "BRCA1"
    return geneinfo.split("|")[0].split(":")[0] if geneinfo else ""


def _consequences_from_mc(mc: str) -> List[str]:
    # "SO:0001587|nonsense,SO:0001583|missense_variant" -> ["nonsense", "missense_variant"]
    terms = []
    for item in mc.split(",") if mc else []:
        _, _, name = item.partition("|")
        if name:
            terms.append(name.strip().lower())
    return terms


def load_clinvar(path: Union[str, Path], wanted: Iterable[Key]) -> Dict[Key, ClinVarRecord]:
    """Stream a ClinVar VCF and keep only records for the variants of interest."""
    wanted = set(wanted)
    records: Dict[Key, ClinVarRecord] = {}
    for v in read_vcf(path, wanted=wanted):
        sig = unquote(v.info.get("CLNSIG", ""))
        rev = unquote(v.info.get("CLNREVSTAT", ""))
        records[v.key] = ClinVarRecord(
            significance=sig,
            category=clinvar_category(sig),
            review_status=rev,
            stars=review_stars(rev),
            gene=_gene_from_geneinfo(v.info.get("GENEINFO", "")),
            disease=unquote(v.info.get("CLNDN", "")).replace("_", " "),
            consequences=_consequences_from_mc(v.info.get("MC", "")),
        )
    return records


def load_frequencies(
    path: Union[str, Path], wanted: Iterable[Key], af_field: str = "AF"
) -> Dict[Key, float]:
    """Stream a population-frequency VCF (e.g. a gnomAD sites subset)."""
    wanted = set(wanted)
    freqs: Dict[Key, float] = {}
    for v in read_vcf(path, wanted=wanted, per_allele=(af_field,)):
        raw = v.info.get(af_field)
        if raw is None:
            continue
        try:
            freqs[v.key] = float(raw)
        except ValueError:
            continue
    return freqs


def _from_ann(variant: Variant):
    """Parse SnpEff-style ANN (Allele|Effect|Impact|Gene|...) if the input VCF carries it."""
    raw = variant.info.get("ANN")
    if not raw:
        return "", []
    entries = [e.split("|") for e in raw.split(",")]
    entries = [e for e in entries if len(e) >= 4]
    if not entries:
        return "", []
    chosen = next((e for e in entries if e[0] == variant.alt), entries[0])
    effects = [t.lower() for t in chosen[1].split("&") if t]
    return chosen[3], effects


def annotate(
    variants: Iterable[Variant],
    clinvar: Dict[Key, ClinVarRecord],
    freqs: Dict[Key, float],
) -> List[Annotation]:
    out: List[Annotation] = []
    for v in variants:
        cv = clinvar.get(v.key)
        gene, effects = _from_ann(v)
        if not gene and cv:
            gene = cv.gene
        if not effects and cv:
            effects = cv.consequences
        out.append(Annotation(v, gene or "", effects, freqs.get(v.key), cv))
    return out
