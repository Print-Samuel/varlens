"""Blind benchmark of varlens' evidence rules against ClinVar-labelled variants.

Why "blind": varlens' normal mode lets ClinVar's own verdict decide the tier
(ClinVar pathogenic with 2+ review stars -> high priority). Scoring that mode on
ClinVar variants would be circular: the answer key is also an input. So this
benchmark hides the label. Each ClinVar variant is run through ``prioritise``
with ``clinvar=None``, using only evidence that does not carry the clinical
verdict:

* population allele frequency, taken from the AF_EXAC / AF_TGP / AF_ESP fields
  that ClinVar itself ships in its VCF (older and smaller cohorts than gnomAD;
  a variant missing from them is treated as "not observed", i.e. rare);
* molecular consequence, from ClinVar's MC field (loss of function or not).

The ClinVar label (pathogenic vs benign, with at least ``--min-stars`` review
stars) is then used only afterwards, as the answer key. Because the "high"
tier requires a ClinVar assertion, it cannot occur here: the question this
benchmark answers is how much of the pathogenic/benign signal the
frequency and consequence rules recover *on their own*.

Reported rates are conditional on the true class (so they do not depend on how
many pathogenic vs benign variants ClinVar happens to contain) and come with
95% Wilson score intervals. Positive predictive value is also given, with the
prevalence it depends on.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from urllib.parse import unquote

from . import __version__
from .annotate import Annotation, _consequences_from_mc, _gene_from_geneinfo, clinvar_category, review_stars
from .rules import TIER_ORDER, Thresholds, is_lof, prioritise
from .vcf import parse_region, read_vcf

AF_FIELDS = ("AF_EXAC", "AF_TGP", "AF_ESP")
FLAG_TIERS = ("high", "review")  # "needs a human to look at it"
LABELS = ("pathogenic", "benign")

NOTES = [
    "Blind mode: the ClinVar label is never given to the rules; it is only the answer key.",
    "Frequencies come from ClinVar's AF_EXAC / AF_TGP / AF_ESP fields (older, smaller cohorts than gnomAD). "
    "A variant missing from them is treated as not observed (rare).",
    "The 'high' tier needs a ClinVar assertion, so it cannot occur in blind mode.",
    "Rates are conditional on the true class, with 95% Wilson score intervals; "
    "positive predictive value depends on the pathogenic/benign mix (prevalence is reported).",
    "Labels are ClinVar's own and are not ground truth: they can be revised, and reflect which genes and variants get submitted.",
]


def wilson_interval(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """95% Wilson score interval for a binomial proportion k/n."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _proportion(k: int, n: int) -> dict:
    lo, hi = wilson_interval(k, n)
    return {"k": k, "n": n, "rate": (k / n if n else None), "ci": [lo, hi]}


def max_population_af(info: Dict[str, str]) -> Optional[float]:
    """Highest allele frequency among ClinVar's AF_* fields, or None if none are present."""
    best: Optional[float] = None
    for name in AF_FIELDS:
        raw = info.get(name)
        if raw is None:
            continue
        try:
            value = float(raw.split(",")[0])
        except ValueError:
            continue
        best = value if best is None else max(best, value)
    return best


def label_for(significance: str) -> Optional[str]:
    category = clinvar_category(significance)
    return category if category in LABELS else None


def run_benchmark(
    clinvar_path: Union[str, Path],
    min_stars: int = 1,
    region=None,
    thresholds: Thresholds = Thresholds(),
    max_examples: int = 15,
) -> dict:
    tiers = {label: Counter() for label in LABELS}
    n = Counter()
    af_seen = Counter()
    lof_seen = Counter()
    dismissed: List[dict] = []

    for v in read_vcf(clinvar_path, region=region, per_allele=()):
        significance = unquote(v.info.get("CLNSIG", ""))
        label = label_for(significance)
        if label is None:
            continue  # VUS, conflicting and other categories are not a clean answer key
        stars = review_stars(unquote(v.info.get("CLNREVSTAT", "")))
        if stars < min_stars:
            continue

        af = max_population_af(v.info)
        consequences = _consequences_from_mc(v.info.get("MC", ""))
        gene = _gene_from_geneinfo(v.info.get("GENEINFO", ""))
        # clinvar=None is what keeps this benchmark blind.
        tier = prioritise(Annotation(v, gene, consequences, af, None), thresholds).tier

        n[label] += 1
        tiers[label][tier] += 1
        if af is not None:
            af_seen[label] += 1
        if is_lof(consequences):
            lof_seen[label] += 1
        if label == "pathogenic" and tier == "low":
            dismissed.append({
                "locus": v.locus,
                "gene": gene or "?",
                "af": af,
                "significance": significance.replace("_", " "),
                "condition": unquote(v.info.get("CLNDN", "")).replace("_", " ").split("|")[0],
            })

    dismissed.sort(key=lambda d: -(d["af"] or 0))
    p, b = n["pathogenic"], n["benign"]
    flag_p = sum(tiers["pathogenic"][t] for t in FLAG_TIERS)
    flag_b = sum(tiers["benign"][t] for t in FLAG_TIERS)

    return {
        "varlens_version": __version__,
        "date": _dt.date.today().isoformat(),
        "clinvar_file": Path(clinvar_path).name,
        "region": None,
        "min_stars": min_stars,
        "mode": "blind",
        "thresholds": {"ba1": thresholds.ba1, "bs1": thresholds.bs1, "pm2": thresholds.pm2},
        "n": {"pathogenic": p, "benign": b},
        "prevalence_pathogenic": (p / (p + b) if p + b else None),
        "tiers": {label: {t: tiers[label][t] for t in TIER_ORDER} for label in LABELS},
        "af_coverage": {label: _proportion(af_seen[label], n[label]) for label in LABELS},
        "lof_fraction": {label: _proportion(lof_seen[label], n[label]) for label in LABELS},
        "metrics": {
            "sensitivity_flag": _proportion(flag_p, p),
            "false_alarm_rate": _proportion(flag_b, b),
            "pathogenic_dismissed": _proportion(tiers["pathogenic"]["low"], p),
            "benign_recognised": _proportion(tiers["benign"]["low"], b),
            "ppv_flag": _proportion(flag_p, flag_p + flag_b),
        },
        "dismissed_pathogenic": {"count": len(dismissed), "examples": dismissed[:max_examples]},
        "notes": NOTES,
    }


def _pct(x: Optional[float]) -> str:
    if x is None:
        return "n/a"
    return f"{100 * x:.3f}%" if 0 < x < 0.001 else f"{100 * x:.1f}%"


def _fmt(m: dict) -> str:
    lo, hi = m["ci"]
    return f"{_pct(m['rate'])} ({m['k']}/{m['n']}, 95% CI {_pct(lo)} to {_pct(hi)})"


def format_summary(r: dict) -> str:
    n, m = r["n"], r["metrics"]
    lines = [
        f"varlens {r['varlens_version']} blind benchmark on {r['clinvar_file']}"
        + (f" (region {r['region']})" if r["region"] else ""),
        f"labelled variants (>= {r['min_stars']} review star(s)): "
        f"{n['pathogenic']} pathogenic, {n['benign']} benign",
        "",
        "tier distribution (count)        high  review     vus     low",
    ]
    for label in LABELS:
        t = r["tiers"][label]
        lines.append(f"  {label:<10}                 {t['high']:>7}{t['review']:>8}{t['vus']:>8}{t['low']:>8}")
    lines += [
        "",
        f"pathogenic flagged for review     {_fmt(m['sensitivity_flag'])}",
        f"benign flagged for review         {_fmt(m['false_alarm_rate'])}",
        f"pathogenic dismissed as low       {_fmt(m['pathogenic_dismissed'])}",
        f"benign recognised as low          {_fmt(m['benign_recognised'])}",
        f"PPV of the review flag            {_fmt(m['ppv_flag'])}  [pathogenic prevalence {_pct(r['prevalence_pathogenic'])}]",
        "",
        f"allele frequency available: pathogenic {_pct(r['af_coverage']['pathogenic']['rate'])}, "
        f"benign {_pct(r['af_coverage']['benign']['rate'])}",
        f"loss-of-function consequence: pathogenic {_pct(r['lof_fraction']['pathogenic']['rate'])}, "
        f"benign {_pct(r['lof_fraction']['benign']['rate'])}",
    ]
    d = r["dismissed_pathogenic"]
    if d["examples"]:
        lines += ["", f"pathogenic variants dismissed as low priority ({d['count']} total); most common:"]
        for e in d["examples"][:5]:
            lines.append(f"  {e['gene']:<10} {e['locus']:<28} AF {e['af']:.4g}  {e['condition']}")
    return "\n".join(lines)


def ppv_at_prevalence(sensitivity: float, false_alarm: float, prevalence: float) -> float:
    """PPV of the review flag if the pathogenic share were `prevalence` instead of ClinVar's."""
    tp = sensitivity * prevalence
    fp = false_alarm * (1 - prevalence)
    return tp / (tp + fp) if (tp + fp) else 0.0


def interpretation(r: dict) -> List[str]:
    """Plain-English reading of the numbers, built from the results so it never drifts from them."""
    n, t, m = r["n"], r["tiers"]["pathogenic"], r["metrics"]
    lof_p, lof_b = r["lof_fraction"]["pathogenic"], r["lof_fraction"]["benign"]
    af_p, af_b = r["af_coverage"]["pathogenic"], r["af_coverage"]["benign"]
    lines = [
        f"The review flag mostly tracks loss of function: {_pct(lof_p['rate'])} of pathogenic variants "
        f"have a loss-of-function consequence, against {_pct(lof_b['rate'])} of benign ones. ClinVar "
        "curators also use consequence as evidence, so part of this agreement is by construction.",
    ]
    if n["pathogenic"]:
        lines.append(
            f"{t['vus']} pathogenic variants ({_pct(t['vus'] / n['pathogenic'])}) stay in the uncertain tier, "
            f"and {n['pathogenic'] - lof_p['k']} pathogenic variants have no loss-of-function consequence. "
            "Frequency and consequence alone cannot separate these from benign variants; this is where "
            "missense-aware evidence (predictors, a trained model) would have to add value."
        )
    sens, fa = m["sensitivity_flag"]["rate"], m["false_alarm_rate"]["rate"]
    if sens is not None and fa is not None:
        lines.append(
            f"The positive predictive value depends on prevalence. ClinVar's answer key is "
            f"{_pct(r['prevalence_pathogenic'])} pathogenic; among the rare variants in a real sample, "
            f"pathogenic ones are far scarcer. With the same sensitivity and false-alarm rate, a prevalence "
            f"of 0.1% would give a PPV of about {_pct(ppv_at_prevalence(sens, fa, 0.001))}."
        )
    lines.append(
        f"Frequency data exists for {_pct(af_p['rate'])} of pathogenic and {_pct(af_b['rate'])} of benign "
        "variants, so the frequency rules rarely get to act on pathogenic variants, and a missing frequency "
        "is treated as rare."
    )
    return lines


def format_markdown(r: dict) -> str:
    n, m = r["n"], r["metrics"]
    rows = "\n".join(
        f"| {label} | {r['tiers'][label]['high']} | {r['tiers'][label]['review']} | "
        f"{r['tiers'][label]['vus']} | {r['tiers'][label]['low']} |"
        for label in LABELS
    )
    notes = "\n".join(f"- {x}" for x in r["notes"])
    reading = "\n".join(f"- {x}" for x in interpretation(r))
    return f"""# varlens blind benchmark

Run on `{r['clinvar_file']}` with varlens {r['varlens_version']} on {r['date']}.
Labelled variants with at least {r['min_stars']} review star(s): {n['pathogenic']} pathogenic, {n['benign']} benign.

| ClinVar label | high | review | uncertain | low |
|---|---|---|---|---|
{rows}

| Question | Result |
|---|---|
| Pathogenic variants flagged for review | {_fmt(m['sensitivity_flag'])} |
| Benign variants flagged for review | {_fmt(m['false_alarm_rate'])} |
| Pathogenic variants dismissed as low priority | {_fmt(m['pathogenic_dismissed'])} |
| Benign variants recognised as low priority | {_fmt(m['benign_recognised'])} |
| Positive predictive value of the review flag | {_fmt(m['ppv_flag'])} |

## How to read these numbers

{reading}

## Method and caveats

{notes}
"""


_EMBED_RE = re.compile(r'(<script id="benchmark-data" type="application/json">)(.*?)(</script>)', re.S)


def embed_results(html_path: Union[str, Path], results: dict) -> None:
    """Write ``results`` into the landing page's <script id="benchmark-data"> block."""
    text = Path(html_path).read_text(encoding="utf-8")
    payload = json.dumps(results, ensure_ascii=False).replace("</", "<\\/")
    new_text, count = _EMBED_RE.subn(lambda m: m.group(1) + payload + m.group(3), text, count=1)
    if count != 1:
        raise ValueError(f'no <script id="benchmark-data" type="application/json"> block found in {html_path}')
    Path(html_path).write_text(new_text, encoding="utf-8")


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        prog="varlens-benchmark",
        description="Blind benchmark of varlens' frequency and consequence rules against ClinVar-labelled variants.",
    )
    p.add_argument("clinvar", help="ClinVar VCF (.vcf or .vcf.gz)")
    p.add_argument("--min-stars", type=int, default=1, help="minimum ClinVar review stars for the answer key (default 1)")
    p.add_argument("--region", help="restrict to a chromosome or window, e.g. '21' or '21:1000000-30000000'")
    p.add_argument("--ba1", type=float, default=0.05)
    p.add_argument("--bs1", type=float, default=0.01)
    p.add_argument("--pm2", type=float, default=0.0001)
    p.add_argument("--out", default="benchmark_results.json", help="JSON results path")
    p.add_argument("--markdown", help="also write a Markdown summary (e.g. BENCHMARK.md)")
    p.add_argument("--embed", help="landing page HTML to update with the results (e.g. docs/index.html)")
    args = p.parse_args(argv)

    if not Path(args.clinvar).exists():
        print(f"error: ClinVar file not found: {args.clinvar}", file=sys.stderr)
        return 2
    region = None
    if args.region:
        try:
            region = parse_region(args.region)
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
    if args.embed and not Path(args.embed).exists():
        print(f"error: landing page not found: {args.embed}", file=sys.stderr)
        return 2

    print(f"scanning {args.clinvar} (this reads the whole file once)...", file=sys.stderr)
    results = run_benchmark(args.clinvar, args.min_stars, region, Thresholds(args.ba1, args.bs1, args.pm2))
    results["region"] = args.region

    if results["n"]["pathogenic"] == 0 or results["n"]["benign"] == 0:
        print("error: need both pathogenic and benign labelled variants; try a lower --min-stars or a wider --region",
              file=sys.stderr)
        return 2
    if results["af_coverage"]["pathogenic"]["k"] + results["af_coverage"]["benign"]["k"] == 0:
        print("warning: no AF_EXAC/AF_TGP/AF_ESP values found; frequency rules had nothing to work with "
              "(is this the full ClinVar VCF?)", file=sys.stderr)

    Path(args.out).write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    if args.markdown:
        Path(args.markdown).write_text(format_markdown(results), encoding="utf-8")
    if args.embed:
        embed_results(args.embed, results)
    print(format_summary(results))
    print(f"\nresults -> {args.out}" + (f", {args.markdown}" if args.markdown else "")
          + (f", embedded in {args.embed}" if args.embed else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
