"""Render annotated, prioritised variants as a self-contained HTML report or a TSV."""
from __future__ import annotations

import csv
import datetime as _dt
from html import escape
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Union

from . import __version__
from .annotate import Annotation
from .rules import TIER_LABEL, TIER_ORDER, Result

Row = Tuple[Annotation, Result]

_CSS = """
:root{--bg:#fff;--fg:#1c1f24;--muted:#5b6470;--line:#d8dde3;--card:#f6f8fa}
@media (prefers-color-scheme:dark){:root{--bg:#14171b;--fg:#e6e9ed;--muted:#9aa4b1;--line:#2c333b;--card:#1b2026}}
body{font:15px/1.5 system-ui,sans-serif;background:var(--bg);color:var(--fg);margin:0;padding:24px;max-width:1200px;margin-inline:auto}
h1{font-size:22px;margin:0 0 4px}.muted{color:var(--muted)}
.warn{background:var(--card);border:1px solid var(--line);border-left:4px solid #c98a00;padding:10px 14px;border-radius:6px;margin:16px 0}
.summary{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}
.summary div{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 16px;min-width:130px}
.summary b{font-size:22px;display:block}
.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%}
th,td{border-bottom:1px solid var(--line);padding:8px 10px;text-align:left;vertical-align:top}
th{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
.tag{display:inline-block;padding:2px 8px;border-radius:99px;font-size:12px;font-weight:600;white-space:nowrap}
.high{background:#f8d7da;color:#842029}.review{background:#fff0c2;color:#664d03}
.vus{background:#e2e3e5;color:#41464b}.low{background:#d1e7dd;color:#0f5132}
ul{margin:0;padding-left:18px}code{font-family:ui-monospace,monospace;font-size:13px}
"""


def sort_rows(rows: Sequence[Row]) -> List[Row]:
    def key(row: Row):
        ann, res = row
        v = ann.variant
        chrom = v.chrom
        chrom_key = (0, int(chrom)) if chrom.isdigit() else (1, chrom)
        return (TIER_ORDER[res.tier], chrom_key, v.pos)

    return sorted(rows, key=key)


def _clinvar_cell(ann: Annotation) -> str:
    cv = ann.clinvar
    if not cv:
        return '<span class="muted">not in ClinVar</span>'
    stars = "&#9733;" * cv.stars if cv.stars else "no stars"
    return f"{escape(cv.significance.replace('_', ' '))}<br><span class='muted'>{stars}</span>"


def write_html(
    rows: Sequence[Row],
    path: Union[str, Path],
    meta: Dict[str, str],
) -> None:
    rows = sort_rows(rows)
    counts = {t: 0 for t in TIER_ORDER}
    for _, res in rows:
        counts[res.tier] += 1

    summary = "".join(
        f"<div><b>{counts[t]}</b>{escape(TIER_LABEL[t])}</div>" for t in TIER_ORDER
    )
    body = []
    for ann, res in rows:
        v = ann.variant
        af = "-" if ann.allele_freq is None else f"{ann.allele_freq:.3g}"
        effects = ", ".join(ann.consequences) or "-"
        disease = escape(ann.clinvar.disease) if ann.clinvar and ann.clinvar.disease else "-"
        evidence = "".join(f"<li>{escape(e)}</li>" for e in res.evidence)
        body.append(
            "<tr>"
            f'<td><span class="tag {res.tier}">{escape(res.label)}</span></td>'
            f"<td><code>{escape(v.locus)}</code></td>"
            f"<td>{escape(ann.gene or '-')}</td>"
            f"<td>{escape(effects)}</td>"
            f"<td>{escape(v.zygosity)}</td>"
            f"<td>{af}</td>"
            f"<td>{_clinvar_cell(ann)}</td>"
            f"<td>{disease}</td>"
            f"<td><ul>{evidence}</ul></td>"
            "</tr>"
        )

    meta_html = " &middot; ".join(f"{escape(k)}: {escape(v)}" for k, v in meta.items())
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>varlens report</title><style>{_CSS}</style></head>
<body>
<h1>varlens variant report</h1>
<div class="muted">{meta_html} &middot; varlens {__version__} &middot; {_dt.date.today().isoformat()}</div>
<div class="warn"><b>Educational use only.</b> This is a rule-based prioritisation, not a clinical
classification. Do not use it for diagnosis or medical decisions.</div>
<div class="summary">{summary}</div>
<div class="scroll"><table>
<thead><tr><th>Tier</th><th>Variant</th><th>Gene</th><th>Consequence</th><th>Zygosity</th>
<th>Allele freq</th><th>ClinVar</th><th>Condition</th><th>Evidence</th></tr></thead>
<tbody>{''.join(body)}</tbody></table></div>
</body></html>"""
    Path(path).write_text(html, encoding="utf-8")


def write_tsv(rows: Sequence[Row], path: Union[str, Path]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["tier", "chrom", "pos", "ref", "alt", "gene", "consequence", "zygosity",
             "allele_freq", "clinvar_significance", "clinvar_stars", "condition", "evidence"]
        )
        for ann, res in sort_rows(rows):
            v, cv = ann.variant, ann.clinvar
            writer.writerow([
                res.tier, v.chrom, v.pos, v.ref, v.alt, ann.gene, ";".join(ann.consequences),
                v.zygosity, "" if ann.allele_freq is None else ann.allele_freq,
                cv.significance if cv else "", cv.stars if cv else "",
                cv.disease if cv else "", " | ".join(res.evidence),
            ])
