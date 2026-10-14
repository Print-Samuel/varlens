"""The landing page re-implements rules.prioritise in JavaScript for its live demo.
This test runs that JS (via Node, if installed) and the Python rules over the same
grid of inputs and checks that they always pick the same tier, so the two can't drift."""
import itertools
import json
import re
import shutil
import subprocess
from pathlib import Path

from varlens.annotate import Annotation, ClinVarRecord
from varlens.rules import prioritise
from varlens.vcf import Variant

PAGE = Path(__file__).resolve().parent.parent / "docs" / "index.html"


def _rules_js() -> str:
    text = PAGE.read_text(encoding="utf-8")
    match = re.search(r"/\* RULES:START.*?\*/(.*?)/\* RULES:END \*/", text, re.S)
    assert match, "RULES markers not found in docs/index.html"
    return match.group(1)


def _grid():
    afs = [None, 0.0, 1e-6, 1e-4, 1.0001e-4, 0.0099, 0.01, 0.0499, 0.05, 0.3]
    cvs = [None] + [(cat, st) for cat in ("pathogenic", "benign", "vus", "conflicting", "other") for st in range(5)]
    for af, cv, lof in itertools.product(afs, cvs, (False, True)):
        yield af, cv, lof


def test_js_and_python_rules_agree_on_every_case():
    node = shutil.which("node")
    if node is None:
        return  # Node not installed: nothing to compare against
    cases = list(_grid())
    payload = [
        {"af": af, "lof": lof,
         "clinvar": None if cv is None else {"category": cv[0], "stars": cv[1], "label": cv[0].title()}}
        for af, cv, lof in cases
    ]
    script = _rules_js() + "\nvar cases = JSON.parse(require('fs').readFileSync(0, 'utf8'));\n" \
                           "console.log(JSON.stringify(cases.map(function (c) { return prioritise(c).tier; })));"
    out = subprocess.run([node, "-e", script], input=json.dumps(payload), capture_output=True, text=True, check=True)
    js_tiers = json.loads(out.stdout)

    for (af, cv, lof), js_tier in zip(cases, js_tiers):
        v = Variant("1", 1, "A", "T")
        record = None if cv is None else ClinVarRecord(cv[0].title(), cv[0], "", cv[1], "G", "", [])
        ann = Annotation(v, "G", ["stop_gained"] if lof else ["missense_variant"], af, record)
        assert prioritise(ann).tier == js_tier, (af, cv, lof)
    assert len(js_tiers) == len(cases) == 520
