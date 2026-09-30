import os
import tempfile
from pathlib import Path

from varlens.cli import main

DATA = Path(__file__).resolve().parent.parent / "data"


def test_end_to_end_on_synthetic_example():
    with tempfile.TemporaryDirectory() as tmp:
        html_path = os.path.join(tmp, "report.html")
        tsv_path = os.path.join(tmp, "report.tsv")
        code = main([
            str(DATA / "example_patient.vcf"),
            "--clinvar", str(DATA / "example_clinvar.vcf"),
            "--freq", str(DATA / "example_freq.vcf"),
            "--out", html_path,
            "--tsv", tsv_path,
        ])
        assert code == 0
        html = Path(html_path).read_text(encoding="utf-8")
        rows = Path(tsv_path).read_text(encoding="utf-8").splitlines()

    assert "Educational use only" in html
    assert len(rows) == 1 + 9  # header + 8 sites, one of them multi-allelic (2 alleles)
    tiers = [r.split("\t")[0] for r in rows[1:]]
    assert tiers[0] == "high"  # sorted by priority
    assert tiers == sorted(tiers, key=["high", "review", "vus", "low"].index)


def test_missing_input_returns_error_code():
    assert main(["does_not_exist.vcf"]) == 2
