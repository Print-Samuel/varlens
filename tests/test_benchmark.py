import json
import os
import tempfile
from pathlib import Path

import pytest

from varlens.benchmark import (
    embed_results,
    format_markdown,
    interpretation,
    ppv_at_prevalence,
    format_summary,
    main,
    max_population_af,
    run_benchmark,
    wilson_interval,
)
from varlens.vcf import parse_region

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
TWO = "criteria_provided,_multiple_submitters,_no_conflicts"
ONE = "criteria_provided,_single_submitter"
NONE = "no_assertion_criteria_provided"


def _row(chrom, pos, ref, alt, sig, rev, mc, extra=""):
    info = f"CLNSIG={sig};CLNREVSTAT={rev};GENEINFO=GENE{pos}:1;CLNDN=Demo_condition;MC={mc}"
    if extra:
        info += ";" + extra
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\n"


NONSENSE = "SO:0001587|nonsense"
MISSENSE = "SO:0001583|missense_variant"
FRAMESHIFT = "SO:0001589|frameshift_variant"


def _clinvar_file() -> str:
    body = "".join([
        _row("1", 100, "G", "A", "Pathogenic", TWO, NONSENSE),
        _row("1", 200, "C", "T", "Pathogenic", TWO, MISSENSE),
        _row("1", 300, "A", "G", "Pathogenic", ONE, MISSENSE, "AF_EXAC=0.06"),
        _row("1", 400, "T", "C", "Benign", ONE, MISSENSE, "AF_TGP=0.2"),
        _row("1", 500, "G", "T", "Benign", ONE, MISSENSE),
        _row("1", 600, "C", "A", "Benign", ONE, NONSENSE),
        _row("1", 700, "A", "T", "Uncertain_significance", ONE, MISSENSE),
        _row("1", 800, "A", "T", "Conflicting_classifications_of_pathogenicity", ONE, MISSENSE),
        _row("1", 900, "A", "T", "Pathogenic", NONE, NONSENSE),
        _row("2", 100, "A", "AT", "Pathogenic", ONE, FRAMESHIFT, "AF_ESP=0.00001"),
    ])
    handle = tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False)
    handle.write(HEADER + body)
    handle.close()
    return handle.name


def test_wilson_interval_known_values():
    lo, hi = wilson_interval(0, 10)
    assert lo == 0.0
    assert hi == pytest.approx(0.2775, rel=1e-3)
    lo, hi = wilson_interval(5, 10)
    assert lo == pytest.approx(0.2366, rel=1e-3)
    assert hi == pytest.approx(0.7634, rel=1e-3)
    assert wilson_interval(0, 0) == (0.0, 1.0)


def test_max_population_af_takes_highest_and_tolerates_missing_or_bad_values():
    assert max_population_af({}) is None
    assert max_population_af({"AF_EXAC": "0.01", "AF_TGP": "0.2", "AF_ESP": "0.003"}) == pytest.approx(0.2)
    assert max_population_af({"AF_EXAC": "not-a-number", "AF_ESP": "0.5"}) == pytest.approx(0.5)


def test_blind_benchmark_counts_tiers_and_metrics():
    path = _clinvar_file()
    try:
        r = run_benchmark(path, min_stars=1)
    finally:
        os.unlink(path)

    assert r["n"] == {"pathogenic": 4, "benign": 3}  # VUS, conflicting and 0-star are excluded
    assert r["tiers"]["pathogenic"] == {"high": 0, "review": 2, "vus": 1, "low": 1}
    assert r["tiers"]["benign"] == {"high": 0, "review": 1, "vus": 1, "low": 1}
    m = r["metrics"]
    assert (m["sensitivity_flag"]["k"], m["sensitivity_flag"]["n"]) == (2, 4)
    assert (m["false_alarm_rate"]["k"], m["false_alarm_rate"]["n"]) == (1, 3)
    assert (m["pathogenic_dismissed"]["k"], m["pathogenic_dismissed"]["n"]) == (1, 4)
    assert (m["benign_recognised"]["k"], m["benign_recognised"]["n"]) == (1, 3)
    assert (m["ppv_flag"]["k"], m["ppv_flag"]["n"]) == (2, 3)
    assert r["prevalence_pathogenic"] == pytest.approx(4 / 7)
    assert (r["af_coverage"]["pathogenic"]["k"], r["af_coverage"]["benign"]["k"]) == (2, 1)
    assert (r["lof_fraction"]["pathogenic"]["k"], r["lof_fraction"]["benign"]["k"]) == (2, 1)


def test_benchmark_is_blind_a_two_star_pathogenic_missense_is_not_high():
    """If the label leaked into the rules, variant #2 (pathogenic, 2 stars) would be 'high'."""
    path = _clinvar_file()
    try:
        r = run_benchmark(path, min_stars=1)
    finally:
        os.unlink(path)
    assert r["tiers"]["pathogenic"]["high"] == 0
    assert r["tiers"]["pathogenic"]["vus"] == 1


def test_dismissed_pathogenic_examples_carry_gene_and_frequency():
    path = _clinvar_file()
    try:
        r = run_benchmark(path, min_stars=1)
    finally:
        os.unlink(path)
    d = r["dismissed_pathogenic"]
    assert d["count"] == 1
    assert d["examples"][0]["gene"] == "GENE300"
    assert d["examples"][0]["af"] == pytest.approx(0.06)


def test_min_stars_and_region_filter_the_answer_key():
    path = _clinvar_file()
    try:
        strict = run_benchmark(path, min_stars=2)
        chr2 = run_benchmark(path, min_stars=1, region=parse_region("2"))
    finally:
        os.unlink(path)
    assert strict["n"] == {"pathogenic": 2, "benign": 0}
    assert strict["metrics"]["benign_recognised"]["rate"] is None
    assert chr2["n"] == {"pathogenic": 1, "benign": 0}


def test_summary_and_markdown_mention_the_key_numbers():
    path = _clinvar_file()
    try:
        r = run_benchmark(path, min_stars=1)
    finally:
        os.unlink(path)
    text = format_summary(r)
    assert "4 pathogenic, 3 benign" in text
    md = format_markdown(r)
    assert "Blind mode" in md and "| pathogenic |" in md


def test_embed_results_replaces_block_and_is_script_safe():
    results = {"note": "contains </script> inside", "n": {"pathogenic": 1, "benign": 1}}
    page = '<html><body><script id="benchmark-data" type="application/json">null</script><p>x</p></body></html>'
    with tempfile.TemporaryDirectory() as tmp:
        html_path = Path(tmp) / "index.html"
        html_path.write_text(page, encoding="utf-8")
        embed_results(html_path, results)
        out = html_path.read_text(encoding="utf-8")
    assert out.count("</script>") == 1  # the payload's own '</script>' was escaped
    start = out.index('type="application/json">') + len('type="application/json">')
    payload = out[start:out.index("</script>")]
    assert json.loads(payload) == results


def test_embed_results_errors_when_marker_missing():
    with tempfile.TemporaryDirectory() as tmp:
        html_path = Path(tmp) / "index.html"
        html_path.write_text("<html></html>", encoding="utf-8")
        with pytest.raises(ValueError):
            embed_results(html_path, {})


def test_cli_writes_json_markdown_and_embeds():
    path = _clinvar_file()
    with tempfile.TemporaryDirectory() as tmp:
        page = Path(tmp) / "index.html"
        page.write_text('<script id="benchmark-data" type="application/json">null</script>', encoding="utf-8")
        out_json = Path(tmp) / "r.json"
        out_md = Path(tmp) / "B.md"
        try:
            code = main([path, "--out", str(out_json), "--markdown", str(out_md), "--embed", str(page)])
        finally:
            os.unlink(path)
        assert code == 0
        assert json.loads(out_json.read_text(encoding="utf-8"))["n"]["pathogenic"] == 4
        assert "varlens blind benchmark" in out_md.read_text(encoding="utf-8")
        assert '"pathogenic": 4' in page.read_text(encoding="utf-8")


def test_cli_errors_cleanly():
    assert main(["does_not_exist.vcf"]) == 2
    path = _clinvar_file()
    try:
        assert main([path, "--region", "1:500-100"]) == 2  # invalid region
        assert main([path, "--min-stars", "4", "--out", os.devnull]) == 2  # no labelled variants left
    finally:
        os.unlink(path)


def test_ppv_at_prevalence_matches_bayes_and_falls_with_prevalence():
    assert ppv_at_prevalence(0.8, 0.0, 0.01) == pytest.approx(1.0)
    assert ppv_at_prevalence(0.5, 0.5, 0.5) == pytest.approx(0.5)
    assert ppv_at_prevalence(0.784, 0.00204, 0.001) == pytest.approx(0.2775, abs=1e-3)
    assert ppv_at_prevalence(0.784, 0.00204, 0.001) < ppv_at_prevalence(0.784, 0.00204, 0.2)


def test_interpretation_and_markdown_explain_how_to_read_the_numbers():
    path = _clinvar_file()
    try:
        r = run_benchmark(path, min_stars=1)
    finally:
        os.unlink(path)
    text = " ".join(interpretation(r))
    assert "loss of function" in text and "prevalence" in text and "Frequency data exists" in text
    assert "## How to read these numbers" in format_markdown(r)
