from varlens.annotate import Annotation, ClinVarRecord, clinvar_category, review_stars
from varlens.rules import Thresholds, prioritise
from varlens.vcf import Variant


def _ann(af=None, cv=None, consequences=None):
    v = Variant("1", 1, "A", "T", genotype="0/1")
    return Annotation(v, "GENE", consequences or [], af, cv)


def _cv(category, stars):
    return ClinVarRecord(category.title(), category, "", stars, "GENE", "", [])


def test_clinvar_category_handles_conflicts_first():
    assert clinvar_category("Conflicting_classifications_of_pathogenicity") == "conflicting"
    assert clinvar_category("Pathogenic/Likely_pathogenic") == "pathogenic"
    assert clinvar_category("Benign/Likely_benign") == "benign"
    assert clinvar_category("Uncertain_significance") == "vus"


def test_review_stars():
    assert review_stars("practice_guideline") == 4
    assert review_stars("reviewed_by_expert_panel") == 3
    assert review_stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert review_stars("criteria_provided,_single_submitter") == 1
    assert review_stars("no_assertion_criteria_provided") == 0


def test_ba1_overrides_clinvar_pathogenic():
    assert prioritise(_ann(af=0.2, cv=_cv("pathogenic", 2))).tier == "low"


def test_pathogenic_two_stars_is_high_but_one_star_needs_review():
    assert prioritise(_ann(af=1e-5, cv=_cv("pathogenic", 2))).tier == "high"
    assert prioritise(_ann(af=1e-5, cv=_cv("pathogenic", 1))).tier == "review"


def test_pathogenic_but_common_is_flagged_for_review():
    res = prioritise(_ann(af=0.03, cv=_cv("pathogenic", 2)))
    assert res.tier == "review"
    assert any("Conflict" in e for e in res.evidence)


def test_novel_rare_lof_needs_review():
    assert prioritise(_ann(af=None, consequences=["frameshift_variant"])).tier == "review"


def test_common_missense_without_clinvar_is_low_and_unknown_is_vus():
    assert prioritise(_ann(af=0.02, consequences=["missense_variant"])).tier == "low"
    assert prioritise(_ann(af=0.0002, consequences=["missense_variant"])).tier == "vus"


def test_thresholds_are_configurable():
    strict = Thresholds(ba1=0.5, bs1=0.5, pm2=0.0001)
    assert prioritise(_ann(af=0.2, consequences=["missense_variant"]), strict).tier == "vus"
