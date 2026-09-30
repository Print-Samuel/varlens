import os
import tempfile

import pytest

from varlens.vcf import parse_region, read_vcf, region_contains

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _write(body: str) -> str:
    handle = tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False)
    handle.write(HEADER + body)
    handle.close()
    return handle.name


def test_parse_region_whole_chromosome():
    assert parse_region("21") == ("21", None, None)
    assert parse_region("chr21") == ("21", None, None)  # chr prefix normalised


def test_parse_region_bounded():
    assert parse_region("21:100-200") == ("21", 100, 200)


def test_parse_region_open_ended():
    assert parse_region("21:100-") == ("21", 100, None)
    assert parse_region("21:-200") == ("21", None, 200)


def test_parse_region_rejects_start_after_end():
    with pytest.raises(ValueError):
        parse_region("21:200-100")


def test_parse_region_rejects_missing_dash():
    with pytest.raises(ValueError):
        parse_region("21:100")


def test_region_contains():
    region = ("21", 100, 200)
    assert region_contains(region, "21", 150)
    assert not region_contains(region, "21", 50)
    assert not region_contains(region, "22", 150)


def test_read_vcf_filters_by_region():
    path = _write(
        "21\t50\t.\tA\tT\t.\t.\t.\n"
        "21\t150\t.\tA\tT\t.\t.\t.\n"
        "21\t250\t.\tA\tT\t.\t.\t.\n"
        "22\t150\t.\tA\tT\t.\t.\t.\n"
    )
    try:
        found = list(read_vcf(path, region=parse_region("21:100-200")))
    finally:
        os.unlink(path)
    assert [v.pos for v in found] == [150]
