import os
import tempfile

from varlens.vcf import normalize_chrom, parse_info, read_vcf

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def _write(body: str) -> str:
    handle = tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False)
    handle.write(HEADER + body)
    handle.close()
    return handle.name


def test_normalize_chrom():
    assert normalize_chrom("chr7") == "7"
    assert normalize_chrom("7") == "7"
    assert normalize_chrom("chrM") == "MT"


def test_parse_info_flags_and_values():
    assert parse_info("AF=0.1;DB") == {"AF": "0.1", "DB": "true"}
    assert parse_info(".") == {}


def test_multiallelic_split_and_genotype():
    path = _write("1\t100\t.\tG\tA,C\t.\t.\tAF=0.1,0.2\tGT\t1/2\n")
    try:
        a, c = list(read_vcf(path))
    finally:
        os.unlink(path)
    assert (a.alt, a.info["AF"], a.zygosity) == ("A", "0.1", "heterozygous")
    assert (c.alt, c.info["AF"], c.zygosity) == ("C", "0.2", "heterozygous")


def test_homozygous_and_chr_prefix():
    path = _write("chr2\t50\t.\tT\tTA\t.\t.\t.\tGT\t1|1\n")
    try:
        (v,) = list(read_vcf(path))
    finally:
        os.unlink(path)
    assert v.chrom == "2"
    assert v.zygosity == "homozygous"


def test_wanted_filter_skips_other_variants():
    path = _write("1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n1\t2\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")
    try:
        found = list(read_vcf(path, wanted={("1", 2, "A", "G")}))
    finally:
        os.unlink(path)
    assert [v.pos for v in found] == [2]
