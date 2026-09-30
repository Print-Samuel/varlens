"""Minimal, dependency-free VCF reader.

Only what varlens needs: sites, INFO fields and the first sample's genotype.
Multi-allelic records are split into one ``Variant`` per ALT allele.
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Set, Tuple, Union

Key = Tuple[str, int, str, str]  # (chrom, pos, ref, alt)


def normalize_chrom(chrom: str) -> str:
    """Strip a leading 'chr' and unify the mitochondrial name ('M' -> 'MT')."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return "MT" if name == "M" else name


# (chrom, start, end) -- start/end are inclusive, 1-based (VCF's own convention).
# Either bound may be None, meaning "no limit on that side".
Region = Tuple[str, Optional[int], Optional[int]]


def parse_region(spec: str) -> Region:
    """Parse '21', '21:1000000-2000000', '21:1000000-' or '21:-2000000' into a Region.

    Raises ValueError on anything else (e.g. missing '-' in a ranged spec, or
    start > end), so callers can turn it into a clean CLI error message.
    """
    if ":" not in spec:
        return normalize_chrom(spec), None, None
    chrom, _, rng = spec.partition(":")
    chrom = normalize_chrom(chrom)
    if "-" not in rng:
        raise ValueError(f"invalid region {spec!r}: expected 'chrom' or 'chrom:start-end'")
    start_s, _, end_s = rng.partition("-")
    try:
        start = int(start_s) if start_s else None
        end = int(end_s) if end_s else None
    except ValueError:
        raise ValueError(f"invalid region {spec!r}: start/end must be integers")
    if start is not None and end is not None and start > end:
        raise ValueError(f"invalid region {spec!r}: start ({start}) > end ({end})")
    return chrom, start, end


def region_contains(region: Region, chrom: str, pos: int) -> bool:
    r_chrom, start, end = region
    if chrom != r_chrom:
        return False
    if start is not None and pos < start:
        return False
    if end is not None and pos > end:
        return False
    return True


def parse_info(raw: str) -> Dict[str, str]:
    """Parse a VCF INFO column into a dict. Flags get the value 'true'."""
    info: Dict[str, str] = {}
    if raw in (".", ""):
        return info
    for item in raw.split(";"):
        if not item:
            continue
        name, sep, value = item.partition("=")
        info[name] = value if sep else "true"
    return info


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    id: str = "."
    qual: str = "."
    filter: str = "."
    info: Dict[str, str] = field(default_factory=dict)
    genotype: Optional[str] = None  # relative to THIS alt: "0"=ref, "1"=this alt, "2"=other alt

    @property
    def key(self) -> Key:
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def locus(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"

    @property
    def zygosity(self) -> str:
        if not self.genotype:
            return "unknown"
        alleles = re.split(r"[/|]", self.genotype)
        if any(a == "." for a in alleles):
            return "unknown"
        n_alt = alleles.count("1")
        if n_alt == 0:
            return "ref"
        if len(alleles) == 1:
            return "hemizygous"
        if n_alt == len(alleles):
            return "homozygous"
        return "heterozygous"


def _open_text(path: Union[str, Path]):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def _remap_genotype(raw_gt: str, alt_index: int) -> str:
    """Rewrite a GT string relative to one ALT allele (1-based ``alt_index``)."""
    sep = "|" if "|" in raw_gt else "/"
    out = []
    for allele in re.split(r"[/|]", raw_gt):
        if allele == ".":
            out.append(".")
        elif allele == "0":
            out.append("0")
        elif allele == str(alt_index):
            out.append("1")
        else:
            out.append("2")
    return sep.join(out)


def read_vcf(
    path: Union[str, Path],
    wanted: Optional[Set[Key]] = None,
    per_allele: Iterable[str] = ("AF", "AC"),
    region: Optional[Region] = None,
) -> Iterator[Variant]:
    """Yield one ``Variant`` per ALT allele.

    ``wanted``: if given, only variants whose key is in the set are yielded
    (lets us stream huge ClinVar/gnomAD files while keeping memory tiny).
    ``per_allele``: INFO keys with one comma-separated value per ALT (Number=A).
    ``region``: if given, skip any line outside it before parsing the rest of
    the line -- lets you scan one chromosome or a coordinate window of a
    whole-genome VCF without reading (or allocating a ``Variant`` for) the
    millions of lines you don't care about.
    """
    per_allele = tuple(per_allele)
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                raise ValueError(f"Malformed VCF line (fewer than 5 columns): {line[:80]!r}")
            chrom = normalize_chrom(cols[0])
            pos = int(cols[1])
            if region is not None and not region_contains(region, chrom, pos):
                continue
            vid, ref, alts = cols[2], cols[3], cols[4].split(",")
            qual = cols[5] if len(cols) > 5 else "."
            flt = cols[6] if len(cols) > 6 else "."
            info = parse_info(cols[7]) if len(cols) > 7 else {}

            raw_gt = None
            if len(cols) > 9:
                fmt = cols[8].split(":")
                if "GT" in fmt:
                    values = cols[9].split(":")
                    idx = fmt.index("GT")
                    if idx < len(values):
                        raw_gt = values[idx]

            for i, alt in enumerate(alts, start=1):
                if alt in ("*", ".") or alt.startswith("<"):
                    continue  # spanning deletions / symbolic alleles are out of scope
                key = (chrom, pos, ref, alt)
                if wanted is not None and key not in wanted:
                    continue
                allele_info = dict(info)
                if len(alts) > 1:
                    for name in per_allele:
                        parts = info.get(name, "").split(",")
                        if len(parts) == len(alts):
                            allele_info[name] = parts[i - 1]
                gt = _remap_genotype(raw_gt, i) if raw_gt else None
                yield Variant(chrom, pos, ref, alt, vid, qual, flt, allele_info, gt)
