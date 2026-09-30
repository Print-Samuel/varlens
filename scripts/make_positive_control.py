"""Build a 'positive control' VCF by injecting one real, known-pathogenic ClinVar
variant (from chromosome 20 or 21, by default) into a copy of a patient VCF.

Why: public reference genomes (like GIAB's HG001/NA12878) are healthy samples
by design, so they won't contain a real pathogenic variant to test the "High
priority" tier against. This script pulls a REAL variant's exact
chrom/pos/ref/alt straight out of YOUR downloaded clinvar.vcf.gz -- no
coordinate is guessed or typed from memory -- and adds it as a heterozygous
call, so you can confirm the pipeline actually flags it correctly.

Usage:
    # 1) See a few candidates first (doesn't write anything):
    python scripts/make_positive_control.py --clinvar clinvar.vcf.gz --chrom 21 --list

    # 2) Inject the Nth candidate (default: the first/best one) into your patient VCF:
    python scripts/make_positive_control.py --clinvar clinvar.vcf.gz --chrom 21 \
        --patient data_real/HG001_GRCh38_chr21.vcf.gz \
        --out data_real/HG001_positive_control.vcf.gz
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import List, NamedTuple, Union

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from varlens.annotate import clinvar_category, review_stars  # noqa: E402
from varlens.vcf import Variant, normalize_chrom, read_vcf  # noqa: E402
from urllib.parse import unquote  # noqa: E402


class Candidate(NamedTuple):
    variant: Variant
    gene: str
    disease: str
    stars: int
    significance: str


def _gene_from_geneinfo(geneinfo: str) -> str:
    return geneinfo.split("|")[0].split(":")[0] if geneinfo else "?"


def find_candidates(clinvar_path: Union[str, Path], chrom: str, min_stars: int) -> List[Candidate]:
    target = normalize_chrom(chrom)
    out: List[Candidate] = []
    for v in read_vcf(clinvar_path):
        if v.chrom != target:
            continue
        sig = unquote(v.info.get("CLNSIG", ""))
        if clinvar_category(sig) != "pathogenic":
            continue
        stars = review_stars(unquote(v.info.get("CLNREVSTAT", "")))
        if stars < min_stars:
            continue
        out.append(Candidate(
            v,
            _gene_from_geneinfo(v.info.get("GENEINFO", "")),
            unquote(v.info.get("CLNDN", "")).replace("_", " "),
            stars,
            sig,
        ))
    out.sort(key=lambda c: -c.stars)
    return out


def _open_read(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path, "rt", encoding="utf-8")


def _open_write(path: Path):
    return gzip.open(path, "wt", encoding="utf-8") if path.suffix == ".gz" else open(path, "wt", encoding="utf-8")


def inject(patient_path: Path, out_path: Path, candidate: Candidate, genotype: str) -> None:
    v = candidate.variant
    new_line = "\t".join([
        v.chrom, str(v.pos), f"pc_{v.id}", v.ref, v.alt, "99", "PASS",
        f"NOTE=synthetic_positive_control_from_clinvar;GENE={candidate.gene}",
        "GT", genotype,
    ]) + "\n"

    saw_header_end = False
    with _open_read(patient_path) as src, _open_write(out_path) as dst:
        for line in src:
            dst.write(line)
            if line.startswith("#CHROM"):
                saw_header_end = True
        if not saw_header_end:
            raise ValueError("patient VCF has no #CHROM header line -- is it a valid VCF?")
        dst.write(new_line)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clinvar", required=True, type=Path)
    p.add_argument("--chrom", default="21")
    p.add_argument("--min-stars", type=int, default=2)
    p.add_argument("--list", action="store_true", help="list candidates and exit, without writing anything")
    p.add_argument("--index", type=int, default=0, help="which candidate to inject (0 = first/best)")
    p.add_argument("--patient", type=Path, help="patient VCF to copy and inject into")
    p.add_argument("--out", type=Path, help="output path for the new positive-control VCF")
    p.add_argument("--genotype", default="0/1", help="GT to assign the injected variant (default 0/1)")
    args = p.parse_args(argv)

    if not args.clinvar.exists():
        print(f"error: ClinVar file not found: {args.clinvar}", file=sys.stderr)
        return 2

    print(f"scanning {args.clinvar} for chromosome {args.chrom} pathogenic variants "
          f"(>= {args.min_stars} stars)...", file=sys.stderr)
    candidates = find_candidates(args.clinvar, args.chrom, args.min_stars)
    if not candidates:
        print(f"error: no pathogenic ClinVar variants found on chromosome {args.chrom} "
              f"with >= {args.min_stars} review stars. Try --chrom 20, or lower --min-stars.",
              file=sys.stderr)
        return 2

    for i, c in enumerate(candidates[:15]):
        marker = "->" if i == args.index else "  "
        print(f"{marker} [{i}] {c.variant.locus}  gene={c.gene}  stars={c.stars}  "
              f"{c.significance}  ({c.disease or 'condition not specified'})")
    if len(candidates) > 15:
        print(f"   ... and {len(candidates) - 15} more")

    if args.list:
        return 0

    if not args.patient or not args.out:
        print("error: --patient and --out are required unless you pass --list", file=sys.stderr)
        return 2
    if not args.patient.exists():
        print(f"error: patient VCF not found: {args.patient}", file=sys.stderr)
        return 2
    if args.index >= len(candidates):
        print(f"error: --index {args.index} out of range (only {len(candidates)} candidates)", file=sys.stderr)
        return 2

    chosen = candidates[args.index]
    inject(args.patient, args.out, chosen, args.genotype)
    print(f"\ninjected {chosen.variant.locus} (gene {chosen.gene}, {chosen.significance}) "
          f"as {args.genotype} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
