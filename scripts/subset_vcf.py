"""Cut a VCF down to a manageable subset: one chromosome and/or a variant cap.

Standalone, dependency-free (stdlib only) -- reads .vcf or .vcf.gz, writes the
same. Useful to shrink a whole-genome VCF (millions of variants) down to a
size that ``varlens --freq-api`` can query without hitting --max-api-variants.

Usage:
    python scripts/subset_vcf.py INPUT OUTPUT [--chrom 21] [--max-variants 1500]

Examples:
    # First 1500 variants of the whole file (fast, arbitrary region: chr1 start)
    python scripts/subset_vcf.py NA12878.vcf.gz small.vcf.gz --max-variants 1500

    # All of chromosome 21 (usually tens of thousands of variants -- combine
    # with --max-variants if you need it under 2000 for --freq-api)
    python scripts/subset_vcf.py NA12878.vcf.gz chr21.vcf.gz --chrom 21 --max-variants 1500
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def _open_read(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else open(path, "rt", encoding="utf-8")


def _open_write(path: Path):
    return gzip.open(path, "wt", encoding="utf-8") if path.suffix == ".gz" else open(path, "wt", encoding="utf-8")


def _matches_chrom(line_chrom: str, wanted_chrom: str) -> bool:
    a = line_chrom[3:] if line_chrom.lower().startswith("chr") else line_chrom
    b = wanted_chrom[3:] if wanted_chrom.lower().startswith("chr") else wanted_chrom
    return a == b


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--chrom", help="keep only this chromosome, e.g. 21 or chr21")
    parser.add_argument("--max-variants", type=int, help="stop after this many data lines")
    args = parser.parse_args(argv)

    if not args.input.exists():
        print(f"error: input not found: {args.input}", file=sys.stderr)
        return 2

    kept = 0
    with _open_read(args.input) as src, _open_write(args.output) as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            if args.max_variants is not None and kept >= args.max_variants:
                break
            chrom = line.split("\t", 1)[0]
            if args.chrom and not _matches_chrom(chrom, args.chrom):
                continue
            dst.write(line)
            kept += 1

    print(f"wrote {kept} variants -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
