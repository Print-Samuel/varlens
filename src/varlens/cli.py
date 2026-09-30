"""Command-line interface: ``varlens sample.vcf --clinvar clinvar.vcf.gz --out report.html``."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__
from .annotate import annotate, load_clinvar, load_frequencies
from .gnomad_api import GnomadAPIError, query_frequencies, variant_id
from .report import write_html, write_tsv
from .rules import TIER_LABEL, TIER_ORDER, Thresholds, prioritise
from .vcf import parse_region, read_vcf


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="varlens",
        description="Annotate a VCF with ClinVar and population frequency and prioritise variants "
        "(educational tool, not for clinical use).",
    )
    p.add_argument("vcf", help="input VCF (.vcf or .vcf.gz); first sample is used")
    p.add_argument(
        "--region",
        help="restrict to one chromosome or a coordinate window, e.g. '21' or "
        "'21:1000000-2000000'. Applied while streaming, so it works directly on a "
        "whole-genome VCF -- no need to pre-cut it with scripts/subset_vcf.py first.",
    )
    p.add_argument("--clinvar", help="ClinVar VCF (same genome build as the input)")
    p.add_argument("--freq", help="population-frequency VCF, e.g. a gnomAD subset")
    p.add_argument("--af-field", default="AF", help="INFO field with allele frequency (default: AF)")
    p.add_argument(
        "--freq-api",
        action="store_true",
        help="look up population frequency from gnomAD's public API instead of a local "
        "file -- no download needed, but slower and requires network access. Not for "
        "large batches (see --max-api-variants).",
    )
    p.add_argument(
        "--dataset",
        default="gnomad_r4",
        help="gnomAD dataset queried by --freq-api (default: gnomad_r4)",
    )
    p.add_argument(
        "--max-api-variants",
        type=int,
        default=4000,
        help="safety cap on the number of variants looked up via --freq-api (default 4000). "
        "Raising this raises how long a run can take, roughly linearly -- see the printed "
        "time estimate before committing to a big one.",
    )
    p.add_argument(
        "--freq-api-delay",
        type=float,
        default=2.0,
        help="starting seconds between --freq-api requests (default 2.0). Automatically "
        "increases (up to --freq-api-max-delay) if gnomAD reports its query rate limit.",
    )
    p.add_argument(
        "--freq-api-max-delay",
        type=float,
        default=30.0,
        help="ceiling in seconds for the adaptive pace above (default 30.0)",
    )
    p.add_argument(
        "--freq-api-cache",
        default=".varlens_gnomad_cache.json",
        help="JSON file caching gnomAD lookups across runs, so re-running (e.g. after "
        "tweaking --ba1/--bs1/--pm2) doesn't repeat old queries (default: "
        ".varlens_gnomad_cache.json in the current folder)",
    )
    p.add_argument("--no-cache", action="store_true", help="disable the gnomAD lookup cache entirely")
    p.add_argument("--out", default="varlens_report.html", help="HTML report path")
    p.add_argument("--tsv", help="also write a TSV table")
    p.add_argument("--ba1", type=float, default=0.05, help="stand-alone benign AF (default 0.05)")
    p.add_argument("--bs1", type=float, default=0.01, help="too-common-for-disease AF (default 0.01)")
    p.add_argument("--pm2", type=float, default=0.0001, help="'very rare' AF (default 0.0001)")
    p.add_argument("--version", action="version", version=f"varlens {__version__}")
    return p


def _load_cache(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"warning: could not read cache {path} ({exc}); starting fresh", file=sys.stderr)
        return {}


def _save_cache(path: Path, cache: dict) -> None:
    try:
        path.write_text(json.dumps(cache), encoding="utf-8")
    except OSError as exc:
        print(f"warning: could not write cache {path}: {exc}", file=sys.stderr)


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.freq and args.freq_api:
        print("error: use either --freq or --freq-api, not both", file=sys.stderr)
        return 2

    for label, path in (("VCF", args.vcf), ("ClinVar", args.clinvar), ("frequency", args.freq)):
        if path and not Path(path).exists():
            print(f"error: {label} file not found: {path}", file=sys.stderr)
            return 2

    region = None
    if args.region:
        try:
            region = parse_region(args.region)
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2

    variants = list(read_vcf(args.vcf, region=region))
    if not variants:
        hint = f" in region {args.region!r}" if args.region else ""
        print(f"error: no usable variants found in the input VCF{hint}", file=sys.stderr)
        return 2
    wanted = {v.key for v in variants}

    clinvar = load_clinvar(args.clinvar, wanted) if args.clinvar else {}

    if args.freq_api:
        cache_path = None if args.no_cache else Path(args.freq_api_cache)
        cache = _load_cache(cache_path) if cache_path else None

        to_fetch = wanted if cache is None else {k for k in wanted if variant_id(k) not in cache}
        if to_fetch:
            est_min = len(to_fetch) * args.freq_api_delay / 60
            cached_note = f" ({len(wanted) - len(to_fetch)} already cached)" if cache else ""
            print(
                f"looking up {len(to_fetch)} variant(s) via gnomAD's API{cached_note} -- "
                f"at least ~{est_min:.1f} min at the current pace, likely more if rate-limited.",
                file=sys.stderr,
            )

        def _progress(done: int, total: int) -> None:
            print(f"\rquerying gnomAD: {done}/{total}", end="", file=sys.stderr, flush=True)

        skipped: list = []

        def _skip(key, reason: str) -> None:
            skipped.append((key, reason))

        try:
            freqs = query_frequencies(
                wanted,
                dataset=args.dataset,
                max_variants=args.max_api_variants,
                delay=args.freq_api_delay,
                max_delay=args.freq_api_max_delay,
                on_progress=_progress,
                on_skip=_skip,
                cache=cache,
            )
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        except GnomadAPIError as exc:
            print(file=sys.stderr)
            print(f"error: {exc}", file=sys.stderr)
            return 2
        finally:
            if cache_path is not None and cache is not None:
                _save_cache(cache_path, cache)
        print(file=sys.stderr)  # end the progress line
        if skipped:
            print(f"warning: {len(skipped)} variant(s) skipped by gnomAD lookup:", file=sys.stderr)
            for key, reason in skipped[:5]:
                chrom, pos, ref, alt = key
                print(f"  {chrom}:{pos} {ref}>{alt} -- {reason}", file=sys.stderr)
            if len(skipped) > 5:
                print(f"  ... and {len(skipped) - 5} more", file=sys.stderr)
    elif args.freq:
        freqs = load_frequencies(args.freq, wanted, args.af_field)
    else:
        freqs = {}

    annotated = annotate(variants, clinvar, freqs)
    th = Thresholds(ba1=args.ba1, bs1=args.bs1, pm2=args.pm2)
    rows = [(a, prioritise(a, th)) for a in annotated]

    meta = {"Input": Path(args.vcf).name}
    if args.region:
        meta["Region"] = args.region
    if args.clinvar:
        meta["ClinVar"] = Path(args.clinvar).name
    if args.freq_api:
        meta["Frequencies"] = f"gnomAD API ({args.dataset})"
    elif args.freq:
        meta["Frequencies"] = Path(args.freq).name
    write_html(rows, args.out, meta)
    if args.tsv:
        write_tsv(rows, args.tsv)

    counts = {t: sum(1 for _, r in rows if r.tier == t) for t in TIER_ORDER}
    print(f"{len(rows)} variants annotated -> {args.out}")
    for tier in TIER_ORDER:
        print(f"  {TIER_LABEL[tier]:<26} {counts[tier]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
