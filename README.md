# varlens

[![tests](https://github.com/Print-Samuel/varlens/actions/workflows/ci.yml/badge.svg)](https://github.com/Print-Samuel/varlens/actions)

A small, dependency-free command-line tool that takes a **VCF**, annotates each variant with
**ClinVar** assertions and **population allele frequency**, applies a transparent set of
ACMG/AMP-inspired rules, and produces a sortable **HTML report** (plus optional TSV).

> **Educational use only.** varlens is a rule-based prioritisation aid built for learning and
> portfolio purposes. It is **not** a clinical classifier and must not be used for diagnosis or
> medical decisions.

[Example report](docs/example_report.html) (generated from synthetic data)

## Why this exists

Clinical genome/exome interpretation pipelines are usually black boxes. varlens keeps every step
readable: each variant's tier comes with the exact evidence lines that produced it, so you can see
*why* a variant was ranked where it was.

## Install and run the demo

```bash
git clone https://github.com/Print-Samuel/varlens.git
cd varlens
pip install -e ".[dev]"

varlens data/example_patient.vcf \
  --clinvar data/example_clinvar.vcf \
  --freq data/example_freq.vcf \
  --out report.html --tsv report.tsv
```

The files in `data/` are **synthetic** (fictional genes `DEMO1`-`DEMO5`, fictional coordinates).
They exist only to exercise the pipeline.

## Scanning part of a genome instead of all of it

A whole genome is millions of variants; `--freq-api` is capped at 4000 by default precisely
because looking each one up over the network, one at a time, doesn't scale to that. Use
`--region` to restrict the *input* VCF itself to one chromosome or a coordinate window --
applied while streaming, so it works directly on a whole-genome file, no pre-processing needed:

```bash
varlens HG001.vcf.gz --region 21 --clinvar clinvar.vcf.gz --freq-api --out report.html
varlens HG001.vcf.gz --region 21:5000000-10000000 --clinvar clinvar.vcf.gz --out report.html
```

`scripts/subset_vcf.py` still exists for writing a smaller VCF *file* to disk (useful if you
want to keep a fixed subset around); `--region` is the quicker option when you just want to
point one run of `varlens` at part of a genome.

## Not re-paying for the same gnomAD lookups every run

`--freq-api` writes every variant it resolves into a local cache file
(`.varlens_gnomad_cache.json` by default, `--freq-api-cache` to change the path, `--no-cache` to
turn it off). Re-running on the same variants -- for example after tweaking `--ba1`/`--bs1`/`--pm2`
and wanting to see how the tiers change -- reuses the cache instead of hitting the network again.
Before starting, `varlens` prints how many variants still need a fresh lookup and a rough time
estimate, so a big run doesn't surprise you.

## Using real data

1. Get ClinVar's VCF for the **same genome build** as your sample (GRCh37 or GRCh38) from the
   [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz) (~80-120 MB,
   check the current path/filename before relying on it).
2. For population frequency, you have two options:
   - **`--freq-api`** (no download): query gnomAD's free public GraphQL API one variant at a
     time. No API key needed. Meant for a single patient's VCF (or one `--region` of it), not a
     whole genome at once -- capped at 4000 variants by default (`--max-api-variants`), paced
     conservatively and slowed down further if gnomAD signals its rate limit (see Limitations).
     Requires network access from wherever you run `varlens`.
     ```bash
     varlens sample.vcf --clinvar clinvar.vcf.gz --freq-api --out report.html
     ```
   - **`--freq path/to/file.vcf.gz`**: use a downloaded gnomAD sites VCF or subset of one
     ([gnomAD downloads](https://gnomad.broadinstitute.org/downloads)). Full gnomAD files are
     tens of gigabytes; restrict to one chromosome (chr21 or chr22 are the smallest) if you
     want to keep the download manageable. Use this for many variants or when you don't have
     network access while running the tool (only the download needs the network, not the run).
3. Normalise your VCF so indels are left-aligned and multi-allelics are split consistently:
   `bcftools norm -m -any -f reference.fa sample.vcf.gz`. Variants are matched by
   `(chrom, pos, ref, alt)`, so representation must agree -- this matters whichever frequency
   source you use.
4. Optional: annotate consequences first (e.g. SnpEff). varlens reads the standard `ANN` field if
   present; otherwise it falls back to ClinVar's `MC` field for known variants.

A real, freely downloadable human genome to try this on: NIST's
[Genome in a Bottle NA12878/HG001](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/),
a well-characterised reference sample (millions of real variants, GRCh38, no login required).
It's a *healthy* reference sample, so it won't contain a real pathogenic variant to check the
"High priority" tier against -- for that, `scripts/make_positive_control.py` pulls a real,
review-starred pathogenic variant straight out of your own downloaded ClinVar file and injects
it into a copy of your patient VCF as a positive control, with the exact coordinates coming from
ClinVar itself rather than anyone's memory:

```bash
python scripts/make_positive_control.py --clinvar clinvar.vcf.gz --chrom 21 --list
python scripts/make_positive_control.py --clinvar clinvar.vcf.gz --chrom 21 \
    --patient HG001.vcf.gz --out HG001_positive_control.vcf.gz
```

Never commit real patient data. `.gitignore` excludes `*.vcf.gz` by default.

## How prioritisation works

| Evidence | Rule (defaults, configurable) | Effect |
|---|---|---|
| BA1 | allele frequency >= 0.05 | Low priority, overrides everything |
| BS1 | allele frequency >= 0.01 | Low priority (or flags a conflict with a ClinVar pathogenic call) |
| PM2 (supporting) | absent or AF <= 0.0001 | Supporting evidence line |
| LoF (PVS1-like) | stop-gain, frameshift, canonical splice, start-lost | Rare LoF with no ClinVar entry -> manual review |
| ClinVar | P/LP with >= 2 review stars | High priority |
| ClinVar | P/LP with < 2 stars, or conflicting | Manual review |
| ClinVar | Benign/LB with >= 1 star | Low priority |

ClinVar is reported as its own evidence line rather than as ACMG PP5/BP6, since ClinGen's
Sequence Variant Interpretation group recommended retiring those criteria.

Tiers: **High priority**, **Needs manual review**, **Uncertain / unannotated**, **Low priority**.

## Limitations (on purpose, and on the roadmap)

- Matching is exact on `(chrom, pos, ref, alt)`; no automatic normalisation.
- Only the first sample in the VCF is used; no trio/segregation logic.
- No gene-disease mechanism data, so LoF is flagged, not asserted.
- No in-silico predictors, conservation, splicing or functional evidence.
- `--freq-api` looks up one variant per HTTP request; fine for a single patient's VCF or one
  `--region`, too slow for a whole genome at once (that's what `--max-api-variants` guards
  against). gnomAD signals its rate limit as an HTTP 400 with a text message rather than the
  standard 429, which this handles, but the exact throughput you'll get on a given day is
  outside this project's control. Streams ClinVar/gnomAD VCF files sequentially otherwise; large
  files are slow (tabix-indexed lookup is planned).
- `--region` only accepts a chromosome or a numeric coordinate window, not a gene symbol yet --
  you need the gene's coordinates from somewhere else (e.g. Ensembl) until gene-symbol lookup
  lands (see Roadmap).
- Not validated against a truth set. A benchmark against ClinVar-labelled variants is planned.

## Development

```bash
pip install -e ".[dev]"
pytest -q
```

CI runs the test suite on Python 3.9, 3.11 and 3.12.

## Roadmap

- [ ] Tabix-indexed ClinVar/gnomAD lookup
- [ ] `bcftools norm` integration / built-in left-alignment
- [ ] Benchmark of the rules on ClinVar-labelled variants (precision/recall per tier)
- [ ] Gene-level LoF-intolerance data (e.g. pLI/LOEUF) to refine the LoF flag
- [ ] Trio mode (de novo / compound heterozygous)
- [ ] `--gene SYMBOL` as a shortcut for `--region`, resolved via a bundled or fetched coordinate table
- [ ] Publish to PyPI (`pip install varlens`) instead of git-clone-and-editable-install
- [ ] Test coverage for `report.py` (HTML/TSV output) -- currently exercised only indirectly via
      the end-to-end pipeline test

## License

MIT
