# varlens blind benchmark

Run on `clinvar.vcf.gz` with varlens 0.1.0 on 2026-10-14.
Labelled variants with at least 1 review star(s): 335306 pathogenic, 1382522 benign.

| ClinVar label | high | review | uncertain | low |
|---|---|---|---|---|
| pathogenic | 0 | 262834 | 72452 | 20 |
| benign | 0 | 2827 | 1251963 | 127732 |

| Question | Result |
|---|---|
| Pathogenic variants flagged for review | 78.4% (262834/335306, 95% CI 78.2% to 78.5%) |
| Benign variants flagged for review | 0.2% (2827/1382522, 95% CI 0.2% to 0.2%) |
| Pathogenic variants dismissed as low priority | 0.006% (20/335306, 95% CI 0.004% to 0.009%) |
| Benign variants recognised as low priority | 9.2% (127732/1382522, 95% CI 9.2% to 9.3%) |
| Positive predictive value of the review flag | 98.9% (262834/265661, 95% CI 98.9% to 99.0%) |

## How to read these numbers

- The review flag mostly tracks loss of function: 78.9% of pathogenic variants have a loss-of-function consequence, against 0.3% of benign ones. ClinVar curators also use consequence as evidence, so part of this agreement is by construction.
- 72452 pathogenic variants (21.6%) stay in the uncertain tier, and 70820 pathogenic variants have no loss-of-function consequence. Frequency and consequence alone cannot separate these from benign variants; this is where missense-aware evidence (predictors, a trained model) would have to add value.
- The positive predictive value depends on prevalence. ClinVar's answer key is 19.5% pathogenic; among the rare variants in a real sample, pathogenic ones are far scarcer. With the same sensitivity and false-alarm rate, a prevalence of 0.1% would give a PPV of about 27.7%.
- Frequency data exists for 5.7% of pathogenic and 38.6% of benign variants, so the frequency rules rarely get to act on pathogenic variants, and a missing frequency is treated as rare.

## Method and caveats

- Blind mode: the ClinVar label is never given to the rules; it is only the answer key.
- Frequencies come from ClinVar's AF_EXAC / AF_TGP / AF_ESP fields (older, smaller cohorts than gnomAD). A variant missing from them is treated as not observed (rare).
- The 'high' tier needs a ClinVar assertion, so it cannot occur in blind mode.
- Rates are conditional on the true class, with 95% Wilson score intervals; positive predictive value depends on the pathogenic/benign mix (prevalence is reported).
- Labels are ClinVar's own and are not ground truth: they can be revised, and reflect which genes and variants get submitted.
