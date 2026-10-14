# Data governance

**No patient-level data enters this repository. Ever.** Git remembers everything that was ever committed, so there is no "temporary" commit.

## What counts as patient-level data

Any file where each row is one person, including:

- identifiers of any kind (patient id, accession or protocol number),
- individual biomarker values, ages, sexes or dates,
- any individual score (biological age, Z score, percentile).

An **aggregate** (mean, SD, n, RMSE, effect size by group, model coefficients, percentile cut points) is not patient-level data. Aggregates and model parameters are what `registry/` accepts. This repository does not ship any study results: outputs are generated locally when you run the pipeline.

## Where things live

| Folder | Content | Published |
|---|---|---|
| `data_raw/` | Raw laboratory exports | No, git-ignored |
| `data_derived/` | Trained models, scored cohort, anything with a patient id, per-patient reference distributions | No, git-ignored (whole folder) |
| `final_output/` | Outputs generated locally by the pipeline | No, git-ignored |
| `registry/` | Model parameters and aggregate performance submitted by laboratories | Yes, validated automatically |

## Rules for contributions

1. Run `python scripts/build_registry.py --check` before opening a pull request. It rejects keys that look like identifiers and long lists of raw numbers.
2. Report group sizes of at least 100. For smaller groups, report nothing.
3. Do not publish an empirical distribution of individual scores. Publish percentile cut points instead.
4. If you believe patient data was committed by mistake, do not just delete the file in a new commit: contact the maintainers so the history can be rewritten.
