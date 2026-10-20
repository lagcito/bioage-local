# BioAge Local: pipeline and registry of locally referenced biological age clocks

An open pipeline to implement a PhenoAge-type biological age clock **without local survival data**, and a public registry where laboratories can publish the clock they use and the local reference they built for it.

- **Pipeline** (`pipeline/`, `R/`): trains a six-marker PhenoAge-derived model with the [BioAge](https://github.com/dayoonkwon/BioAge) toolkit in NHANES III, evaluates it in NHANES IV, applies it to a local cohort, builds a local age- and sex-specific reference (Method B), and checks the local score against biomarkers the clock never used.
- **Registry and site** (`registry/`, `docs/`, `scripts/`): a static website, served with GitHub Pages, listing the clocks submitted by laboratories, with an in-browser calculator. No server, no database, no accounts.
- **No data of any kind is shipped.** Pipeline outputs are written to `final_output/` and `data_derived/` when you run it; both are git-ignored. See [`DATA_GOVERNANCE.md`](DATA_GOVERNANCE.md).

> **Research use only.** Nothing here is a medical device or clinical advice. A local score is a position among local peers for a given age and sex. It does not predict mortality and it is not a measure of the pace of aging. Registry entries are self-reported unless marked verified.

## Run the pipeline

```r
install.packages(c("dplyr", "jsonlite"))
devtools::install_github("dayoonkwon/BioAge")
# put your laboratory export in data_raw/ (never committed), then, from the repository root:
source("pipeline/run_all.R")        # or: Rscript pipeline/run_all.R
```

It takes about half a minute. Edit `R/config.R` to change the file pattern, the identifier column or the clock markers.

### Input format

One row per laboratory record, as CSV. Required columns and units:

| Column | Unit |
|---|---|
| `Protocolo` (set `LOCAL_ID_COLUMN`) | any unique record identifier |
| `age` | years |
| `sex` | `F` or `M` |
| `glucose` | mg/dL |
| `creatinine` | mg/dL |
| `lymphocyte` | % |
| `mcv` | fL |
| `rdw` | % |
| `wbc` | 10^3/µL (cells/µL is detected and converted) |

Optional columns: `albumin` (g/dL), `alp`, `crp` (mg/L), `neutrophils`, `insulin`, `lymphocyte_abs`, and any other analyte listed in `EXTERNAL_CANDIDATES`, which are used for external validation.

### Steps

| Step | Script | What it does | Writes to `final_output/` |
|---|---|---|---|
| 1 | `01_train.R` | Trains the local clock and a nine-marker comparator in NHANES III; evaluates both in NHANES IV, on complete cases and on the participants scored by both | `nhanes_performance.csv` |
| 2 | `02_score_local.R` | Scores the local cohort; agreement with chronological age, RMSE by age group, Bland-Altman limits, scored vs unscored patients | `local_performance.csv`, `local_rmse_by_age.csv`, `local_bland_altman.csv`, `local_scored_vs_unscored.csv` |
| 3 | `03_recenter.R` | Method B: `lm(BA ~ age + sex)`, standardized residual (Z), censoring of extreme Z, empirical percentile cut points and five bands | `recentering_summary.csv` |
| 4 | `04_validate.R` | External construct validity: Cohen d (P5 vs P1) with CI, Spearman correlation, Kruskal-Wallis across bands, Benjamini-Hochberg adjustment; sensitivity without hematology-flagged records; by sex; selection-bias table | `external_validity*.csv`, `ascertainment.csv` |
| 5 | `05_registry_entry.R` | Draft registry entry with the model parameters and aggregate statistics | `registry_entry_draft.json` |

Patient-level intermediates (including the log of censored cases) stay in `data_derived/`.

## The method in four lines

1. The survival model comes from NHANES (where outcomes exist); locally there is no follow-up, so it is not retrained.
2. For each patient, expected biological age is predicted from age and sex in the local cohort; the standardized difference is the local Z score.
3. Percentiles come from the observed distribution of Z, not from a normal curve. P1 is the lowest 10%, P5 the highest 10%.
4. Without outcomes, validity is assessed against biomarkers the clock did not use. P1 and P5 are extreme groups, so Cohen d is an upper bound; the Spearman correlation does not select extremes. Tests ordered because a clinician is concerned can inflate effects, which is why `ascertainment.csv` compares how many patients were measured in each group.

## The registry site

```bash
python scripts/build_registry.py     # validate entries and build docs/registry.json
cd docs && python -m http.server     # open http://localhost:8000
```

Once published it lives at `https://lagcito.github.io/bioage-local/`; the workflow in `.github/workflows/pages.yml` rebuilds the registry and deploys `docs/` on every push to `main`.

### Add your laboratory

1. Run the pipeline and take `final_output/registry_entry_draft.json`, or copy `registry/template.json`.
2. Save it as `registry/clocks/<your-id>.json` and fill in the laboratory fields.
3. Run `python scripts/build_registry.py --check`. It rejects anything that looks like patient-level data (identifier-like keys, long lists of raw numbers, non-finite numbers, oversized values, groups smaller than 100).
4. Open a pull request. A check runs automatically and maintainers review it.

## Calculators

`docs/clock.js` (website) and `clinical/score_patient.R` (R) compute biological age, the local Z score, an approximate percentile and the P1 to P5 group from a registry entry, with the same maths. They agree within rounding.

```r
source("clinical/score_patient.R")
entry <- load_entry("final_output/registry_entry_draft.json")
score_patient(entry, age = 55, sex = "M", glucose = 95, creatinine = 0.9,
              lymph = 30, mcv = 88, rdw = 13.2, wbc = 6.5)
```

## References

- Levine ME et al. An epigenetic biomarker of aging for lifespan and healthspan. *Aging* 2018;10:573-591.
- Liu Z et al. A new aging measure captures morbidity and mortality risk across diverse subpopulations from NHANES IV. *PLoS Med* 2018;15:e1002718.
- Kwon D, Belsky DW. A toolkit for quantification of biological age from blood chemistry and organ function test data: BioAge. *GeroScience* 2021;43:2795-2808.

## License

Code: MIT (`LICENSE`). Registry entries carry the license declared in each file.
