# BioAge Local: pipeline and registry of locally referenced biological age clocks

An open pipeline to implement a PhenoAge-type biological age clock **without local survival data**, plus a public registry where laboratories can publish their clock and the local reference they built for it.

- **Pipeline** (`pipeline/`, `R/`, `clinical/`): trains a six-marker PhenoAge-derived model (`pheno_max_v2`) with the [BioAge](https://github.com/dayoonkwon/BioAge) toolkit on NHANES III, evaluates it on NHANES IV, applies it to a local cohort, recenters it for age and sex (Method B), and validates the local score against biomarkers the clock never used.
- **Registry and site** (`registry/`, `docs/`, `scripts/`): a static web page, served with GitHub Pages, that lists clocks submitted by laboratories and includes an in-browser calculator. No server, no database, no accounts.
- **No data of any kind is shipped.** Outputs of the pipeline are written to `final_output/` when you run it and are git-ignored.

> **Research use only.** Nothing here is a medical device or clinical advice. A local score describes a position among local peers. It does not predict mortality and it is not a measure of the pace of aging. Registry entries are self-reported unless marked verified.

## Patient data is not in this repository

Raw laboratory data and every per-patient artifact stay outside git. See [`DATA_GOVERNANCE.md`](DATA_GOVERNANCE.md). To run the pipeline you need your own laboratory export in `data_raw/` (a template for the clinical calculator input is in `clinical/example_patients_template.csv`).

## The registry site

```bash
python scripts/build_registry.py     # validate entries and build docs/registry.json
cd docs && python -m http.server     # open http://localhost:8000
```

Once published it lives at `https://lagcito.github.io/bioage-local/` (repository: `https://github.com/lagcito/bioage-local`). To publish it: in the repository settings, set **Pages → Source: GitHub Actions**. The workflow in `.github/workflows/pages.yml` rebuilds the registry and deploys `docs/` on every push to `main`.

### Add your laboratory

1. Copy `registry/template.json` to `registry/clocks/<your-id>.json` and fill it in.
2. Run `python scripts/build_registry.py --check`. It rejects anything that looks like patient-level data (identifier-like keys, long lists of raw numbers, groups smaller than 100).
3. Open a pull request. A check runs automatically and maintainers review it.

Each entry holds: the biomarkers and units, the Gompertz model parameters, the local recentering coefficients (`BA ~ age + sex`, residual SD), empirical Z cut points, and aggregate performance. The entry named `example-demo-lab` is synthetic and exists only to demonstrate the interface; copy `registry/template.json` for your own.

## The calculator

`docs/clock.js` computes biological age, the local Z score, an approximate percentile and the P1 to P5 group in the browser. It mirrors `clinical/score_patient.R`; on four test patients the two agree within the R script's two-decimal rounding.

## Running the pipeline

```r
install.packages(c("dplyr", "purrr", "survival", "ggplot2", "tidyr",
                   "tibble", "gridExtra", "corrplot", "flexsurv"))
devtools::install_github("dayoonkwon/BioAge")
# put your laboratory CSV exports in data_raw/  (not versioned)
source("pipeline/run_all.R")
```

| Step | Script | What it does |
|---|---|---|
| 1 | `pipeline/01_train.R` | Trains PhenoAge, KDM, HD and sub-clocks on NHANES III |
| 2 | `pipeline/02_score_argentina.R` | Scores the local cohort, RMSE by age and sex, Bland-Altman |
| 3 | `pipeline/03_recenter.R` | Method B recentering, |Z| > 5 censoring, empirical cut points |
| 4 | `pipeline/04_validate.R` | External construct validity, P1 to P5 profiles, sex replication |
| 5 | `pipeline/05_extract_manuscript.R` | Writes the numeric digest to `final_output/` |

Steps share the R session state, so run them in order (`run_all.R` does).

## What the first entry is, and is not

`pheno_max_v2` reuses the training method of Levine et al. (2018) but keeps six of the nine markers (glucose, creatinine, lymphocyte %, MCV, RDW, WBC); albumin, CRP and alkaline phosphatase were dropped for low coverage in the source laboratory. It is a derived instrument, not the original PhenoAge. The survival parameters come from NHANES; the local recentering defines a reference and does not recalibrate mortality risk. The source cohort is a single-laboratory referral sample and national representativeness was not assessed.

## References

- Levine ME et al. An epigenetic biomarker of aging for lifespan and healthspan. *Aging* 2018;10:573-591.
- Liu Z et al. A new aging measure captures morbidity and mortality risk across diverse subpopulations from NHANES IV. *PLoS Med* 2018;15:e1002718.
- Kwon D, Belsky DW. A toolkit for quantification of biological age from blood chemistry and organ function test data: BioAge. *GeroScience* 2021;43:2795-2808.

## License

Code: MIT (`LICENSE`). Registry entries carry the license declared in each file.
