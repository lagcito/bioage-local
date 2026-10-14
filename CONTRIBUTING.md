# Contributing

Two kinds of contribution are welcome.

**1. Add your laboratory's clock to the registry.** Copy `registry/template.json` to `registry/clocks/<your-id>.json`, fill it in, run `python scripts/build_registry.py --check`, and open a pull request. Read `DATA_GOVERNANCE.md` first: only model parameters and aggregates are accepted. New entries appear as *self-reported*; only maintainers set *verified*.

**2. Improve the pipeline or the site.** The R pipeline lives in `pipeline/` and `R/`; the static site is in `docs/` (plain HTML, CSS and JavaScript, no build step). `docs/clock.js` holds the scoring maths and mirrors `clinical/score_patient.R`; if you change one, change the other and compare results.

Scope reminder: a local score is a position among local peers. Contributions must not describe it as a mortality prediction or a measure of the pace of aging unless outcome data support that.
