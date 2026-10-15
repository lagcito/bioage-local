###############################################################################
# run_all.R - run the whole pipeline, from the repository root:
#   Rscript pipeline/run_all.R
#
# Needs a local laboratory export in data_raw/ (never committed; see README).
# Writes aggregate tables to final_output/ and patient-level intermediates to
# data_derived/. Both folders are git-ignored.
###############################################################################

PROJECT_ROOT <- normalizePath(".")
if (!file.exists(file.path(PROJECT_ROOT, "R", "config.R"))) {
  stop("Run this script from the repository root.")
}
dir.create(file.path(PROJECT_ROOT, "final_output"), showWarnings = FALSE)
dir.create(file.path(PROJECT_ROOT, "data_derived"), showWarnings = FALSE)

for (step in c("01_train.R", "02_score_local.R", "03_recenter.R",
               "04_validate.R", "05_registry_entry.R")) {
  message("\n== ", step, " ==")
  source(file.path(PROJECT_ROOT, "pipeline", step), local = new.env())
}
message("\nDone. Tables are in final_output/; the draft registry entry is final_output/registry_entry_draft.json")
