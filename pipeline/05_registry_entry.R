# Step 5. Draft registry entry (parameters and aggregates only) for this laboratory.
suppressPackageStartupMessages(library(BioAge))
for (f in c("config", "registry_entry")) source(file.path(PROJECT_ROOT, "R", paste0(f, ".R")))

models <- readRDS(file.path(PROJECT_ROOT, "data_derived", "models.rds"))
rc     <- readRDS(file.path(PROJECT_ROOT, "data_derived", "local_recentered.rds"))
perf   <- read.csv(file.path(PROJECT_ROOT, "final_output", "local_performance.csv"))

entry <- build_registry_entry(models$local, rc$rec, perf, rc$df)
path  <- file.path(PROJECT_ROOT, "final_output", "registry_entry_draft.json")
jsonlite::write_json(entry, path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(">> Wrote ", path, "\n   Review it, fill in the laboratory fields, then validate it with:\n",
        "   python scripts/build_registry.py --check")
