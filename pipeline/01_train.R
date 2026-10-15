# Step 1. Train the local clock and the nine-marker comparator in NHANES III,
# evaluate both in NHANES IV.
suppressPackageStartupMessages({ library(BioAge); library(dplyr) })
source(file.path(PROJECT_ROOT, "R", "config.R"))
source(file.path(PROJECT_ROOT, "R", "train.R"))

local     <- train_clock(PANEL_LOCAL)
reference <- train_clock(PANEL_REFERENCE)

perf <- evaluate_nhanes(local, reference)
print(perf, row.names = FALSE)
write.csv(perf, file.path(PROJECT_ROOT, "final_output", "nhanes_performance.csv"), row.names = FALSE)
saveRDS(list(local = local$fit, reference = reference$fit),
        file.path(PROJECT_ROOT, "data_derived", "models.rds"))
