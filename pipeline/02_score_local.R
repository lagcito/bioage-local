# Step 2. Score the local cohort and describe how the clock behaves in it.
suppressPackageStartupMessages({ library(BioAge); library(dplyr) })
for (f in c("config", "train", "local_data", "performance")) source(file.path(PROJECT_ROOT, "R", paste0(f, ".R")))

models <- readRDS(file.path(PROJECT_ROOT, "data_derived", "models.rds"))
prep <- prepare_local_data(read_local_data())
df <- prep$df
message(">> ", nrow(df), " eligible records")

df$ba  <- score_clock(df, PANEL_LOCAL, models$local)
scored <- is.finite(df$ba)
message(">> ", sum(scored), " scored (", round(100 * mean(scored), 1), "%)")

out <- function(x, f) write.csv(x, file.path(PROJECT_ROOT, "final_output", f), row.names = FALSE)
perf <- cbind(clock = CLOCK_NAME, age_metrics(df$ba, df$age))
print(perf, row.names = FALSE)
out(perf,                                  "local_performance.csv")
out(rmse_by_age_group(df$ba, df$age),      "local_rmse_by_age.csv")
out(bland_altman(df$ba, df$age),           "local_bland_altman.csv")
out(scored_vs_unscored(df, scored),        "local_scored_vs_unscored.csv")

saveRDS(list(df = df, raw = prep$raw), file.path(PROJECT_ROOT, "data_derived", "local_scored.rds"))
