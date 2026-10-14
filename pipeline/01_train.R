##############################################################################
# run_train — Train in NHANES III pheno_max_v2 (Paper-Ready)
#
# UNIFIED PIPELINE: Train → Score → Recenter → Discovery
#
# This script runs the entire analysis in a single source() call:
#   STEP 1. Train clocks on NHANES III → validate on NHANES IV
#
# FOCUS CLOCKS (by Argentine data coverage):
#   - pheno_pheno_max_v2             (57% coverage — 6/9 Levine)
# USAGE:  source("pipeline/run_all.R")           — pipeline completo
#         source("pipeline/01_train.R")           — este paso solo
###############################################################################

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  BIOLOGICAL AGING — COMPLETE ANALYSIS PIPELINE          ║\n")
cat("║  Train (NHANES)this script → Score (Argentina) second → Recenter → Profile║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# ══════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════

# PROJECT_ROOT puede venir ya seteado por pipeline/run_all.R; si este script
# se corre solo (source("pipeline/01_train.R") desde la raíz del repo), se
# calcula subiendo un nivel desde pipeline/.
if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- tryCatch({
    normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
  }, error = function(e) normalizePath(".."))
}
setwd(PROJECT_ROOT)
message(">> Project root: ", PROJECT_ROOT)

DATA_DERIVED_DIR <- file.path(PROJECT_ROOT, "data_derived")
if (!dir.exists(DATA_DERIVED_DIR)) dir.create(DATA_DERIVED_DIR, recursive = TRUE)

suppressPackageStartupMessages({
  library(BioAge)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(flexsurv)
})

has_corrplot <- requireNamespace("corrplot", quietly = TRUE)
if (has_corrplot) library(corrplot)

# Source modules
for (mf in c("R/config.R", "R/train.R", "R/assemble.R", "R/residualize.R",
             "R/export.R", "R/qc_recenter.R")) {
  source(file.path(PROJECT_ROOT, mf), local = FALSE)
}
message(">> All modules loaded")

# Output directory
OUTPUT_DIR <- file.path(PROJECT_ROOT, "final_output")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Focus clocks — the ones with meaningful Argentine coverage
FOCUS_CLOCKS <- c("pheno_max_v2")

# ══════════════════════════════════════════════════════════════════════════
# STEP 1. TRAIN CLOCKS ON NHANES III → VALIDATE ON NHANES IV
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 1: TRAIN CLOCKS ON NHANES                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

message(">> Training all clocks (NHANES III → NHANES IV)...")

clocks <- train_all_clocks(
  panels_sub = PANELS_SUB,
  train_kdm  = TRUE,
  train_hd   = TRUE,
  verbose    = TRUE
)

# Assemble NHANES IV dataset
data_nhanes <- assemble_data(clocks)

# ── 1a. NHANES IV Performance Report ──

cat("\n--- NHANES IV Clock Performance ---\n\n")

nhanes_perf <- list()

# BA columns for all trained clocks
ba_cols_nhanes <- list_ba_columns(data_nhanes)

for (bc in ba_cols_nhanes) {
  vals <- data_nhanes[[bc]]
  ages <- data_nhanes$age
  ok   <- !is.na(vals) & !is.na(ages) & is.finite(vals)
  n_ok <- sum(ok)
  if (n_ok < 30) next

  ct     <- cor.test(ages[ok], vals[ok])
  lm_fit <- lm(vals[ok] ~ ages[ok])
  s      <- summary(lm_fit)

  nhanes_perf[[bc]] <- data.frame(
    clock     = bc,
    n         = n_ok,
    r         = round(ct$estimate, 4),
    r_sq      = round(s$r.squared, 4),
    slope     = round(coef(lm_fit)[2], 4),
    intercept = round(coef(lm_fit)[1], 2),
    RMSE      = round(sqrt(mean(s$residuals^2)), 2),
    mean_adv  = round(mean(vals[ok] - ages[ok]), 2),
    sd_adv    = round(sd(vals[ok] - ages[ok]), 2),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

if (length(nhanes_perf) > 0) {
  nhanes_table <- do.call(rbind, nhanes_perf)
  rownames(nhanes_table) <- NULL
  print(nhanes_table, row.names = FALSE)
  write.csv(nhanes_table, file.path(OUTPUT_DIR, "step1_nhanes_performance.csv"),
            row.names = FALSE)
  message(">> Saved: step1_nhanes_performance.csv")
}

# ── 1b. NHANES IV scatter plots (BA vs CA) for focus clocks ──

cat("\n--- NHANES IV: BA vs CA Scatter Plots ---\n\n")

focus_ba <- paste0("pheno_", FOCUS_CLOCKS)
focus_ba <- intersect(focus_ba, ba_cols_nhanes)

for (bc in focus_ba) {
  vals <- data_nhanes[[bc]]
  ok   <- !is.na(vals) & is.finite(vals)
  if (sum(ok) < 30) next

  r_val <- cor(data_nhanes$age[ok], vals[ok])
  label <- gsub("pheno_", "", bc)

  p <- ggplot(data_nhanes[ok, ], aes(x = age, y = .data[[bc]])) +
    geom_point(alpha = 0.08, size = 0.5, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("NHANES IV: BA vs CA — ", label),
         subtitle = sprintf("r = %.4f, n = %s", r_val, format(sum(ok), big.mark = ",")),
         x = "Chronological Age (years)",
         y = "Biological Age (PhenoAge)")
  print(p)

  ggsave(file.path(OUTPUT_DIR, paste0("step1_nhanes_scatter_", label, ".png")),
         p, width = 8, height = 6, dpi = 150)
}

# ── 1c. Save deployment bundle ──

# Z-score and residualize for bundle export
ba_all <- list_ba_columns(data_nhanes)
data_nhanes <- add_zscores(data_nhanes, ba_all)
adv_cols_nhanes <- list_advance_columns(data_nhanes)
data_nhanes <- add_residuals(data_nhanes, adv_cols_nhanes)

bundle <- build_deployment_bundle(
  clocks     = clocks,
  data       = data_nhanes,
  clock_vars = adv_cols_nhanes,
  hr_results = NULL,
  qc         = NULL
)
save_bundle(bundle, file.path(DATA_DERIVED_DIR, "bioage_deployment_bundle.rds"))

message("\n>> STEP 1 COMPLETE: ", length(bundle$meta$clocks_trained), " clocks trained")