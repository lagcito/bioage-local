# run_score_argentina.R
# Script to score Argentine clinical cohort with focus clocks and evaluate performance.
# ══════════════════════════════════════════════════════════════════════════
# STEP 2. SCORE ARGENTINE CLINICAL COHORT
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 2: SCORE ARGENTINE COHORT                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

message(">> Loading clinical data...")
f1 <- file.path(PROJECT_ROOT, "data_raw", "phenoage_master_PARTE1.csv")
f2 <- file.path(PROJECT_ROOT, "data_raw", "phenoage_master_PARTE2.csv")
cohort <- read.csv(f1, stringsAsFactors = FALSE)
if (file.exists(f2)) cohort <- rbind(cohort, read.csv(f2, stringsAsFactors = FALSE))
message(">> Rows: ", nrow(cohort))

# ── 2a. Transform biomarkers ──

df <- cohort %>%
  mutate(
    patient_id  = Protocolo,
    age         = as.numeric(age),
    albumin     = as.numeric(albumin) * 10,         # g/dL → g/L
    alp         = as.numeric(alp),
    glucose_mg  = as.numeric(glucose),              # keep mg/dL for HOMA-IR
    glucose     = glucose_mg / 18.0,                # mg/dL → mmol/L
    crp_raw     = as.numeric(crp),
    crp         = crp_raw / 10,                     # mg/L → mg/dL
    lncrp       = ifelse(is.na(crp) | crp <= 0, NA, log(crp)),
    creat       = as.numeric(creatinine) * 88.4,    # mg/dL → μmol/L
    lncreat     = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),
    lymph       = as.numeric(lymphocyte),
    mcv         = as.numeric(mcv),
    rdw         = as.numeric(rdw),
    wbc         = ifelse(as.numeric(wbc) > 300, as.numeric(wbc) / 1000, as.numeric(wbc)),
    rbc         = as.numeric(rbc),
    rbc         = ifelse(rbc > 100, rbc / 1e6, rbc),
    ggt         = as.numeric(ggt),
    insulin     = as.numeric(insulin),
    trig        = as.numeric(triglycerides),
    totchol     = as.numeric(cholesterol),
    hba1c       = as.numeric(hba1c),
    vitaminB12  = as.numeric(vitamin_b12),
    bun         = as.numeric(bun),
    uap         = as.numeric(uric_acid),
    gender      = ifelse(sex == "M", 1, 2),
    lnglucose   = ifelse(is.na(glucose) | glucose <= 0, NA, log(glucose))
  )

df <- df[!is.na(df$age) & df$age >= 18 & df$age <= 120, ]
message(">> After age filter: ", nrow(df), " rows")
message(">> Age: ", round(mean(df$age, na.rm = TRUE), 1), " ± ",
        round(sd(df$age, na.rm = TRUE), 1), " years")

# ── 2a-bis. Derived biomarkers: NLR and HOMA-IR ──

# NLR = Neutrophil-to-Lymphocyte Ratio (works with % or absolute)
df$nlr <- ifelse(!is.na(df$neutrophils) & !is.na(df$lymphocyte) &
                 df$lymphocyte > 0,
                 as.numeric(df$neutrophils) / as.numeric(df$lymphocyte), NA)

# HOMA-IR = (glucose_mg_dL × insulin_μU_mL) / 405
# glucose_mg was preserved before unit conversion in the mutate above
df$homa_ir <- ifelse(!is.na(df$glucose_mg) & !is.na(df$insulin) &
                     df$glucose_mg > 0 & df$insulin > 0,
                     (df$glucose_mg * df$insulin) / 405, NA)

n_nlr  <- sum(!is.na(df$nlr))
n_homa <- sum(!is.na(df$homa_ir))
message(">> NLR:     ", format(n_nlr, big.mark = ","), " valid values (",
        round(100 * n_nlr / nrow(df), 1), "%)")
message(">> HOMA-IR: ", format(n_homa, big.mark = ","), " valid values (",
        round(100 * n_homa / nrow(df), 1), "%)")

# ── 2a-ter. QC: Bandera de posible patología hematológica ──
#
# Pacientes con leucocitosis extrema o linfocitosis marcada pueden corresponder
# a leucemias o linfomas (LLC, LMC, linfoma leucemizado).
# Su inclusión puede sesgar artificialmente los perfiles de envejecimiento,
# especialmente el efecto de linfocitos_abs y WBC.
#
# Criterios de flagueado (conservadores, OR entre criterios):
#   - WBC > 30 ×10³/µL   → leucocitosis severa (sugestivo de LMC/LLC)
#   - linfocitos% > 70%  → linfocitosis relativa marcada
#   - linfocitos_abs > 10 ×10³/µL → linfocitosis absoluta clínica
#
# Los pacientes son MARCADOS pero NO EXCLUIDOS del scoring.
# El análisis de constructo externo se repite excluyéndolos (análisis
# de sensibilidad) para evaluar si el resultado de linfocitos cambia.

cat("\n--- QC: Posible Patología Hematológica ---\n\n")

wbc_num  <- suppressWarnings(as.numeric(df$wbc))
lymp_num <- suppressWarnings(as.numeric(df$lymphocyte))
lymabs_num <- if ("lymphocyte_abs" %in% names(df))
                suppressWarnings(as.numeric(df$lymphocyte_abs)) else
                rep(NA_real_, nrow(df))

flag_wbc    <- !is.na(wbc_num)    & wbc_num    > 30
flag_lymp   <- !is.na(lymp_num)   & lymp_num   > 70
flag_lymabs <- !is.na(lymabs_num) & lymabs_num > 10

df$flag_hematol <- flag_wbc | flag_lymp | flag_lymabs

n_fw  <- sum(flag_wbc,    na.rm = TRUE)
n_fl  <- sum(flag_lymp,   na.rm = TRUE)
n_fla <- sum(flag_lymabs, na.rm = TRUE)
n_ft  <- sum(df$flag_hematol, na.rm = TRUE)

cat(sprintf("  WBC > 30 ×10³/µL         : %d pacientes\n", n_fw))
cat(sprintf("  Linfocitos%% > 70%%         : %d pacientes\n", n_fl))
cat(sprintf("  Linfocitos abs > 10 ×10³  : %d pacientes\n", n_fla))
cat(sprintf("  Total flagueados (OR)     : %d (%.2f%% de la cohorte)\n",
            n_ft, 100 * n_ft / nrow(df)))

if (n_ft > 0) {
  fl_idx <- df$flag_hematol & !is.na(df$flag_hematol)
  cat(sprintf("\n  WBC en flagueados        : mediana=%.1f  máx=%.1f ×10³/µL\n",
              median(wbc_num[fl_idx], na.rm = TRUE),
              max(wbc_num[fl_idx],    na.rm = TRUE)))
  cat(sprintf("  Linfocitos%% en flagueados: mediana=%.1f  máx=%.1f\n",
              median(lymp_num[fl_idx], na.rm = TRUE),
              max(lymp_num[fl_idx],    na.rm = TRUE)))
  if (any(!is.na(lymabs_num[fl_idx])))
    cat(sprintf("  Linf abs en flagueados   : mediana=%.2f  máx=%.2f ×10³/µL\n",
                median(lymabs_num[fl_idx], na.rm = TRUE),
                max(lymabs_num[fl_idx],    na.rm = TRUE)))
}

# ── 2b. Biomarker availability ──

cat("\n--- Biomarker Availability ---\n\n")

all_bm <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
            "lymph", "mcv", "rdw", "wbc", "rbc", "ggt", "insulin",
            "trig", "totchol", "vitaminB12", "hba1c", "bun", "uap", "lnglucose")
for (v in all_bm) {
  if (v %in% names(df)) {
    n_ok <- sum(!is.na(df[[v]]))
    cat(sprintf("  %-14s %6d / %d  (%5.1f%%)\n",
                v, n_ok, nrow(df), 100 * n_ok / nrow(df)))
  }
}

# ── 2c. Score focus clocks ──

cat("\n--- Scoring Focus Clocks ---\n\n")

clocks_to_score <- list()

# Add focus sub-clocks from bundle
for (nm in names(bundle$subclock_models)) {
  if (nm == "pheno_max_v2") {
    clocks_to_score[[paste0("pheno_", nm)]] <- list(
      biomarkers = bundle$panels$subclocks[[nm]],
      fit        = bundle$subclock_models[[nm]]
    )
  }
}

# Master results
results <- data.frame(
  patient_id   = df$patient_id,
  age          = df$age,
  sex          = df$sex,
  nlr          = df$nlr,
  homa_ir      = df$homa_ir,
  flag_hematol = df$flag_hematol,
  stringsAsFactors = FALSE
)

score_summary <- list()

for (clock_name in names(clocks_to_score)) {
  clock <- clocks_to_score[[clock_name]]
  bm    <- clock$biomarkers
  fit   <- clock$fit

  bm_plus_age <- c(bm, "age")
  missing_cols <- setdiff(bm_plus_age, names(df))
  if (length(missing_cols) > 0) next

  complete_mask <- complete.cases(df[, bm_plus_age, drop = FALSE])
  n_complete <- sum(complete_mask)
  if (n_complete == 0) {
    message(">> ", clock_name, ": SKIP — 0 complete cases")
    next
  }

  message(">> ", clock_name, ": ", n_complete, " complete cases (",
          round(100 * n_complete / nrow(df), 1), "%)")

  dat_subset <- df[complete_mask, bm_plus_age, drop = FALSE]

  calc_result <- tryCatch({
    phenoage_calc(data = dat_subset, biomarkers = bm, fit = fit)
  }, error = function(e) {
    message("   ERROR: ", e$message)
    NULL
  })
  if (is.null(calc_result)) next

  pa_vals <- calc_result$data$phenoage
  n_valid <- sum(!is.na(pa_vals) & is.finite(pa_vals))
  message("   -> ", n_valid, " valid PhenoAge values")

  if (n_valid > 0) {
    ba_col  <- clock_name
    adv_col <- paste0(clock_name, "_advance")
    results[[ba_col]]  <- NA_real_
    results[[adv_col]] <- NA_real_
    results[[ba_col]][complete_mask]  <- pa_vals
    results[[adv_col]][complete_mask] <- calc_result$data$phenoage_advance

    score_summary[[clock_name]] <- data.frame(
      clock    = clock_name,
      n_scored = n_valid,
      pct      = round(100 * n_valid / nrow(df), 1),
      stringsAsFactors = FALSE
    )
  }
}

# ── 2c-bis. Missing data characterisation: scored vs unscored ──
#
# The 43% of patients without a score may differ systematically from those scored.
# A simple age + sex comparison flags potential selection bias.

cat("\n--- Missing Data: Scored vs Unscored Comparison ---\n\n")

scored_mask <- !is.na(results$pheno_pheno_max_v2) &
               is.finite(results$pheno_pheno_max_v2)

scored_age   <- results$age[scored_mask]
unscored_age <- results$age[!scored_mask]
scored_sex   <- results$sex[scored_mask]
unscored_sex <- results$sex[!scored_mask]

n_scored   <- sum(scored_mask)
n_unscored <- sum(!scored_mask)

cat(sprintf("  Scored:   n = %s  (%.1f%%)\n",
            format(n_scored,   big.mark = ","), 100 * n_scored / nrow(results)))
cat(sprintf("  Unscored: n = %s  (%.1f%%)\n",
            format(n_unscored, big.mark = ","), 100 * n_unscored / nrow(results)))

# Age comparison
age_t <- tryCatch(t.test(scored_age, unscored_age), error = function(e) NULL)
cat(sprintf("\n  Age — Scored: %.1f ± %.1f  |  Unscored: %.1f ± %.1f",
            mean(scored_age, na.rm = TRUE), sd(scored_age, na.rm = TRUE),
            mean(unscored_age, na.rm = TRUE), sd(unscored_age, na.rm = TRUE)))
if (!is.null(age_t)) {
  smd_age <- (mean(scored_age, na.rm = TRUE) - mean(unscored_age, na.rm = TRUE)) /
             sqrt((var(scored_age, na.rm = TRUE) + var(unscored_age, na.rm = TRUE)) / 2)
  flag_age <- if (abs(smd_age) > 0.1) "  FLAG (SMD > 0.1)" else "  OK"
  cat(sprintf("  SMD = %.3f  p = %.4f%s\n", smd_age, age_t$p.value, flag_age))
} else cat("\n")

# Sex comparison
if (all(c("F", "M") %in% c(scored_sex, unscored_sex))) {
  pct_F_scored   <- 100 * mean(scored_sex   == "F", na.rm = TRUE)
  pct_F_unscored <- 100 * mean(unscored_sex == "F", na.rm = TRUE)
  sex_chisq <- tryCatch(
    chisq.test(rbind(table(scored_sex), table(unscored_sex))),
    error = function(e) NULL
  )
  cat(sprintf("  Sex %%F — Scored: %.1f%%  |  Unscored: %.1f%%",
              pct_F_scored, pct_F_unscored))
  if (!is.null(sex_chisq)) {
    flag_sex <- if (sex_chisq$p.value < 0.05) "  FLAG (p < 0.05)" else "  OK"
    cat(sprintf("  p = %.4f%s\n", sex_chisq$p.value, flag_sex))
  } else cat("\n")
}

# Save summary
missing_summary <- data.frame(
  group        = c("scored", "unscored"),
  n            = c(n_scored, n_unscored),
  pct          = round(100 * c(n_scored, n_unscored) / nrow(results), 1),
  mean_age     = round(c(mean(scored_age, na.rm = TRUE),
                         mean(unscored_age, na.rm = TRUE)), 1),
  sd_age       = round(c(sd(scored_age, na.rm = TRUE),
                         sd(unscored_age, na.rm = TRUE)), 1),
  pct_female   = round(c(100 * mean(scored_sex   == "F", na.rm = TRUE),
                         100 * mean(unscored_sex == "F", na.rm = TRUE)), 1),
  stringsAsFactors = FALSE
)
write.csv(missing_summary,
          file.path(OUTPUT_DIR, "step2_missing_data_comparison.csv"),
          row.names = FALSE)
message(">> Saved: step2_missing_data_comparison.csv")

# ── 2d. Argentine scoring performance ──

cat("\n--- Argentine Cohort: Clock Performance ---\n\n")

ba_cols_arg <- names(score_summary)
arg_perf <- list()

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ages <- results$age
  ok   <- !is.na(vals) & !is.na(ages) & is.finite(vals)
  n_ok <- sum(ok)
  if (n_ok < 10) next

  ct     <- cor.test(ages[ok], vals[ok])
  lm_fit <- lm(vals[ok] ~ ages[ok])
  s      <- summary(lm_fit)

  arg_perf[[bc]] <- data.frame(
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

  cat(sprintf("  %-40s r = %.4f  RMSE = %.2f  n = %s\n",
              bc, ct$estimate, sqrt(mean(s$residuals^2)),
              format(n_ok, big.mark = ",")))
}

if (length(arg_perf) > 0) {
  arg_table <- do.call(rbind, arg_perf)
  rownames(arg_table) <- NULL
  cat("\n")
  print(arg_table, row.names = FALSE)
  write.csv(arg_table, file.path(OUTPUT_DIR, "step2_argentina_performance.csv"),
            row.names = FALSE)
  message(">> Saved: step2_argentina_performance.csv")
}

# ── 2d-bis. RMSE by Age Group (Argentine Cohort) ──

cat("\n--- RMSE by Age Group (Argentine Cohort) ---\n\n")

age_breaks  <- c(0, 30, 40, 50, 60, 70, 80, 120)
age_labels  <- c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
results$age_group_rmse <- cut(results$age, breaks = age_breaks, labels = age_labels,
                              right = FALSE)

rmse_by_age <- list()
for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ages <- results$age
  grp  <- results$age_group_rmse
  for (ag in age_labels) {
    ok <- !is.na(vals) & is.finite(vals) & !is.na(ages) & !is.na(grp) & grp == ag
    n_ok <- sum(ok)
    if (n_ok < 10) next
    resid <- vals[ok] - ages[ok]
    rmse_val  <- sqrt(mean(resid^2))
    mae_val   <- mean(abs(resid))
    r_val     <- cor(ages[ok], vals[ok])
    bias_val  <- mean(resid)
    rmse_by_age[[paste0(bc, "_", ag)]] <- data.frame(
      clock     = bc,
      age_group = ag,
      n         = n_ok,
      RMSE      = round(rmse_val, 2),
      MAE       = round(mae_val, 2),
      bias      = round(bias_val, 2),
      r         = round(r_val, 4),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
}

if (length(rmse_by_age) > 0) {
  rmse_age_table <- do.call(rbind, rmse_by_age)
  rownames(rmse_age_table) <- NULL

  # Print per clock
  for (bc in ba_cols_arg) {
    sub_t <- rmse_age_table[rmse_age_table$clock == bc, ]
    if (nrow(sub_t) == 0) next
    cat(sprintf("\n  %s:\n", bc))
    cat(sprintf("    %-8s %6s %7s %7s %7s %7s\n", "Age", "n", "RMSE", "MAE", "Bias", "r"))
    for (i in seq_len(nrow(sub_t))) {
      r <- sub_t[i, ]
      cat(sprintf("    %-8s %6d %7.2f %7.2f %7.2f %7.4f\n",
                  r$age_group, r$n, r$RMSE, r$MAE, r$bias, r$r))
    }
    # Highlight worst RMSE group
    worst <- sub_t[which.max(sub_t$RMSE), ]
    best  <- sub_t[which.min(sub_t$RMSE), ]
    cat(sprintf("    >> Worst RMSE: %s (%.2f)  |  Best RMSE: %s (%.2f)\n",
                worst$age_group, worst$RMSE, best$age_group, best$RMSE))
  }

  write.csv(rmse_age_table, file.path(OUTPUT_DIR, "step2_rmse_by_age_group.csv"),
            row.names = FALSE)
  message(">> Saved: step2_rmse_by_age_group.csv")

  # Plot: RMSE by age group per clock
  for (bc in unique(rmse_age_table$clock)) {
    sub_t <- rmse_age_table[rmse_age_table$clock == bc, ]
    if (nrow(sub_t) < 3) next
    label <- gsub("pheno_", "", bc)
    sub_t$age_group <- factor(sub_t$age_group, levels = age_labels)

    p_rmse <- ggplot(sub_t, aes(x = age_group, y = RMSE)) +
      geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
      geom_text(aes(label = sprintf("%.1f", RMSE)), vjust = -0.5, size = 3.5) +
      geom_line(aes(x = as.numeric(age_group), y = RMSE),
                colour = "darkred", linewidth = 0.8, linetype = "dashed") +
      theme_minimal(base_size = 13) +
      labs(title = paste0("RMSE by Age Group — ", label),
           subtitle = sprintf("n per group shown below bars"),
           x = "Age Group", y = "RMSE (years)") +
      geom_text(aes(label = paste0("n=", format(n, big.mark = ",")), y = 0),
                vjust = 1.5, size = 2.8, colour = "grey40")
    print(p_rmse)

    fn <- file.path(OUTPUT_DIR, paste0("step2_rmse_by_age_", label, ".png"))
    ggsave(fn, p_rmse, width = 8, height = 5, dpi = 150)
  }
  message(">> Saved: RMSE by age group plots")
}

# ── 2d-ter. RMSE estratificado: por sexo, por quintil BAA, global con bootstrap CI ──

cat("\n--- RMSE Stratified Analysis ---\n\n")

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ages <- results$age
  ok   <- !is.na(vals) & is.finite(vals) & !is.na(ages)
  if (sum(ok) < 30) next
  label <- gsub("pheno_", "", bc)
  adv_col <- paste0(bc, "_advance")

  cat(sprintf("\n  %s\n", label))
  cat(sprintf("  %s\n", strrep("-", nchar(label) + 2)))

  # (a) Global RMSE with 1000-resample bootstrap 95% CI
  resid_global <- vals[ok] - ages[ok]
  rmse_global  <- sqrt(mean(resid_global^2))
  set.seed(42)
  boot_rmse <- replicate(1000, {
    s <- sample(resid_global, length(resid_global), replace = TRUE)
    sqrt(mean(s^2))
  })
  ci_lo <- quantile(boot_rmse, 0.025)
  ci_hi <- quantile(boot_rmse, 0.975)
  cat(sprintf("    Global RMSE: %.2f  [95%% CI: %.2f – %.2f]  n = %s\n",
              rmse_global, ci_lo, ci_hi, format(sum(ok), big.mark = ",")))

  # (b) By sex
  cat(sprintf("    By sex:\n"))
  for (sx in c("F", "M")) {
    idx <- ok & !is.na(results$sex) & results$sex == sx
    if (sum(idx) < 10) next
    r_sx <- vals[idx] - ages[idx]
    cat(sprintf("      %s: RMSE = %.2f  MAE = %.2f  bias = %+.2f  n = %d\n",
                sx, sqrt(mean(r_sx^2)), mean(abs(r_sx)), mean(r_sx), sum(idx)))
  }

  # (c) By BAA quintile (P1–P5) — uses advancement column if available
  if (adv_col %in% names(results)) {
    adv_vals <- results[[adv_col]]
    ok_adv   <- ok & !is.na(adv_vals) & is.finite(adv_vals)
    if (sum(ok_adv) >= 50) {
      qts_adv <- quantile(adv_vals[ok_adv], probs = c(0.10, 0.25, 0.75, 0.90))
      results$baa_quintile <- NA_character_
      results$baa_quintile[ok_adv] <- as.character(cut(
        adv_vals[ok_adv],
        breaks = c(-Inf, qts_adv[1], qts_adv[2], qts_adv[3], qts_adv[4], Inf),
        labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
      ))
      cat(sprintf("    By BAA quintile:\n"))
      for (grp in c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")) {
        idx_g <- ok & !is.na(results$baa_quintile) & results$baa_quintile == grp
        if (sum(idx_g) < 10) next
        r_g <- vals[idx_g] - ages[idx_g]
        cat(sprintf("      %-15s RMSE = %.2f  MAE = %.2f  bias = %+.2f  n = %d\n",
                    grp, sqrt(mean(r_g^2)), mean(abs(r_g)), mean(r_g), sum(idx_g)))
      }
    }
  }
}

# ── 2e. Argentine scatter plots (BA vs CA) ──

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ok   <- !is.na(vals) & is.finite(vals)
  if (sum(ok) < 30) next

  r_val <- cor(results$age[ok], vals[ok])
  label <- gsub("pheno_", "", bc)

  p <- ggplot(results[ok, ], aes(x = age, y = .data[[bc]])) +
    geom_point(alpha = 0.05, size = 0.3, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "red", linewidth = 0.7) +
    geom_smooth(method = "lm", colour = "navy", se = TRUE, linewidth = 1) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Argentina: BA vs CA — ", label),
         subtitle = sprintf("r = %.4f, n = %s", r_val, format(sum(ok), big.mark = ",")),
         x = "Chronological Age (years)",
         y = "Biological Age (PhenoAge)")
  print(p)

  ggsave(file.path(OUTPUT_DIR, paste0("step2_arg_scatter_", label, ".png")),
         p, width = 8, height = 6, dpi = 150)
}

# ── 2f. BA vs CA by sex ──

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ok <- !is.na(vals) & is.finite(vals) & results$sex %in% c("F", "M")
  if (sum(ok) < 50) next
  label <- gsub("pheno_", "", bc)

  p_sex <- ggplot(results[ok, ], aes(x = age, y = .data[[bc]], colour = sex)) +
    geom_point(alpha = 0.03, size = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
    scale_colour_manual(values = c("F" = "firebrick", "M" = "steelblue")) +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Argentina — BA vs CA by Sex: ", label),
         x = "Chronological Age", y = "Biological Age", colour = "Sex")
  print(p_sex)

  ggsave(file.path(OUTPUT_DIR, paste0("step2_arg_sex_", label, ".png")),
         p_sex, width = 8, height = 6, dpi = 150)
}

# ── 2g. Bland-Altman: BA vs CA agreement ──

cat("\n--- Bland-Altman Agreement Analysis ---\n\n")

ba_summary_rows <- list()

for (bc in ba_cols_arg) {
  vals <- results[[bc]]
  ok   <- !is.na(vals) & is.finite(vals) & !is.na(results$age)
  if (sum(ok) < 30) next
  label <- gsub("pheno_", "", bc)

  mean_ba_ca <- (vals[ok] + results$age[ok]) / 2
  diff_ba_ca <- vals[ok] - results$age[ok]
  mean_diff  <- mean(diff_ba_ca)
  sd_diff    <- sd(diff_ba_ca)
  loa_lo     <- mean_diff - 1.96 * sd_diff
  loa_hi     <- mean_diff + 1.96 * sd_diff
  pct_within <- 100 * mean(abs(diff_ba_ca) <= 5)

  # LoA empíricos (percentiles 2.5/97.5 de la diferencia BA-CA), sin asumir
  # normalidad — más apropiados que los paramétricos cuando la distribución
  # es asimétrica (ver sección de distribución de Z, típicamente right-skewed
  # en este reloj). Se reportan ambos para que se pueda ver la divergencia.
  loa_emp_lo <- as.numeric(quantile(diff_ba_ca, 0.025))
  loa_emp_hi <- as.numeric(quantile(diff_ba_ca, 0.975))

  cat(sprintf("  %s  |  Bias: %+.2f yr  |  SD: %.2f  |  LoA parametrico: [%+.2f, %+.2f]  |  LoA empirico: [%+.2f, %+.2f]  |  within ±5yr: %.1f%%\n",
              label, mean_diff, sd_diff, loa_lo, loa_hi, loa_emp_lo, loa_emp_hi, pct_within))

  ba_summary_rows[[bc]] <- data.frame(
    clock          = bc,
    n              = sum(ok),
    bias           = round(mean_diff, 3),
    sd_diff        = round(sd_diff, 3),
    loa_lo         = round(loa_lo, 3),
    loa_hi         = round(loa_hi, 3),
    loa_lo_empirical = round(loa_emp_lo, 3),
    loa_hi_empirical = round(loa_emp_hi, 3),
    pct_within_5yr = round(pct_within, 1),
    stringsAsFactors = FALSE
  )

  # Plot
  ba_df <- data.frame(mean_ba_ca = mean_ba_ca, diff_ba_ca = diff_ba_ca)

  p_ba <- ggplot(ba_df, aes(x = mean_ba_ca, y = diff_ba_ca)) +
    geom_point(alpha = 0.04, size = 0.3, colour = "steelblue") +
    geom_hline(yintercept = mean_diff,  colour = "black",   linewidth = 0.8) +
    geom_hline(yintercept = loa_lo,     colour = "firebrick", linewidth = 0.7, linetype = "dashed") +
    geom_hline(yintercept = loa_hi,     colour = "firebrick", linewidth = 0.7, linetype = "dashed") +
    geom_hline(yintercept = loa_emp_lo, colour = "darkorange", linewidth = 0.7, linetype = "dotdash") +
    geom_hline(yintercept = loa_emp_hi, colour = "darkorange", linewidth = 0.7, linetype = "dotdash") +
    geom_smooth(method = "loess", colour = "darkgreen", se = FALSE,
                linewidth = 0.7, linetype = "dotted") +
    annotate("text", x = Inf, y = mean_diff, hjust = 1.1, vjust = -0.5, size = 3,
             label = sprintf("Bias = %+.2f", mean_diff)) +
    annotate("text", x = Inf, y = loa_hi, hjust = 1.1, vjust = -0.5, size = 3,
             colour = "firebrick", label = sprintf("LoA parametrico = %+.2f", loa_hi)) +
    annotate("text", x = Inf, y = loa_lo, hjust = 1.1, vjust = 1.5, size = 3,
             colour = "firebrick", label = sprintf("LoA parametrico = %+.2f", loa_lo)) +
    annotate("text", x = Inf, y = loa_emp_hi, hjust = 1.1, vjust = 1.5, size = 3,
             colour = "darkorange", label = sprintf("LoA empirico = %+.2f", loa_emp_hi)) +
    annotate("text", x = Inf, y = loa_emp_lo, hjust = 1.1, vjust = -0.5, size = 3,
             colour = "darkorange", label = sprintf("LoA empirico = %+.2f", loa_emp_lo)) +
    theme_minimal(base_size = 13) +
    labs(title    = paste0("Bland-Altman: BA vs CA — ", label),
         subtitle = sprintf("Bias = %+.2f yr | LoA parametrico [%+.2f, %+.2f] | LoA empirico [%+.2f, %+.2f] | n = %s",
                            mean_diff, loa_lo, loa_hi, loa_emp_lo, loa_emp_hi, format(sum(ok), big.mark = ",")),
         x = "Mean of BA and CA (years)",
         y = "BA − CA (years)")
  print(p_ba)

  ggsave(file.path(OUTPUT_DIR, paste0("step2_bland_altman_", label, ".png")),
         p_ba, width = 8, height = 6, dpi = 150)
}

if (length(ba_summary_rows) > 0) {
  ba_summary_df <- do.call(rbind, ba_summary_rows)
  rownames(ba_summary_df) <- NULL
  write.csv(ba_summary_df,
            file.path(OUTPUT_DIR, "step2_bland_altman_summary.csv"),
            row.names = FALSE)
  message(">> Saved: step2_bland_altman_summary.csv")
}

# Save scored data
write.csv(results, file.path(DATA_DERIVED_DIR, "clinical_scored_lite.csv"),
          row.names = FALSE)
message(">> Saved: clinical_scored_lite.csv")

message("\n>> STEP 2 COMPLETE: ", length(score_summary), " clocks scored on ",
        format(nrow(results), big.mark = ","), " patients")
