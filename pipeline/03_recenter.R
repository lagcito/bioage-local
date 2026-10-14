#run_recentering.R
# ══════════════════════════════════════════════════════════════════════════
# STEP 3. POPULATION RECENTERING
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 3: POPULATION RECENTERING                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Work on a copy of results for recentering
rc <- results

# Identify advancement columns with enough data
adv_cols <- grep("_advance$", names(rc), value = TRUE)
viable_adv <- character(0)
for (ac in adv_cols) {
  n_ok <- sum(!is.na(rc[[ac]]) & is.finite(rc[[ac]]))
  if (n_ok >= 100) viable_adv <- c(viable_adv, ac)
}
message(">> Viable clocks for recentering: ", paste(sub("_advance$", "", viable_adv),
                                                     collapse = ", "))

# ── 3b. Method B: Sex-adjusted residualization (BA ~ age + sex) ──
#
# Fit lm(BA ~ age + sex) once per clock; reuse fit for residuals, Z-scores,
# and exported coefficients (avoids fitting the same model twice).

cat("\n--- Method B: Sex-Adjusted Residualization (BA ~ age + sex) ---\n\n")

arg_method_b_coefs <- list()

for (ac in viable_adv) {
  ba_col  <- sub("_advance$", "", ac)
  ba_vals <- rc[[ba_col]]
  ok <- !is.na(ba_vals) & is.finite(ba_vals) &
        !is.na(rc$age) & rc$sex %in% c("F", "M")
  if (sum(ok) < 30) next

  fit_sex    <- lm(ba_vals[ok] ~ rc$age[ok] + rc$sex[ok])
  coefs      <- coef(fit_sex)
  resid_vals <- residuals(fit_sex)

  resid_col <- paste0(ba_col, "_sexadj_resid")
  z_col     <- paste0(ba_col, "_sexadj_z")

  rc[[resid_col]] <- NA_real_
  rc[[resid_col]][ok] <- resid_vals
  rc[[z_col]] <- NA_real_
  valid_resid <- !is.na(rc[[resid_col]])
  if (sum(valid_resid) > 1)
    rc[[z_col]][valid_resid] <- as.numeric(scale(rc[[resid_col]][valid_resid]))

  # ── Eliminación de outliers extremos en Z-score (|Z| > 5) ──
  #
  # Valores |Z| > 5 corresponden a patología severa o errores analíticos.
  # Con |Z| > 5: BA > CA + 5*SD_resid → biológicamente implausible en ausencia
  # de enfermedad grave. Se anulan para no distorsionar los análisis de
  # construcción de perfiles y los cortes de percentiles normativos.
  # Se imprimen los IDs para auditoría.

  z_extreme_idx <- !is.na(rc[[z_col]]) & abs(rc[[z_col]]) > 5
  n_extreme <- sum(z_extreme_idx)
  if (n_extreme > 0) {
    cat(sprintf("\n    [OUTLIERS] |Z| > 5 en %s: %d caso(s)\n", ba_col, n_extreme))

    # ── Cruce con flag_hematol: ¿cuánto de lo censurado ya tenía una causa
    # hematológica identificada, y cuánto queda sin explicación clínica? ──
    # Sin este cruce, censurar |Z|>5 podría estar descartando señal real
    # (patología no hematológica: oncológica, hepática aguda, sepsis) en vez
    # de sólo ruido/error analítico.
    if ("flag_hematol" %in% names(rc)) {
      fh <- rc$flag_hematol[z_extreme_idx]
      n_explained   <- sum(fh, na.rm = TRUE)
      n_unexplained <- sum(!fh | is.na(fh))
      cat(sprintf("    Explicados por flag_hematol : %d (%.0f%%)\n",
                  n_explained, 100 * n_explained / n_extreme))
      cat(sprintf("    SIN explicación hematológica: %d (%.0f%%)  — posible señal clínica real descartada\n",
                  n_unexplained, 100 * n_unexplained / n_extreme))
    }

    if (n_extreme <= 30) {
      sub_ext <- rc[z_extreme_idx, ]
      for (ie in seq_len(nrow(sub_ext))) {
        cat(sprintf("      ID=%-10s  Z=%+7.3f  BA=%5.1f yr  CA=%5.1f yr  sex=%s  hematol=%s\n",
                    sub_ext$patient_id[ie],
                    rc[[z_col]][z_extreme_idx][ie],
                    sub_ext[[ba_col]][ie],
                    sub_ext$age[ie],
                    sub_ext$sex[ie],
                    if ("flag_hematol" %in% names(sub_ext)) isTRUE(sub_ext$flag_hematol[ie]) else NA))
      }
    }

    # Guardar el detalle completo (no sólo los primeros 30) para auditoría.
    # A NIVEL DE PACIENTE (patient_id, edad, sexo) -> va a data_derived/, NUNCA
    # a final_output/, porque final_output/ se versiona en git y este detalle
    # es identificable. Ver hallazgo de auditoría sobre datos de pacientes.
    out_cols <- intersect(c("patient_id", "age", "sex", ba_col, "flag_hematol"), names(rc))
    extreme_log <- rc[z_extreme_idx, out_cols, drop = FALSE]
    extreme_log$z_value <- rc[[z_col]][z_extreme_idx]
    extreme_log$clock <- ba_col
    write.csv(extreme_log,
              file.path(DATA_DERIVED_DIR, paste0("step3_zoutliers_censored_", ba_col, ".csv")),
              row.names = FALSE)
    message(">> Guardado detalle de outliers censurados (data_derived/, no versionado): step3_zoutliers_censored_", ba_col, ".csv")

    # Resumen AGREGADO (sin patient_id) — este sí es seguro para final_output/
    write.csv(
      data.frame(
        clock                 = ba_col,
        n_censurados          = n_extreme,
        n_explicado_hematol   = if ("flag_hematol" %in% names(extreme_log)) sum(extreme_log$flag_hematol, na.rm = TRUE) else NA_integer_,
        n_sin_explicacion     = if ("flag_hematol" %in% names(extreme_log)) sum(!extreme_log$flag_hematol | is.na(extreme_log$flag_hematol)) else NA_integer_,
        stringsAsFactors = FALSE
      ),
      file.path(OUTPUT_DIR, paste0("step3_zoutliers_summary_", ba_col, ".csv")),
      row.names = FALSE
    )
    message(">> Guardado resumen agregado (final_output/, sin datos de paciente): step3_zoutliers_summary_", ba_col, ".csv")

    rc[[z_col]][z_extreme_idx] <- NA_real_
    cat(sprintf("    → Z-scores extremos reemplazados por NA. n válidos restantes: %d\n",
                sum(!is.na(rc[[z_col]]))))
  }

  # ── Empirical percentile (ECDF) — avoids Gaussian assumption ──
  pct_col <- paste0(ba_col, "_pct_empirical")
  rc[[pct_col]] <- NA_real_
  if (sum(valid_resid) > 1) {
    ecdf_fn <- ecdf(rc[[resid_col]][valid_resid])
    rc[[pct_col]][valid_resid] <- ecdf_fn(rc[[resid_col]][valid_resid])
  }

  # ── Sex-stratified Z-score — separate SD per sex ──
  # Accounts for males having wider spread (variance ratio M/F ≈ 1.14)
  zsex_col <- paste0(ba_col, "_sexstrat_z")
  rc[[zsex_col]] <- NA_real_
  for (sx in c("F", "M")) {
    idx <- valid_resid & !is.na(rc$sex) & rc$sex == sx
    if (sum(idx) > 1)
      rc[[zsex_col]][idx] <- as.numeric(scale(rc[[resid_col]][idx]))
  }

  cat(sprintf("  %-40s slope_age = %.4f  slope_sexM = %+.3f  n = %d\n",
              ba_col, coefs[2], ifelse(length(coefs) >= 3, coefs[3], NA), sum(ok)))

  # Export coefficients for single-patient scoring (Argentine reference)
  arg_method_b_coefs[[ba_col]] <- list(
    clock      = ba_col,
    intercept  = coefs[1],
    slope_age  = coefs[2],
    slope_sexM = ifelse(length(coefs) >= 3, coefs[3], 0),
    resid_sd   = sd(resid_vals),
    resid_mean = mean(resid_vals),   # should be ~0
    n          = sum(ok),
    sex_levels = c("F", "M"),        # F = reference, M = + slope_sexM
    method     = "lm(BA ~ age + sex) on Argentine cohort"
  )
}

# Save to RDS for use by score_single_patient_v2.R / v3.R
arg_coefs_path <- file.path(DATA_DERIVED_DIR, "argentine_method_b_coefs.rds")
saveRDS(arg_method_b_coefs, arg_coefs_path)
message(">> Saved Argentine Method B coefficients: ", arg_coefs_path)
message("   Clocks: ", paste(names(arg_method_b_coefs), collapse = ", "))

# Save residual distributions for ECDF (used by score_single_patient_v3.R)
# These raw residuals allow exact empirical percentile lookup per patient.
ecdf_ref <- list()
for (ac in viable_adv) {
  ba_col    <- sub("_advance$", "", ac)
  resid_col <- paste0(ba_col, "_sexadj_resid")
  if (resid_col %in% names(rc)) {
    valid_r        <- !is.na(rc[[resid_col]]) & is.finite(rc[[resid_col]])
    ecdf_ref[[ba_col]] <- rc[[resid_col]][valid_r]
  }
}
ecdf_ref_path <- file.path(DATA_DERIVED_DIR, "argentine_ecdf_reference.rds")
saveRDS(ecdf_ref, ecdf_ref_path)
message(">> Saved Argentine ECDF reference distributions: ", ecdf_ref_path)
message("   n por reloj: ", paste(
  paste0(names(ecdf_ref), "=", format(sapply(ecdf_ref, length), big.mark = ",")),
  collapse = ", "))


# ══════════════════════════════════════════════════════════════════════════
# STEP 3-VALIDATION: Statistical Validation of Recentered Z-Scores
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  STEP 3-VALIDATION: Z-Score Diagnostics\n")
cat("══════════════════════════════════════════════════════════════\n\n")

sexadj_z_cols_val <- paste0(sub("_advance$", "", viable_adv), "_sexadj_z")
sexadj_z_cols_val <- intersect(sexadj_z_cols_val, names(rc))

validation_results     <- list()
validation_results_sex <- list()   # separate list — different schema

for (zc in sexadj_z_cols_val) {

  clock_label <- sub("_sexadj_z$", "", zc)
  z_vals <- rc[[zc]]

  ok_z <- !is.na(z_vals) & is.finite(z_vals)
  z_clean <- z_vals[ok_z]
  n_z <- length(z_clean)

  if (n_z < 50) next

  cat("\n------------------------------------------------------------\n")
  cat("  VALIDATION:", clock_label, "| n =", format(n_z, big.mark = ","), "\n")
  cat("------------------------------------------------------------\n")

  # ── BASIC STATS ──
  z_min <- min(z_clean)
  z_max <- max(z_clean)
  z_iqr <- IQR(z_clean)
  z_p01 <- quantile(z_clean, 0.01)
  z_p99 <- quantile(z_clean, 0.99)

  # ── OUTLIERS ──
  n_gt3 <- sum(abs(z_clean) > 3)
  n_gt4 <- sum(abs(z_clean) > 4)
  n_gt5 <- sum(abs(z_clean) > 5)

  cat(sprintf("\n  |Z| > 3: %d (%.2f%%)\n", n_gt3, 100 * n_gt3 / n_z))
  cat(sprintf("  |Z| > 4: %d (%.2f%%)\n", n_gt4, 100 * n_gt4 / n_z))
  cat(sprintf("  |Z| > 5: %d (%.2f%%)\n", n_gt5, 100 * n_gt5 / n_z))

  cat(sprintf("\n  Range: [%.3f, %.3f]\n", z_min, z_max))
  cat(sprintf("  IQR: %.3f | P01: %.3f | P99: %.3f\n", z_iqr, z_p01, z_p99))

  # ── EMPIRICAL CUTS (pooled) ──
  qts <- quantile(z_clean, probs = c(0.10, 0.25, 0.75, 0.90))

  cat("\n  Empirical vs Gaussian cut points (pooled):\n")
  cat(sprintf("    %-5s  %8s  %8s  %+9s\n", "Cut", "Empirical", "Gaussian", "Delta"))
  probs_labels <- c("P10", "P25", "P75", "P90")
  gauss_vals   <- qnorm(c(0.10, 0.25, 0.75, 0.90))
  for (i in seq_along(qts)) {
    delta <- as.numeric(qts[i]) - gauss_vals[i]
    flag  <- if (abs(delta) > 0.1) " FLAG" else ""
    cat(sprintf("    %-5s  %8.4f  %8.4f  %+9.4f%s\n",
                probs_labels[i], qts[i], gauss_vals[i], delta, flag))
  }

  # ── EMPIRICAL CUTS BY SEX ──
  sex_cut_rows <- list()
  for (sx in c("F", "M")) {
    sx_idx <- ok_z & !is.na(rc$sex) & rc$sex == sx
    z_sx   <- z_vals[sx_idx]
    n_sx   <- length(z_sx)
    if (n_sx < 30) next
    qts_sx <- quantile(z_sx, probs = c(0.10, 0.25, 0.75, 0.90))
    cat(sprintf("\n  Empirical cuts — %s (n = %s):\n", sx, format(n_sx, big.mark = ",")))
    for (i in seq_along(qts_sx)) {
      delta_sx <- as.numeric(qts_sx[i]) - gauss_vals[i]
      flag_sx  <- if (abs(delta_sx) > 0.1) " FLAG" else ""
      cat(sprintf("    %-5s  %8.4f  (Gaussian: %+.4f  delta: %+.4f%s)\n",
                  probs_labels[i], qts_sx[i], gauss_vals[i], delta_sx, flag_sx))
    }
    sex_cut_rows[[sx]] <- data.frame(
      clock = clock_label, sex = sx, n = n_sx,
      sd_z  = round(sd(z_sx), 4),
      p10_empirical = round(as.numeric(qts_sx[1]), 4),
      p25_empirical = round(as.numeric(qts_sx[2]), 4),
      p75_empirical = round(as.numeric(qts_sx[3]), 4),
      p90_empirical = round(as.numeric(qts_sx[4]), 4),
      p10_gaussian  = round(gauss_vals[1], 4),
      p90_gaussian  = round(gauss_vals[4], 4),
      stringsAsFactors = FALSE
    )
  }

  # ── STORE ──
  validation_results[[clock_label]] <- data.frame(
    clock         = clock_label,
    n             = n_z,
    z_min         = z_min,
    z_max         = z_max,
    z_iqr         = z_iqr,
    z_p01         = as.numeric(z_p01),
    z_p99         = as.numeric(z_p99),
    n_gt3         = n_gt3,
    n_gt4         = n_gt4,
    n_gt5         = n_gt5,
    p10_empirical = round(as.numeric(qts[1]), 4),
    p25_empirical = round(as.numeric(qts[2]), 4),
    p75_empirical = round(as.numeric(qts[3]), 4),
    p90_empirical = round(as.numeric(qts[4]), 4),
    p10_gaussian  = round(gauss_vals[1], 4),
    p25_gaussian  = round(gauss_vals[2], 4),
    p75_gaussian  = round(gauss_vals[3], 4),
    p90_gaussian  = round(gauss_vals[4], 4),
    stringsAsFactors = FALSE
  )

  if (length(sex_cut_rows) > 0)
    validation_results_sex[[clock_label]] <- do.call(rbind, sex_cut_rows)
}

# ── SAVE ──
if (length(validation_results) > 0) {
  val_df <- do.call(rbind, validation_results)
  rownames(val_df) <- NULL
  write.csv(val_df,
            file.path(OUTPUT_DIR, "step3_validation_summary.csv"),
            row.names = FALSE)
  message(">> Saved: step3_validation_summary.csv")

  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("  VALIDATION SUMMARY (pooled)\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  print(val_df, row.names = FALSE)
}

if (length(validation_results_sex) > 0) {
  val_sex_df <- do.call(rbind, validation_results_sex)
  rownames(val_sex_df) <- NULL
  write.csv(val_sex_df,
            file.path(OUTPUT_DIR, "step3_validation_summary_by_sex.csv"),
            row.names = FALSE)
  message(">> Saved: step3_validation_summary_by_sex.csv")

  cat("\n══════════════════════════════════════════════════════════════\n")
  cat("  VALIDATION SUMMARY (by sex)\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  print(val_sex_df, row.names = FALSE)
}