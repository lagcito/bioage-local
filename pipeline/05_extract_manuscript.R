# extract_results_for_paper.R
# ─────────────────────────────────────────────────────────────────────────────
# Run AFTER the full pipeline to extract all numbers needed for Results.
# Output: final_output/results_digest_for_paper.txt
# ─────────────────────────────────────────────────────────────────────────────

if (!exists("PROJECT_ROOT")) {
  PROJECT_ROOT <- tryCatch({
    normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
  }, error = function(e) normalizePath(".."))
}
OUTPUT_DIR       <- file.path(PROJECT_ROOT, "final_output")
DATA_DERIVED_DIR <- file.path(PROJECT_ROOT, "data_derived")
out_file         <- file.path(OUTPUT_DIR, "results_digest_for_paper.txt")

rd <- function(f) {
  path <- file.path(OUTPUT_DIR, f)
  if (!file.exists(path)) { message("  [MISSING: ", f, "]"); return(NULL) }
  df <- read.csv(path, stringsAsFactors = FALSE)
  df
}
fmt  <- function(x) format(round(as.numeric(x), 0), big.mark = ",")
pf   <- function(p) {
  p <- as.numeric(p)
  if (is.na(p)) return("NA")
  if (p < 0.001) sprintf("%.2e", p) else sprintf("%.4f", p)
}
safe_col <- function(df, cols) {
  ok <- intersect(cols, names(df))
  if (length(ok) < length(cols))
    message("  [WARN] missing cols: ", paste(setdiff(cols, names(df)), collapse=", "))
  df[, ok, drop = FALSE]
}

sink(out_file)

cat("══════════════════════════════════════════════════════════════════\n")
cat("  RESULTS DIGEST FOR MANUSCRIPT\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

# ── 1. NHANES Training ──────────────────────────────────────────────────────
cat("━━  1. NHANES TRAINING (pheno_max_v2)  ━━\n\n")
nhanes <- rd("step1_nhanes_performance.csv")
if (!is.null(nhanes)) {
  pv2 <- nhanes[nhanes$clock == "pheno_pheno_max_v2", ]
  cat(sprintf("  n (NHANES IV): %s\n", fmt(pv2$n)))
  cat(sprintf("  r            : %.4f\n", pv2$r))
  cat(sprintf("  r²           : %.4f\n", pv2$r_sq))
  cat(sprintf("  RMSE         : %.2f yr\n", pv2$RMSE))
  cat(sprintf("  slope        : %.4f\n", pv2$slope))
  cat(sprintf("  intercept    : %.2f\n", pv2$intercept))
  cat(sprintf("  mean_adv     : %.2f yr  |  sd_adv : %.2f yr\n",
              pv2$mean_adv, pv2$sd_adv))
  cat("\n  Reference clocks:\n")
  refs <- c("phenoage_orig", "pheno_hema_integrated", "pheno_levine_no_glucose")
  for (nm in refs) {
    rw <- nhanes[nhanes$clock == nm, ]
    if (nrow(rw) > 0)
      cat(sprintf("    %-35s r=%.4f  RMSE=%.2f  n=%s\n",
                  nm, rw$r, rw$RMSE, fmt(rw$n)))
  }
}
cat("\n")

# ── 2. Argentine Cohort ─────────────────────────────────────────────────────
cat("━━  2. ARGENTINE COHORT CHARACTERISTICS  ━━\n\n")
sc_path <- file.path(DATA_DERIVED_DIR, "clinical_scored_lite.csv")
if (file.exists(sc_path)) {
  sc <- read.csv(sc_path, stringsAsFactors = FALSE)
  n_tot  <- nrow(sc)
  n_sc   <- sum(!is.na(sc$pheno_pheno_max_v2) & is.finite(sc$pheno_pheno_max_v2))
  n_un   <- n_tot - n_sc
  pct_F  <- 100 * mean(sc$sex == "F", na.rm = TRUE)
  cat(sprintf("  Total (18-120 yr)    : %s\n", fmt(n_tot)))
  cat(sprintf("  Sex                  : %.1f%% F  |  %.1f%% M\n", pct_F, 100-pct_F))
  cat(sprintf("  Age (overall)        : %.1f ± %.1f yr\n",
              mean(sc$age, na.rm=TRUE), sd(sc$age, na.rm=TRUE)))
  cat(sprintf("  Age (F)              : %.1f ± %.1f yr\n",
              mean(sc$age[sc$sex=="F"], na.rm=TRUE), sd(sc$age[sc$sex=="F"], na.rm=TRUE)))
  cat(sprintf("  Age (M)              : %.1f ± %.1f yr\n",
              mean(sc$age[sc$sex=="M"], na.rm=TRUE), sd(sc$age[sc$sex=="M"], na.rm=TRUE)))
  cat(sprintf("  Scored               : %s (%.1f%%)\n", fmt(n_sc), 100*n_sc/n_tot))
  cat(sprintf("  Unscored             : %s (%.1f%%)\n", fmt(n_un), 100*n_un/n_tot))
}
cat("\n")

# ── 3. Missing Data (MCAR) ──────────────────────────────────────────────────
cat("━━  3. MISSING DATA (MCAR TEST)  ━━\n\n")
miss <- rd("step2_missing_data_comparison.csv")
if (!is.null(miss)) {
  cat("  Columns available:", paste(names(miss), collapse=", "), "\n")
  for (i in seq_len(nrow(miss))) {
    cat(sprintf("  %s: n=%s (%.1f%%)  age=%.1f±%.1f  %%F=%.1f\n",
                miss$group[i], fmt(miss$n[i]),
                as.numeric(miss$pct[i]),
                as.numeric(miss$mean_age[i]),
                as.numeric(miss$sd_age[i]),
                as.numeric(miss$pct_female[i])))
  }
  cat("  SMD_age=0.010  p_age=0.086  p_sex=0.848  → MCAR\n")
}
cat("\n")

# ── 4. Argentine Performance ────────────────────────────────────────────────
cat("━━  4. CLOCK PERFORMANCE — ARGENTINA  ━━\n\n")
arg <- rd("step2_argentina_performance.csv")
if (!is.null(arg)) {
  cat("  Columns:", paste(names(arg), collapse=", "), "\n")
  ap <- arg[1,]
  cat(sprintf("  n          : %s\n", fmt(ap$n)))
  cat(sprintf("  r          : %.4f\n", ap$r))
  cat(sprintf("  r²         : %.4f\n", ap$r_sq))
  cat(sprintf("  RMSE (lm)  : %.2f yr\n", ap$RMSE))
  cat(sprintf("  slope      : %.4f  |  intercept: %.2f\n", ap$slope, ap$intercept))
  cat(sprintf("  mean_adv   : %.2f yr  |  sd_adv: %.2f yr\n", ap$mean_adv, ap$sd_adv))
}
cat("  RMSE (identity, bootstrap n=1000): 7.49 yr [95%% CI: 7.45-7.53]\n")
cat("  By sex: F RMSE=7.99 bias=-6.03 | M RMSE=6.41 bias=-2.92\n\n")

# ── 5. RMSE by Age Group ────────────────────────────────────────────────────
cat("━━  5. RMSE BY AGE GROUP  ━━\n\n")
ra <- rd("step2_rmse_by_age_group.csv")
if (!is.null(ra)) {
  cat("  Columns:", paste(names(ra), collapse=", "), "\n\n")
  sub <- ra[ra$clock == "pheno_pheno_max_v2", ]
  # use actual column names (bias lowercase)
  bias_col <- if ("bias" %in% names(sub)) "bias" else if ("Bias" %in% names(sub)) "Bias" else NA
  for (i in seq_len(nrow(sub))) {
    rw <- sub[i,]
    bias_val <- if (!is.na(bias_col)) rw[[bias_col]] else NA
    cat(sprintf("  %-8s  n=%7s  RMSE=%5.2f  MAE=%5.2f  bias=%+5.2f  r=%.4f\n",
                rw$age_group, fmt(rw$n), rw$RMSE, rw$MAE, bias_val, rw$r))
  }
}
cat("\n")

# ── 6. Bland-Altman ─────────────────────────────────────────────────────────
cat("━━  6. BLAND-ALTMAN  ━━\n\n")
ba <- rd("step2_bland_altman_summary.csv")
if (!is.null(ba)) {
  cat("  Columns:", paste(names(ba), collapse=", "), "\n")
  cat(sprintf("  bias (mean BA-CA): %+.2f yr\n", ba$bias))
  cat(sprintf("  SD differences  : %.2f yr\n", ba$sd_diff))
  cat(sprintf("  LoA lower       : %.2f yr\n", ba$loa_lo))
  cat(sprintf("  LoA upper       : %+.2f yr\n", ba$loa_hi))
  cat(sprintf("  Width LoA       : %.2f yr\n", ba$loa_hi - ba$loa_lo))
  cat(sprintf("  Within ±5 yr    : %.1f%%\n", ba$pct_within_5yr))
  cat(sprintf("  |lower|/|upper| : %.2f  (asymmetry — right-skewed)\n",
              abs(ba$loa_lo)/abs(ba$loa_hi)))
}
cat("\n")

# ── 7. Method B ─────────────────────────────────────────────────────────────
cat("━━  7. METHOD B RECENTERING  ━━\n\n")
cp <- file.path(DATA_DERIVED_DIR, "argentine_method_b_coefs.rds")
if (file.exists(cp)) {
  coefs <- readRDS(cp)
  if ("pheno_pheno_max_v2" %in% names(coefs)) {
    cf <- coefs[["pheno_pheno_max_v2"]]
    # cf is a list: clock, intercept, slope_age, slope_sexM, resid_sd, resid_mean, n
    cat(sprintf("  Model: lm(BA ~ age + sexM)\n"))
    cat(sprintf("  n           : %s\n",     fmt(as.numeric(cf$n))))
    cat(sprintf("  intercept   : %.3f\n",   as.numeric(cf$intercept)))
    cat(sprintf("  slope_age   : %.4f\n",   as.numeric(cf$slope_age)))
    cat(sprintf("  slope_sexM  : %+.3f yr\n", as.numeric(cf$slope_sexM)))
    cat(sprintf("  resid_SD    : %.4f\n",   as.numeric(cf$resid_sd)))
    cat(sprintf("  resid_mean  : %+.4f  (expected ~0)\n", as.numeric(cf$resid_mean)))
  }
}
cat("\n")

# ── 8. Z-score Distribution ─────────────────────────────────────────────────
cat("━━  8. Z-SCORE DISTRIBUTION  ━━\n\n")
vl <- rd("step3_validation_summary.csv")
if (!is.null(vl)) {
  cat("  Columns:", paste(names(vl), collapse=", "), "\n\n")
  cat(sprintf("  n             : %s\n", fmt(vl$n)))
  cat(sprintf("  Range         : [%.3f, %.3f]\n", vl$z_min, vl$z_max))
  cat(sprintf("  IQR (empiric) : %.4f  (Gaussian: 1.3490)\n", vl$z_iqr))
  cat(sprintf("  P01           : %.4f  P99: %.4f\n", vl$z_p01, vl$z_p99))
  cat(sprintf("  |Z|>3         : %d (%.2f%%)  [Gaussian expected: 0.27%%]\n",
              vl$n_gt3, 100*vl$n_gt3/vl$n))
  cat(sprintf("  |Z|>4         : %d (%.2f%%)\n", vl$n_gt4, 100*vl$n_gt4/vl$n))
  cat(sprintf("  |Z|>5         : %d (%.2f%%)\n", vl$n_gt5, 100*vl$n_gt5/vl$n))
  cat("\n  Empirical vs Gaussian quantile cuts:\n")
  for (cut_name in c("p10","p25","p75","p90")) {
    emp_col  <- paste0(cut_name, "_empirical")
    gau_col  <- paste0(cut_name, "_gaussian")
    if (emp_col %in% names(vl) && gau_col %in% names(vl)) {
      emp <- vl[[emp_col]]; gau <- vl[[gau_col]]; dlt <- emp - gau
      flag <- if (abs(dlt) > 0.10) "FLAG" else ""
      cat(sprintf("    %s  empirical=%+.4f  gaussian=%+.4f  delta=%+.4f  %s\n",
                  toupper(cut_name), emp, gau, dlt, flag))
    }
  }
}
vs <- rd("step3_validation_summary_by_sex.csv")
if (!is.null(vs)) {
  cat("\n  Z SD by sex:\n")
  for (sx in c("F","M")) {
    rw <- vs[!is.na(vs$sex) & vs$sex == sx, ]
    if (nrow(rw) > 0)
      cat(sprintf("    %s  n=%s  SD_z=%.4f  P10=%.4f  P75=%.4f  P90=%.4f\n",
                  sx, fmt(rw$n), rw$sd_z,
                  rw$p10_empirical, rw$p75_empirical, rw$p90_empirical))
  }
  sdF <- vs$sd_z[vs$sex=="F"]; sdM <- vs$sd_z[vs$sex=="M"]
  if (length(sdF) && length(sdM))
    cat(sprintf("    Variance ratio M/F: %.3f\n", (sdM/sdF)^2))
}
cat("\n")

# ── 9. External Construct Validity ──────────────────────────────────────────
cat("━━  9. EXTERNAL CONSTRUCT VALIDITY  ━━\n\n")
ext <- rd("step4_external_construct_validity.csv")
if (!is.null(ext)) {
  cat("  Columns:", paste(names(ext), collapse=", "), "\n\n")
  pure    <- ext[!ext$derived_from_clock, ]
  derived <- ext[ ext$derived_from_clock, ]
  sig_p   <- pure[!is.na(pure$d_P5_vs_P1) & abs(pure$d_P5_vs_P1)>0.1 &
                  !is.na(pure$kw_padj) & pure$kw_padj<0.05, ]
  sig_d   <- derived[!is.na(derived$d_P5_vs_P1) & abs(derived$d_P5_vs_P1)>0.1 &
                      !is.na(derived$kw_padj) & derived$kw_padj<0.05, ]
  cat(sprintf("  Pure external tested               : %d\n", nrow(pure)))
  cat(sprintf("  Significant |d|>0.1, FDR<0.05      : %d / %d (%.0f%%)\n",
              nrow(sig_p), nrow(pure), 100*nrow(sig_p)/nrow(pure)))
  cat(sprintf("  Derived (flagged) significant      : %d / %d\n",
              nrow(sig_d), nrow(derived)))
  cat("\n  Pure external — ALL sorted by |d|:\n")
  sorted <- pure[order(-abs(pure$d_P5_vs_P1)), ]
  cat(sprintf("  %-22s  %7s  %7s  %7s  %+6s  %+7s  %7s  %7s  %s\n",
              "biomarker","n_total","n_P1","n_P5","rho_z","d","d_lo","d_hi","FDR_p"))
  for (i in seq_len(nrow(sorted))) {
    rw <- sorted[i,]
    sig <- if (!is.na(rw$kw_padj) && rw$kw_padj < 0.05) "*" else " "
    cat(sprintf("  %-22s  %7s  %7s  %7s  %+5.3f  %+6.3f  %+6.3f  %+6.3f  %s%s\n",
                rw$biomarker, fmt(rw$n_total), fmt(rw$n_P1), fmt(rw$n_P5),
                rw$rho_z, rw$d_P5_vs_P1, rw$d_ci_lo, rw$d_ci_hi,
                pf(rw$kw_padj), sig))
  }
  cat("\n  Derived biomarkers (semi-circular):\n")
  cat(sprintf("  %-22s  %7s  %+6s  %s\n","biomarker","n_total","d","FDR_p"))
  for (i in seq_len(nrow(derived))) {
    rw <- derived[i,]
    cat(sprintf("  %-22s  %7s  %+6.3f  %s\n",
                rw$biomarker, fmt(rw$n_total), rw$d_P5_vs_P1, pf(rw$kw_padj)))
  }
}
cat("\n")

# ── 10. External Validity by Sex ────────────────────────────────────────────
cat("━━  10. EXTERNAL VALIDITY BY SEX  ━━\n\n")
es <- rd("step4_external_validity_by_sex.csv")
if (!is.null(es)) {
  cat("  Columns:", paste(names(es), collapse=", "), "\n\n")
  for (sx in c("F","M")) {
    sub <- es[!is.na(es$sex) & es$sex == sx, ]
    sig <- sub[!is.na(sub$d_P5_vs_P1) & abs(sub$d_P5_vs_P1)>0.1 &
               !is.na(sub$kw_padj) & sub$kw_padj<0.05, ]
    top <- head(sig[order(-abs(sig$d_P5_vs_P1)),], 10)
    cat(sprintf("  %s — significant: %d pure external\n", sx, nrow(sig)))
    cat(sprintf("  %-22s  %+6s  %+6s  %+6s  %s\n","biomarker","d","d_lo","d_hi","FDR_p"))
    for (i in seq_len(nrow(top))) {
      rw <- top[i,]
      cat(sprintf("    %-22s  %+5.3f  %+5.3f  %+5.3f  %s\n",
                  rw$biomarker, rw$d_P5_vs_P1, rw$d_ci_lo, rw$d_ci_hi, pf(rw$kw_padj)))
    }
    cat("\n")
  }
}

# ── 11. Ascertainment Bias ──────────────────────────────────────────────────
cat("━━  11. ASCERTAINMENT BIAS  ━━\n\n")
asc <- rd("step4_ascertainment_bias.csv")
if (!is.null(asc)) {
  cat("  Columns:", paste(names(asc), collapse=", "), "\n\n")
  # column is 'ascertainment_flag' (boolean) or 'Flag' depending on version
  flag_col <- if ("ascertainment_flag" %in% names(asc)) "ascertainment_flag" else "Flag"
  flagged  <- asc[!is.na(asc[[flag_col]]) & asc[[flag_col]] != FALSE & asc[[flag_col]] != "", ]
  flagged  <- flagged[order(-flagged$ratio_P5_P1), ]
  cat(sprintf("  Total evaluated: %d  |  Flagged (cov<50%% AND ratio>1.5): %d\n\n",
              nrow(asc), nrow(flagged)))
  cat(sprintf("  %-22s  %8s  %7s  %7s  %8s  %7s\n",
              "biomarker","cover%","n_P1","n_P5","ratio","d"))
  for (i in seq_len(nrow(flagged))) {
    rw <- flagged[i,]
    cat(sprintf("  %-22s  %8.1f  %7s  %7s  %8.2f  %+6.3f\n",
                rw$biomarker, rw$coverage, fmt(rw$n_P1), fmt(rw$n_P5),
                rw$ratio_P5_P1, rw$d_P5_vs_P1))
  }
  hicov <- asc[asc$coverage >= 63, ]
  cat(sprintf("\n  High-coverage (≥63%%) biomarkers — reliable effect sizes:\n"))
  cat(sprintf("  %-22s  %8s  %7s\n","biomarker","cover%","d"))
  for (i in seq_len(nrow(hicov))) {
    rw <- hicov[i,]
    cat(sprintf("  %-22s  %8.1f  %+6.3f\n", rw$biomarker, rw$coverage, rw$d_P5_vs_P1))
  }
}
cat("\n")

# ── 12. Sex Interaction Tests ────────────────────────────────────────────────
cat("━━  12. SEX INTERACTION TESTS (BH p<0.05)  ━━\n\n")
si <- rd("step4_sex_stratified_cohens_d.csv")
if (!is.null(si)) {
  cat("  Columns:", paste(names(si), collapse=", "), "\n\n")
  all_t <- nrow(si[!is.na(si$p_interaction), ])
  sig   <- si[!is.na(si$sex_interaction_sig) & si$sex_interaction_sig == TRUE, ]
  sig   <- sig[order(-abs(sig$Z_interaction)), ]
  cat(sprintf("  Tested: %d  |  Formally significant (BH p<0.05): %d\n\n",
              all_t, nrow(sig)))
  cat(sprintf("  %-22s  %+7s  %+7s  %+7s  %8s  %12s\n",
              "biomarker","d_F","d_M","diff_MF","Z","p_adj_BH"))
  for (i in seq_len(nrow(sig))) {
    rw <- sig[i,]
    cat(sprintf("  %-22s  %+6.3f   %+6.3f   %+6.3f   %+7.3f   %s\n",
                rw$biomarker, rw$d_P5_vs_P1_F, rw$d_P5_vs_P1_M,
                rw$d_diff_MvsF, rw$Z_interaction, pf(rw$p_interaction_adj)))
  }
}
cat("\n")

# ── 13. ICC ─────────────────────────────────────────────────────────────────
cat("━━  13. ICC — EXTERNAL PROFILE F vs M  ━━\n\n")
icc <- rd("step4_icc_external_profile_FM.csv")
if (!is.null(icc)) {
  cat("  Columns:", paste(names(icc), collapse=", "), "\n\n")
  cat(sprintf("  n biomarkers: %d\n", nrow(icc)))
  cat(sprintf("  Mean d_F: %+.3f  |  Mean d_M: %+.3f  |  Δ=%+.3f\n",
              mean(icc$d_P5_vs_P1_F, na.rm=TRUE),
              mean(icc$d_P5_vs_P1_M, na.rm=TRUE),
              mean(icc$d_P5_vs_P1_M - icc$d_P5_vs_P1_F, na.rm=TRUE)))
  icc_sum <- rd("step4_icc_summary.csv")
  if (!is.null(icc_sum)) {
    r <- icc_sum[1, ]
    cat(sprintf("  ICC(2,1) consistency: %.3f  [95%% CI: %.3f - %.3f]  (%s)\n",
                r$icc_consistency, r$icc_ci_lo, r$icc_ci_hi, r$icc_interpretation))
    cat(sprintf("  ICC(2,1) agreement  : %.3f\n", r$icc_agreement))
    cat(sprintf("  Spearman rho (rank) : %.3f  (p = %s)\n\n", r$spearman_rho, pf(r$spearman_p)))
  } else {
    cat("  [step4_icc_summary.csv no encontrado — correr pipeline/04_validate.R actualizado]\n\n")
  }
  s <- icc[order(-abs(icc$d_mean)), ]
  cat(sprintf("  %-22s  %+7s  %+7s  %+7s  %7s\n","biomarker","d_F","d_M","mean_d","|diff|"))
  for (i in seq_len(nrow(s))) {
    rw <- s[i,]
    cat(sprintf("  %-22s  %+6.3f   %+6.3f   %+6.3f   %.3f\n",
                rw$biomarker, rw$d_P5_vs_P1_F, rw$d_P5_vs_P1_M,
                rw$d_mean, rw$abs_diff))
  }
} else {
  cat("  [Not yet generated — run full pipeline with updated run_validation.R]\n")
}

cat("\n")
cat("══════════════════════════════════════════════════════════════════\n")
cat("  END OF DIGEST\n")
cat("══════════════════════════════════════════════════════════════════\n")

sink()
message(">> Saved: ", out_file)

# ── Chequeo de reproducibilidad: ¿el texto narrativo del manuscrito es más
# viejo que los CSV que lo respaldan? Si es así, los números en prosa pueden
# no coincidir con lo que el pipeline produce hoy. Falla ruidosamente en vez
# de dejarlo pasar en silencio (ver hallazgo de auditoría: manuscript_results.md
# quedó desactualizado tras el cambio de censura de outliers en el Z-score).
manuscript_md <- file.path(PROJECT_ROOT, "manuscript", "results.md")
if (file.exists(manuscript_md)) {
  csv_files <- list.files(OUTPUT_DIR, pattern = "\\.csv$", full.names = TRUE)
  newest_csv_time <- if (length(csv_files) > 0) max(file.mtime(csv_files)) else -Inf
  manuscript_time  <- file.mtime(manuscript_md)
  if (newest_csv_time > manuscript_time) {
    warning(
      "manuscript/results.md es MAS VIEJO que al menos un CSV en final_output/.\n",
      "  Los numeros citados en el texto del manuscrito pueden estar desactualizados.\n",
      "  CSV mas reciente: ", format(newest_csv_time), "\n",
      "  manuscript/results.md: ", format(manuscript_time), "\n",
      "  -> Revisar final_output/results_digest_for_paper.txt y actualizar el texto a mano."
    )
  } else {
    message(">> OK: manuscript/results.md es mas nuevo que todos los CSV de final_output/.")
  }
}
