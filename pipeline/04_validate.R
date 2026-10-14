# run_validation.R
#  ══════════════════════════════════════════════════════════════════════════
# STEP 4. FAST/SLOW AGER PROFILING & BIOMARKER SIGNATURES
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  STEP 4: FAST vs SLOW AGER PROFILING                   ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

# Load profiling module
source(file.path(PROJECT_ROOT, "R/profile.R"), local = FALSE)

# Use sex-adjusted Z-scores for profiling (best Z, adjusted for age + sex)
z_cols <- grep("_sexadj_z$", names(rc), value = TRUE)
if (length(z_cols) == 0) z_cols <- grep("_local_z$", names(rc), value = TRUE)

# Pick highest-coverage Z-score as primary
z_coverage <- sapply(z_cols, function(zc) sum(!is.na(rc[[zc]])))
primary_z <- z_cols[which.max(z_coverage)]
message(">> Primary Z-score for profiling: ", primary_z,
        " (n = ", max(z_coverage), ")")

# Load raw clinical data for ALL biomarkers
raw <- cohort
raw$patient_id <- raw$Protocolo
raw$age <- as.numeric(raw$age)
raw <- raw[!is.na(raw$age) & raw$age >= 18 & raw$age <= 120, ]

# Merge recentered scores with raw biomarkers
prof_df <- merge(rc, raw[, c("patient_id", setdiff(names(raw),
            c(names(rc), "Protocolo", "source_file",
              "date_admission", "date_result")))],
            by = "patient_id", all.x = TRUE)

# Convert biomarkers to numeric
all_biomarkers <- c(
  "hemoglobin", "hematocrit", "platelets", "mch", "mchc",
  "neutrophils", "eosinophils", "basophils", "monocytes",
  "lymphocyte_abs", "band_neutrophils",
  "cholesterol", "hdl", "ldl", "non_hdl", "triglycerides", "lp_a",
  "glucose", "hba1c", "insulin",
  "ast", "alt", "ggt", "ldh", "albumin", "alp",
  "ferritin", "transferrin", "transferrin_saturation",
  "creatinine", "bun", "uric_acid", "eGFR",
  "crp", "fibrinogen", "homocysteine",
  "TSH", "free_t4",
  "vitamin_d", "vitamin_b12",
  "wbc", "rbc", "mcv", "rdw", "lymphocyte",
  "nlr", "homa_ir"
)
all_biomarkers <- intersect(all_biomarkers, names(prof_df))
for (bm in all_biomarkers) prof_df[[bm]] <- suppressWarnings(as.numeric(prof_df[[bm]]))

# ── 4a. Classify fast vs slow agers (quintile-based) ──

ok      <- !is.na(prof_df[[primary_z]]) & is.finite(prof_df[[primary_z]])
n_valid <- sum(ok)
message(">> Valid subjects for profiling: ", n_valid)

# Derive ager_group here for distribution display (compute_biomarker_profile
# recomputes it internally — single source of truth via the function)
quants <- quantile(prof_df[[primary_z]][ok], probs = c(0.10, 0.25, 0.75, 0.90))
prof_df$ager_group <- NA_character_
prof_df$ager_group[ok] <- as.character(cut(
  prof_df[[primary_z]][ok],
  breaks = c(-Inf, quants[1], quants[2], quants[3], quants[4], Inf),
  labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
))
df_valid <- prof_df[ok, ]

cat("\n--- Ager Group Distribution ---\n\n")
tbl <- table(df_valid$ager_group)
for (g in names(tbl))
  cat(sprintf("   %-15s n = %6d (%5.1f%%)\n", g, tbl[g], 100 * tbl[g] / n_valid))

# ── 4b. Biomarker profiling across ager groups ──

cat("\n--- Biomarker Profiles: Fast vs Slow Agers ---\n\n")

usable_bm <- all_biomarkers[sapply(all_biomarkers, function(bm) {
  bm %in% names(prof_df) &&
  sum(!is.na(prof_df[[bm]]) & is.finite(prof_df[[bm]])) >= 50
})]
message(">> Usable biomarkers: ", length(usable_bm))

profile_table <- compute_biomarker_profile(prof_df, primary_z, usable_bm)

if (!is.null(profile_table)) {
  sig <- profile_table[!is.na(profile_table$d_P5_vs_P1) &
                       abs(profile_table$d_P5_vs_P1) > 0.1, ]
  if (nrow(sig) > 0) {
    cat("  Significantly differentiating biomarkers (|d| > 0.1):\n\n")
    print(sig[, c("biomarker", "n_total", "n_P1", "n_P5", "rho_z", "d_P5_vs_P1",
                   "d_ci_lo", "d_ci_hi", "mean_P1", "mean_P5", "direction", "kw_padj")],
          row.names = FALSE)
  }
  write.csv(profile_table, file.path(OUTPUT_DIR, "step4_ager_biomarker_profiles.csv"),
            row.names = FALSE)
  message(">> Saved: step4_ager_biomarker_profiles.csv")

  # ── Ascertainment bias table ──
  # Biomarkers with low coverage may show inflated Cohen's d because
  # physicians order them selectively for symptomatic (fast-aging) patients.
  # Flag: coverage < 50% AND n_P5/n_P1 ratio > 1.5
  cat("\n--- Ascertainment Bias Check ---\n\n")
  cat(sprintf("  %-28s %6s %6s %6s %7s %8s %s\n",
              "Biomarker", "Cover%", "n_P1", "n_P5", "Ratio", "d", "Flag"))
  asc_rows <- list()
  for (bm in profile_table$biomarker) {
    row <- profile_table[profile_table$biomarker == bm, ]
    n_total_bm <- sum(!is.na(prof_df[[bm]]) & is.finite(prof_df[[bm]]))
    coverage   <- round(100 * n_total_bm / nrow(prof_df), 1)
    ratio      <- if (!is.na(row$n_P1) && row$n_P1 > 0)
                    round(row$n_P5 / row$n_P1, 2) else NA_real_
    flag <- if (!is.na(coverage) && !is.na(ratio) &&
                coverage < 50 && ratio > 1.5) "FLAG" else ""
    cat(sprintf("  %-28s %6.1f %6d %6d %7.2f %8.3f %s\n",
                bm, coverage, row$n_P1, row$n_P5,
                ifelse(is.na(ratio), NA, ratio),
                ifelse(is.na(row$d_P5_vs_P1), NA, row$d_P5_vs_P1),
                flag))
    asc_rows[[bm]] <- data.frame(
      biomarker  = bm,
      coverage   = coverage,
      n_P1       = row$n_P1,
      n_P5       = row$n_P5,
      ratio_P5_P1 = ratio,
      d_P5_vs_P1  = row$d_P5_vs_P1,
      ascertainment_flag = nchar(flag) > 0,
      stringsAsFactors = FALSE
    )
  }
  asc_table <- do.call(rbind, asc_rows)
  rownames(asc_table) <- NULL
  write.csv(asc_table, file.path(OUTPUT_DIR, "step4_ascertainment_bias.csv"),
            row.names = FALSE)
  n_flagged_asc <- sum(asc_table$ascertainment_flag, na.rm = TRUE)
  message(sprintf(">> Saved: step4_ascertainment_bias.csv  (%d biomarkers flagged)",
                  n_flagged_asc))
}

# ── 4c. Aging signature heatmap ──

if (!is.null(profile_table) && nrow(profile_table) > 0) {
  n_sig_bm <- sum(!is.na(profile_table$d_P5_vs_P1) &
                  abs(profile_table$d_P5_vs_P1) > 0.15 &
                  !is.na(profile_table$kw_padj) &
                  profile_table$kw_padj < 0.05)
  cat(sprintf("\n>> Building aging signature from %d significant biomarkers\n\n", n_sig_bm))

  heat <- build_aging_heatmap(profile_table, prof_df, primary_z)

  if (!is.null(heat)) {
    p_heat <- heat$plot +
      ggplot2::labs(title = "Accelerated Aging Biomarker Signature",
                    subtitle = paste0("Based on ", primary_z, " quintiles"))

    ggsave(file.path(OUTPUT_DIR, "step4_aging_signature_heatmap.png"),
           p_heat, width = 10, height = max(6, n_sig_bm * 0.35 + 2), dpi = 150)
    print(p_heat)
    message(">> Saved: step4_aging_signature_heatmap.png")

    write.csv(heat$wide, file.path(OUTPUT_DIR, "step4_aging_signature_data.csv"),
              row.names = FALSE)
    cat("\n--- Accelerated Aging Signature ---\n\n")
    print(heat$wide, row.names = FALSE)

    # Enriched CSV: merge per-group means with Cohen's d and 95% CIs
    heatmap_enriched <- merge(
      heat$data,
      profile_table[, c("biomarker", "n_total", "n_P1", "n_P5",
                        "d_P5_vs_P1", "d_ci_lo", "d_ci_hi", "rho_z", "kw_padj")],
      by = "biomarker", all.x = TRUE
    )
    heatmap_enriched <- heatmap_enriched[order(heatmap_enriched$biomarker,
                                               heatmap_enriched$group), ]
    write.csv(heatmap_enriched, file.path(OUTPUT_DIR, "step4_heatmap_with_cohens_d.csv"),
              row.names = FALSE)
    message(">> Saved: step4_heatmap_with_cohens_d.csv")

    cat("\n--- Cohen's d 95% CIs for Heatmap Biomarkers ---\n\n")
    d_summary <- unique(heatmap_enriched[, c("biomarker", "n_P1", "n_P5",
                                              "d_P5_vs_P1", "d_ci_lo", "d_ci_hi")])
    d_summary <- d_summary[order(-abs(d_summary$d_P5_vs_P1)), ]
    print(d_summary, row.names = FALSE)
  }
}

# ── 4d. Per-clock fast/slow ager profiling ──

cat("\n--- Per-Clock Ager Classification ---\n\n")
cat("  Comparing ager groups defined by EACH sub-clock independently\n\n")

for (zc in z_cols) {
  clock_label <- sub("_sexadj_z$|_local_z$", "", zc)
  z_vals <- prof_df[[zc]]
  ok_z <- !is.na(z_vals) & is.finite(z_vals)
  if (sum(ok_z) < 100) next

  # Classify using this clock's Z-score
  qts <- quantile(z_vals[ok_z], probs = c(0.25, 0.75))
  slow <- ok_z & z_vals <= qts[1]
  fast <- ok_z & z_vals >= qts[2]

  cat(sprintf("  %s: slow (Q1) = %d, fast (Q4) = %d\n",
              clock_label, sum(slow), sum(fast)))

  # For each usable biomarker, compare fast vs slow
  clock_diffs <- list()
  for (bm in usable_bm) {
    x <- prof_df[[bm]]
    x_slow <- x[slow & !is.na(x) & is.finite(x)]
    x_fast <- x[fast & !is.na(x) & is.finite(x)]
    if (length(x_slow) < 20 || length(x_fast) < 20) next

    ps <- sqrt((var(x_slow) * (length(x_slow) - 1) + var(x_fast) * (length(x_fast) - 1)) /
               (length(x_slow) + length(x_fast) - 2))
    d <- if (ps > 0) (mean(x_fast) - mean(x_slow)) / ps else NA

    clock_diffs[[bm]] <- data.frame(
      clock = clock_label, biomarker = bm,
      d = round(d, 3),
      mean_slow = round(mean(x_slow), 3),
      mean_fast = round(mean(x_fast), 3),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  if (length(clock_diffs) > 0) {
    cd <- do.call(rbind, clock_diffs)
    cd <- cd[order(-abs(cd$d)), ]
    # Show top 5
    top <- head(cd, 5)
    for (i in seq_len(nrow(top))) {
      cat(sprintf("    %-25s d = %+.3f  (slow: %.1f → fast: %.1f)\n",
                  top$biomarker[i], top$d[i], top$mean_slow[i], top$mean_fast[i]))
    }
    cat("\n")
  }
}

# ── 4d-bis. FULL biomarker profile using pheno_max_v2 as classifier ──

cat("\n##########################################################\n")
cat("#  FULL BIOMARKER PROFILE — pheno_max_v2 classifier      #\n")
cat("##########################################################\n\n")

maxv2_z <- "pheno_pheno_max_v2_sexadj_z"
if (!maxv2_z %in% names(prof_df)) maxv2_z <- "pheno_pheno_max_v2_local_z"

if (maxv2_z %in% names(prof_df)) {
  z_v2 <- prof_df[[maxv2_z]]
  ok_v2 <- !is.na(z_v2) & is.finite(z_v2)
  n_v2 <- sum(ok_v2)
  message(">> pheno_max_v2 classifier: ", n_v2, " subjects")

  # Quintile classification
  qts_v2 <- quantile(z_v2[ok_v2], probs = c(0.10, 0.25, 0.75, 0.90), na.rm = TRUE)
  prof_df$ager_maxv2 <- NA_character_
  prof_df$ager_maxv2[ok_v2] <- as.character(cut(z_v2[ok_v2],
    breaks = c(-Inf, qts_v2[1], qts_v2[2], qts_v2[3], qts_v2[4], Inf),
    labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
  ))

  cat("--- Ager Groups (pheno_max_v2) ---\n\n")
  tbl_v2 <- table(prof_df$ager_maxv2[ok_v2])
  for (g in names(tbl_v2))
    cat(sprintf("   %-15s n = %6d (%5.1f%%)\n", g, tbl_v2[g], 100 * tbl_v2[g] / n_v2))

  # Full profiling via shared function
  v2_table <- compute_biomarker_profile(prof_df[ok_v2, ], maxv2_z, usable_bm)

  if (!is.null(v2_table)) {
    sig_v2 <- v2_table[!is.na(v2_table$d_P5_vs_P1) & abs(v2_table$d_P5_vs_P1) > 0.1, ]
    if (nrow(sig_v2) > 0) {
      cat("\n  Biomarkers differentiating fast vs slow agers (pheno_max_v2, |d| > 0.1):\n\n")
      print(sig_v2[, c("biomarker", "n_total", "n_P1", "n_P5", "rho_z", "d_P5_vs_P1",
                        "d_ci_lo", "d_ci_hi", "mean_P1", "mean_P5", "direction", "kw_padj")],
            row.names = FALSE)
    }
    write.csv(v2_table, file.path(OUTPUT_DIR, "step4_maxv2_biomarker_profiles.csv"),
              row.names = FALSE)
    message(">> Saved: step4_maxv2_biomarker_profiles.csv")

    # Heatmap for pheno_max_v2
    n_sig_v2 <- sum(!is.na(v2_table$d_P5_vs_P1) & abs(v2_table$d_P5_vs_P1) > 0.15 &
                    !is.na(v2_table$kw_padj) & v2_table$kw_padj < 0.05)
    heat_v2 <- build_aging_heatmap(v2_table, prof_df[ok_v2, ], maxv2_z)

    if (!is.null(heat_v2)) {
      p_v2 <- heat_v2$plot +
        ggplot2::labs(title = "Accelerated Aging Signature (pheno_max_v2 classifier)",
                      subtitle = sprintf("6/9 Levine biomarkers — n = %s",
                                         format(n_v2, big.mark = ",")))
      ggsave(file.path(OUTPUT_DIR, "step4_maxv2_aging_signature_heatmap.png"),
             p_v2, width = 10, height = max(6, n_sig_v2 * 0.35 + 2), dpi = 150)
      print(p_v2)
      message(">> Saved: step4_maxv2_aging_signature_heatmap.png")

      write.csv(heat_v2$wide, file.path(OUTPUT_DIR, "step4_maxv2_aging_signature.csv"),
                row.names = FALSE)
      cat("\n--- Accelerated Aging Signature (pheno_max_v2) ---\n\n")
      print(heat_v2$wide, row.names = FALSE)
    }
  }
} else {
  message(">> pheno_max_v2 Z-score not found — skipping")
}

# ── 4d-ter. External construct validity by sex ──
#
# Reproduce Cohen's d (P5 vs P1) for males and females separately.
# Flag biomarkers where |d_M - d_F| > 0.2 (clinically meaningful sex interaction).
# Expected candidates: ferritin, hemoglobin, hematocrit, MCHC, transferrin_saturation.

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  STEP 4d-ter: CONSTRUCT VALIDITY BY SEX\n")
cat("══════════════════════════════════════════════════════════════\n\n")

if (maxv2_z %in% names(prof_df) && !is.null(usable_bm) && length(usable_bm) > 0) {

  sex_profiles <- list()
  for (sx in c("F", "M")) {
    idx_sx <- !is.na(prof_df[[maxv2_z]]) & is.finite(prof_df[[maxv2_z]]) &
              !is.na(prof_df$sex) & prof_df$sex == sx
    n_sx <- sum(idx_sx)
    if (n_sx < 100) next
    message(">> Sex-stratified profiling — ", sx, " (n = ", n_sx, ")")
    tbl_sx <- compute_biomarker_profile(prof_df[idx_sx, ], maxv2_z, usable_bm)
    if (!is.null(tbl_sx)) {
      tbl_sx$sex <- sx
      sex_profiles[[sx]] <- tbl_sx
    }
  }

  if (length(sex_profiles) == 2) {
    tbl_F <- sex_profiles[["F"]][, c("biomarker", "d_P5_vs_P1", "d_ci_lo", "d_ci_hi",
                                      "n_P1", "n_P5")]
    tbl_M <- sex_profiles[["M"]][, c("biomarker", "d_P5_vs_P1", "d_ci_lo", "d_ci_hi",
                                      "n_P1", "n_P5")]
    names(tbl_F)[-1] <- paste0(names(tbl_F)[-1], "_F")
    names(tbl_M)[-1] <- paste0(names(tbl_M)[-1], "_M")

    sex_compare <- merge(tbl_F, tbl_M, by = "biomarker", all = TRUE)
    sex_compare$d_diff_MvsF <- sex_compare$d_P5_vs_P1_M - sex_compare$d_P5_vs_P1_F
    sex_compare$sex_flag    <- abs(sex_compare$d_diff_MvsF) > 0.2

    # ── Formal Z-test for sex interaction in Cohen's d ──
    # SE_d via Hedges-Olkin formula; Z = (d_M - d_F) / sqrt(SE_M^2 + SE_F^2)
    # BH correction across all tested biomarkers
    compute_se_d <- function(d, n1, n5) {
      if (any(is.na(c(d, n1, n5))) || n1 < 5 || n5 < 5) return(NA_real_)
      sqrt((n1 + n5) / (n1 * n5) + d^2 / (2 * (n1 + n5)))
    }
    sex_compare$se_d_F <- mapply(compute_se_d,
                                  sex_compare$d_P5_vs_P1_F,
                                  sex_compare$n_P1_F,
                                  sex_compare$n_P5_F)
    sex_compare$se_d_M <- mapply(compute_se_d,
                                  sex_compare$d_P5_vs_P1_M,
                                  sex_compare$n_P1_M,
                                  sex_compare$n_P5_M)
    sex_compare$se_diff  <- sqrt(sex_compare$se_d_F^2 + sex_compare$se_d_M^2)
    sex_compare$Z_interaction <- sex_compare$d_diff_MvsF / sex_compare$se_diff
    sex_compare$p_interaction <- 2 * pnorm(-abs(sex_compare$Z_interaction))

    valid_p <- !is.na(sex_compare$p_interaction)
    sex_compare$p_interaction_adj <- NA_real_
    sex_compare$p_interaction_adj[valid_p] <-
      p.adjust(sex_compare$p_interaction[valid_p], method = "BH")

    sex_compare$sex_interaction_sig <- !is.na(sex_compare$p_interaction_adj) &
                                        sex_compare$p_interaction_adj < 0.05

    sex_compare <- sex_compare[order(-abs(sex_compare$d_diff_MvsF)), ]

    # Print flagged biomarkers (exploratory |Δd|>0.2)
    flagged <- sex_compare[!is.na(sex_compare$sex_flag) & sex_compare$sex_flag, ]
    if (nrow(flagged) > 0) {
      cat(sprintf("  %d biomarkers with |d_M - d_F| > 0.2 (exploratory):\n\n",
                  nrow(flagged)))
      print(flagged[, c("biomarker", "d_P5_vs_P1_F", "d_P5_vs_P1_M",
                        "d_diff_MvsF", "Z_interaction", "p_interaction_adj",
                        "sex_interaction_sig")],
            row.names = FALSE)
    } else {
      cat("  No biomarkers flagged (|d_M - d_F| <= 0.2 for all).\n")
    }

    # Print formally significant interactions (BH-corrected p < 0.05)
    formal_sig <- sex_compare[!is.na(sex_compare$sex_interaction_sig) &
                               sex_compare$sex_interaction_sig, ]
    cat(sprintf("\n  Formally significant sex interactions (BH p < 0.05): %d\n\n",
                nrow(formal_sig)))
    if (nrow(formal_sig) > 0) {
      print(formal_sig[, c("biomarker", "d_P5_vs_P1_F", "d_P5_vs_P1_M",
                           "d_diff_MvsF", "Z_interaction", "p_interaction_adj")],
            row.names = FALSE)
    }

    write.csv(sex_compare,
              file.path(OUTPUT_DIR, "step4_sex_stratified_cohens_d.csv"),
              row.names = FALSE)
    message(">> Saved: step4_sex_stratified_cohens_d.csv")

    # Save per-sex full profile tables
    for (sx in names(sex_profiles)) {
      write.csv(sex_profiles[[sx]],
                file.path(OUTPUT_DIR, paste0("step4_biomarker_profile_sex_", sx, ".csv")),
                row.names = FALSE)
      message(">> Saved: step4_biomarker_profile_sex_", sx, ".csv")
    }

    # Heatmaps per sex
    for (sx in names(sex_profiles)) {
      idx_sx <- !is.na(prof_df[[maxv2_z]]) & is.finite(prof_df[[maxv2_z]]) &
                !is.na(prof_df$sex) & prof_df$sex == sx
      heat_sx <- build_aging_heatmap(sex_profiles[[sx]], prof_df[idx_sx, ], maxv2_z)
      if (!is.null(heat_sx)) {
        n_sig_sx <- sum(!is.na(sex_profiles[[sx]]$d_P5_vs_P1) &
                        abs(sex_profiles[[sx]]$d_P5_vs_P1) > 0.15 &
                        !is.na(sex_profiles[[sx]]$kw_padj) &
                        sex_profiles[[sx]]$kw_padj < 0.05)
        p_sx <- heat_sx$plot +
          ggplot2::labs(
            title    = sprintf("Aging Signature — %s (pheno_max_v2)",
                               ifelse(sx == "F", "Females", "Males")),
            subtitle = sprintf("n = %s", format(sum(idx_sx), big.mark = ","))
          )
        fn_sx <- file.path(OUTPUT_DIR,
                           paste0("step4_maxv2_aging_heatmap_sex_", sx, ".png"))
        ggsave(fn_sx, p_sx,
               width = 10, height = max(6, n_sig_sx * 0.35 + 2), dpi = 150)
        print(p_sx)
        message(">> Saved: ", basename(fn_sx))
      }
    }
  }
} else {
  message(">> Skipping sex-stratified profiling — maxv2_z not found or no usable biomarkers")
}

# ── 4d-quat. External construct validity — exclude clock biomarkers ──
#
# pheno_max_v2 uses: glucose, lncreat (→creatinine), lymph (→lymphocyte),
#                   mcv, rdw, wbc
#
# Including those in the profiling is circular: the clock was trained to
# predict age via these biomarkers, so of course they discriminate
# fast/slow agers best. This section repeats the profiling using ONLY
# biomarkers NOT in the clock, which constitutes true external
# construct validity.
#
# Additionally flagged as "derived" (semi-circular):
#   - nlr        (= neutrophils / lymphocyte — lymphocyte IS in clock)
#   - eGFR       (derived from creatinine  — creatinine IS in clock)
#   - homa_ir    (derived from glucose     — glucose IS in clock)
# These are kept but annotated so reviewers can interpret with nuance.

cat("\n")
cat("══════════════════════════════════════════════════════════════\n")
cat("  STEP 4d-quat: EXTERNAL CONSTRUCT VALIDITY\n")
cat("  (clock biomarkers excluded — true out-of-sample validation)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

if (maxv2_z %in% names(prof_df) && !is.null(usable_bm) && length(usable_bm) > 0) {

  # ── Define sets ──
  clock_bm_direct <- c("rdw", "wbc", "mcv", "glucose", "creatinine", "lymphocyte")
  clock_bm_derived <- c("nlr", "eGFR", "homa_ir")   # derived from clock inputs

  external_bm <- setdiff(usable_bm, clock_bm_direct)
  pure_external_bm <- setdiff(external_bm, clock_bm_derived)

  cat(sprintf("  Clock biomarkers (excluded): %s\n",
              paste(clock_bm_direct, collapse = ", ")))
  cat(sprintf("  Derived biomarkers (kept, flagged): %s\n",
              paste(intersect(clock_bm_derived, external_bm), collapse = ", ")))
  cat(sprintf("  Pure external biomarkers: %d\n\n", length(pure_external_bm)))

  # ── Sensibilidad: distribución de pacientes hematológicos por grupo ──
  if ("flag_hematol" %in% names(prof_df)) {
    n_hematol <- sum(prof_df$flag_hematol, na.rm = TRUE)
    cat(sprintf("\n  Pacientes con posible patología hematológica: %d (%.2f%%)\n",
                n_hematol, 100 * n_hematol / nrow(prof_df)))
    if (n_hematol > 0 && "ager_group" %in% names(prof_df)) {
      tbl_fl <- table(prof_df$ager_group[prof_df$flag_hematol & !is.na(prof_df$flag_hematol)])
      cat("  Distribución en grupos de envejecimiento:\n")
      for (g in names(tbl_fl))
        cat(sprintf("    %-18s %d\n", g, tbl_fl[g]))
    }
  }

  # ── Full external profile (includes derived, flagged) ──
  ext_profile <- compute_biomarker_profile(prof_df, maxv2_z, external_bm)

  if (!is.null(ext_profile) && nrow(ext_profile) > 0) {
    # Flag derived biomarkers
    ext_profile$derived_from_clock <- ext_profile$biomarker %in% clock_bm_derived

    # Significant external biomarkers (|d| > 0.1, FDR < 0.05)
    ext_sig <- ext_profile[!is.na(ext_profile$d_P5_vs_P1) &
                           abs(ext_profile$d_P5_vs_P1) > 0.1 &
                           !is.na(ext_profile$kw_padj) &
                           ext_profile$kw_padj < 0.05, ]

    cat(sprintf("  External biomarkers significantly differentiating fast vs slow agers:\n"))
    cat(sprintf("  (|d| > 0.1, FDR < 0.05, n = %d / %d)\n\n",
                nrow(ext_sig), length(external_bm)))

    print_cols <- c("biomarker", "derived_from_clock", "n_total", "n_P1", "n_P5",
                    "rho_z", "d_P5_vs_P1", "d_ci_lo", "d_ci_hi",
                    "mean_P1", "mean_P5", "direction", "kw_padj")
    print_cols <- intersect(print_cols, names(ext_sig))
    print(ext_sig[, print_cols], row.names = FALSE)

    # Summary stats
    n_pure_sig <- sum(ext_sig$biomarker %in% pure_external_bm)
    n_derived_sig <- sum(ext_sig$biomarker %in% clock_bm_derived)
    cat(sprintf("\n  Pure external significant: %d\n", n_pure_sig))
    cat(sprintf("  Derived (semi-circular) significant: %d\n\n", n_derived_sig))

    # Save
    ext_profile$clock <- "pheno_max_v2"
    write.csv(ext_profile,
              file.path(OUTPUT_DIR, "step4_external_construct_validity.csv"),
              row.names = FALSE)
    message(">> Saved: step4_external_construct_validity.csv")

    # ── Análisis de sensibilidad: excluir pacientes con posible patología hematológica ──
    #
    # Objetivo: verificar si el efecto de linfocitos/WBC en envejecedores rápidos
    # se mantiene al excluir pacientes que podrían tener leucemia/linfoma.
    # Biomarcadores de interés: lymphocyte_abs, lymphocyte (%), wbc, neutrophils, nlr.

    if ("flag_hematol" %in% names(prof_df)) {
      n_hematol_prof <- sum(prof_df$flag_hematol, na.rm = TRUE)
      if (n_hematol_prof > 0) {
        cat(sprintf("\n  ── Sensibilidad: excluir %d pacientes con flag_hematol ──\n\n",
                    n_hematol_prof))
        prof_df_clean <- prof_df[!prof_df$flag_hematol | is.na(prof_df$flag_hematol), ]
        ext_profile_clean <- compute_biomarker_profile(prof_df_clean, maxv2_z, external_bm)

        if (!is.null(ext_profile_clean) && nrow(ext_profile_clean) > 0) {
          ext_profile_clean$derived_from_clock <- ext_profile_clean$biomarker %in% clock_bm_derived
          ext_profile_clean$clock <- "pheno_max_v2_clean"
          write.csv(ext_profile_clean,
                    file.path(OUTPUT_DIR, "step4_external_validity_sensitivity.csv"),
                    row.names = FALSE)
          message(">> Saved: step4_external_validity_sensitivity.csv")

          # Comparar efecto de biomarcadores hematológicos con/sin flagueados
          bm_hematol_check <- c("lymphocyte_abs", "lymphocyte", "wbc", "neutrophils", "nlr",
                                 "basophils", "monocytes", "eosinophils")
          bm_present <- intersect(bm_hematol_check, ext_profile$biomarker)
          if (length(bm_present) > 0) {
            cat(sprintf("  %-24s  %8s  %8s  %8s  %s\n",
                        "biomarker", "d_total", "d_clean", "delta", "direction"))
            for (bm in bm_present) {
              row_all   <- ext_profile[ext_profile$biomarker == bm, ]
              row_clean <- ext_profile_clean[ext_profile_clean$biomarker == bm, ]
              if (nrow(row_all) > 0 && nrow(row_clean) > 0) {
                d_all   <- row_all$d_P5_vs_P1
                d_clean <- row_clean$d_P5_vs_P1
                delta   <- d_clean - d_all
                flag_change <- if (abs(delta) > 0.10) " *** CAMBIO" else ""
                cat(sprintf("  %-24s  %+7.3f  %+7.3f  %+7.3f%s\n",
                            bm, d_all, d_clean, delta, flag_change))
              }
            }
          }
        }
      } else {
        message(">> No hay pacientes flagueados — análisis de sensibilidad no aplica.")
      }
    }

    # Heatmap — pure external only
    pure_ext_profile <- ext_profile[ext_profile$biomarker %in% pure_external_bm, ]
    heat_ext <- build_aging_heatmap(pure_ext_profile, prof_df, maxv2_z, d_thresh = 0.15)
    if (!is.null(heat_ext)) {
      n_sig_ext <- sum(!is.na(pure_ext_profile$d_P5_vs_P1) &
                       abs(pure_ext_profile$d_P5_vs_P1) > 0.15 &
                       !is.na(pure_ext_profile$kw_padj) &
                       pure_ext_profile$kw_padj < 0.05)
      p_ext <- heat_ext$plot +
        ggplot2::labs(
          title    = "External Construct Validity — pheno_max_v2",
          subtitle = sprintf("Clock biomarkers excluded (rdw, wbc, mcv, glucose, creatinine, lymphocyte) | n = %s",
                             format(sum(!is.na(prof_df[[maxv2_z]])), big.mark = ","))
        )
      ggsave(file.path(OUTPUT_DIR, "step4_external_construct_heatmap.png"),
             p_ext,
             width = 10,
             height = max(6, n_sig_ext * 0.35 + 2),
             dpi = 150)
      print(p_ext)
      message(">> Saved: step4_external_construct_heatmap.png")
    }

    # Sex-stratified external construct validity
    cat("\n  External construct validity by sex:\n\n")
    sex_ext_profiles <- list()
    for (sx in c("F", "M")) {
      idx_sx <- !is.na(prof_df[[maxv2_z]]) & is.finite(prof_df[[maxv2_z]]) &
                !is.na(prof_df$sex) & prof_df$sex == sx
      if (sum(idx_sx) < 100) next
      tbl_ext_sx <- compute_biomarker_profile(prof_df[idx_sx, ], maxv2_z, pure_external_bm)
      if (!is.null(tbl_ext_sx)) {
        tbl_ext_sx$sex <- sx
        sex_ext_profiles[[sx]] <- tbl_ext_sx
        n_sig_sx_ext <- sum(!is.na(tbl_ext_sx$d_P5_vs_P1) & abs(tbl_ext_sx$d_P5_vs_P1) > 0.1 &
                            !is.na(tbl_ext_sx$kw_padj) & tbl_ext_sx$kw_padj < 0.05)
        cat(sprintf("    %s (n=%s): %d pure-external biomarkers |d|>0.1, FDR<0.05\n",
                    sx, format(sum(idx_sx), big.mark = ","), n_sig_sx_ext))
        top5 <- head(tbl_ext_sx[order(-abs(tbl_ext_sx$d_P5_vs_P1 %/% 1 * 1000 +
                                           abs(tbl_ext_sx$d_P5_vs_P1))), ], 5)
        top5 <- head(tbl_ext_sx[order(-abs(tbl_ext_sx$d_P5_vs_P1)), ], 5)
        cat(sprintf("      Top 5: %s\n",
                    paste(sprintf("%s(d=%+.2f)", top5$biomarker, top5$d_P5_vs_P1),
                          collapse = ", ")))
      }
    }

    if (length(sex_ext_profiles) == 2) {
      sex_ext_all <- do.call(rbind, sex_ext_profiles)
      write.csv(sex_ext_all,
                file.path(OUTPUT_DIR, "step4_external_validity_by_sex.csv"),
                row.names = FALSE)
      message(">> Saved: step4_external_validity_by_sex.csv")

      # ── ICC: replication of external d-profile between sexes ──
      #
      # ICC(2,1) measures whether the SAME biomarkers consistently
      # differentiate fast vs slow agers in both females and males.
      # ICC consistency = rank agreement ignoring mean sex shift in d.
      # ICC agreement   = penalises mean sex shift (stricter criterion).
      # 95% CI via F-distribution (Shrout & Fleiss 1979).

      cat("\n  --- ICC: External Construct Validity Profile (F vs M) ---\n\n")

      df_F  <- sex_ext_profiles[["F"]][, c("biomarker", "d_P5_vs_P1", "n_P1", "n_P5")]
      df_M  <- sex_ext_profiles[["M"]][, c("biomarker", "d_P5_vs_P1", "n_P1", "n_P5")]
      names(df_F)[-1] <- paste0(names(df_F)[-1], "_F")
      names(df_M)[-1] <- paste0(names(df_M)[-1], "_M")

      d_icc <- merge(df_F, df_M, by = "biomarker")
      d_icc <- d_icc[complete.cases(d_icc[, c("d_P5_vs_P1_F", "d_P5_vs_P1_M")]), ]
      n_bm  <- nrow(d_icc)
      k     <- 2

      mat      <- as.matrix(d_icc[, c("d_P5_vs_P1_F", "d_P5_vs_P1_M")])
      row_m    <- rowMeans(mat)
      col_m    <- colMeans(mat)
      grand_m  <- mean(mat)

      SS_bm    <- k    * sum((row_m  - grand_m)^2)  # between biomarkers
      SS_sex   <- n_bm * sum((col_m  - grand_m)^2)  # between sexes
      SS_total <- sum((mat - grand_m)^2)
      SS_error <- SS_total - SS_bm - SS_sex

      MS_bm    <- SS_bm    / (n_bm - 1)
      MS_sex   <- SS_sex   / (k - 1)
      MS_error <- SS_error / ((n_bm - 1) * (k - 1))

      # ICC(2,1) consistency — does NOT penalise mean d shift between sexes
      icc_c <- (MS_bm - MS_error) / (MS_bm + (k - 1) * MS_error)

      # ICC(2,1) agreement — DOES penalise mean d shift between sexes
      denom_a <- MS_bm + (k - 1) * MS_error + k * (MS_sex - MS_error) / n_bm
      icc_a   <- (MS_bm - MS_error) / denom_a

      # 95% CI for ICC consistency via F-distribution
      F_val  <- MS_bm / MS_error
      df1    <- n_bm - 1
      df2    <- (n_bm - 1) * (k - 1)
      FL     <- F_val / qf(0.975, df1, df2)
      FU     <- F_val * qf(0.975, df2, df1)
      ci_lo  <- (FL - 1) / (FL + k - 1)
      ci_hi  <- (FU - 1) / (FU + k - 1)

      spear_test <- cor.test(d_icc$d_P5_vs_P1_F, d_icc$d_P5_vs_P1_M,
                             method = "spearman", exact = FALSE)

      cat(sprintf("    n biomarkers: %d\n", n_bm))
      cat(sprintf("    ICC(2,1) consistency: %.3f  [95%% CI: %.3f \u2013 %.3f]\n",
                  icc_c, ci_lo, ci_hi))
      cat(sprintf("    ICC(2,1) agreement:   %.3f\n", icc_a))
      cat(sprintf("    Spearman rho (rank):  %.3f  (p = %.2e)\n",
                  spear_test$estimate, spear_test$p.value))
      cat(sprintf("    Mean d_F: %+.3f  |  Mean d_M: %+.3f  |  \u0394 = %+.3f\n",
                  mean(d_icc$d_P5_vs_P1_F), mean(d_icc$d_P5_vs_P1_M),
                  mean(d_icc$d_P5_vs_P1_M) - mean(d_icc$d_P5_vs_P1_F)))
      cat(sprintf("    MS_between_bm: %.4f  |  MS_sex: %.4f  |  MS_error: %.4f\n",
                  MS_bm, MS_sex, MS_error))

      # Interpretation guide
      icc_label <- dplyr::case_when(
        icc_c >= 0.90 ~ "Excellent",
        icc_c >= 0.75 ~ "Good",
        icc_c >= 0.50 ~ "Moderate",
        TRUE          ~ "Poor"
      )
      cat(sprintf("    Interpretation (ICC consistency): %s\n\n", icc_label))

      # Guardar el resumen ICC en CSV — antes sólo se imprimía a consola/log,
      # así que extract_results_for_paper.R nunca podía citar el número real
      # (el manuscrito quedaba con placeholders tipo "reported from model
      # output"). Ahora queda disponible para el paso 5 del pipeline.
      icc_summary <- data.frame(
        clock            = "pheno_pheno_max_v2",
        n_biomarkers     = n_bm,
        icc_consistency  = round(icc_c, 4),
        icc_ci_lo        = round(ci_lo, 4),
        icc_ci_hi        = round(ci_hi, 4),
        icc_agreement    = round(icc_a, 4),
        spearman_rho     = round(as.numeric(spear_test$estimate), 4),
        spearman_p       = spear_test$p.value,
        mean_d_F         = round(mean(d_icc$d_P5_vs_P1_F), 4),
        mean_d_M         = round(mean(d_icc$d_P5_vs_P1_M), 4),
        icc_interpretation = icc_label,
        stringsAsFactors = FALSE
      )
      write.csv(icc_summary,
                file.path(OUTPUT_DIR, "step4_icc_summary.csv"),
                row.names = FALSE)
      message(">> Saved: step4_icc_summary.csv")

      # Per-biomarker agreement: (d_F + d_M)/2 ± abs(d_F - d_M)
      d_icc$d_mean  <- (d_icc$d_P5_vs_P1_F + d_icc$d_P5_vs_P1_M) / 2
      d_icc$d_diff  <- d_icc$d_P5_vs_P1_M - d_icc$d_P5_vs_P1_F
      d_icc$abs_diff <- abs(d_icc$d_diff)
      d_icc <- d_icc[order(-abs(d_icc$d_mean)), ]

      cat("    Per-biomarker (sorted by |mean d|):\n")
      cat(sprintf("      %-24s  %7s  %7s  %7s  %7s\n",
                  "biomarker", "d_F", "d_M", "mean_d", "|diff|"))
      for (i in seq_len(min(nrow(d_icc), 20))) {
        rw <- d_icc[i, ]
        cat(sprintf("      %-24s  %+6.3f   %+6.3f   %+6.3f   %.3f\n",
                    rw$biomarker, rw$d_P5_vs_P1_F, rw$d_P5_vs_P1_M,
                    rw$d_mean, rw$abs_diff))
      }

      write.csv(d_icc,
                file.path(OUTPUT_DIR, "step4_icc_external_profile_FM.csv"),
                row.names = FALSE)
      message(">> Saved: step4_icc_external_profile_FM.csv")
    }

  } else {
    message(">> No external biomarkers available for construct validity")
  }

} else {
  message(">> Skipping external construct validity — maxv2_z not found")
}

# ── 4e. Biomarker profiles by age decade and sex ──

cat("\n--- Biomarker Trends by Age Decade & Sex ---\n\n")

prof_df$age_decade <- cut(prof_df$age,
  breaks = c(0, 30, 40, 50, 60, 70, 80, 120),
  labels = c("<30", "30s", "40s", "50s", "60s", "70s", "80+"),
  right = FALSE
)

# Select key biomarkers that differentiate fast/slow agers
key_bm <- character(0)
if (!is.null(profile_table) && nrow(profile_table) > 0) {
  key_bm <- profile_table$biomarker[
    !is.na(profile_table$d_P5_vs_P1) &
    abs(profile_table$d_P5_vs_P1) > 0.15 &
    !is.na(profile_table$kw_padj) &
    profile_table$kw_padj < 0.05
  ]
  key_bm <- head(key_bm, 12)  # Top 12 for readability
}

if (length(key_bm) >= 3) {
  for (bm in key_bm) {
    if (!bm %in% names(prof_df)) next
    x <- suppressWarnings(as.numeric(prof_df[[bm]]))
    ok <- !is.na(x) & is.finite(x) & !is.na(prof_df$age_decade) &
          prof_df$sex %in% c("F", "M")
    if (sum(ok) < 100) next

    plot_df <- data.frame(
      value = x[ok],
      age_decade = prof_df$age_decade[ok],
      sex = prof_df$sex[ok]
    )

    p_bm <- ggplot(plot_df, aes(x = age_decade, y = value, fill = sex)) +
      geom_boxplot(alpha = 0.6, outlier.size = 0.3, outlier.alpha = 0.1) +
      scale_fill_manual(values = c("F" = "firebrick", "M" = "steelblue")) +
      theme_minimal(base_size = 12) +
      labs(title = paste0("Biomarker by Age & Sex: ", bm),
           x = "Age Decade", y = bm, fill = "Sex")
    print(p_bm)

    ggsave(file.path(OUTPUT_DIR, paste0("step4_biomarker_age_sex_", bm, ".png")),
           p_bm, width = 9, height = 5, dpi = 150)
  }
  message(">> Saved biomarker × age × sex plots for ", length(key_bm), " biomarkers")
}

# ── 4f. Advancement by age decade (recentered) for all clocks ──

cat("\n--- Advancement Distribution by Age Decade (Recentered) ---\n\n")

for (zc in z_cols) {
  clock_label <- sub("_sexadj_z$", "", zc)
  short_label <- gsub("pheno_", "", clock_label)

  ok <- !is.na(rc[[zc]]) & is.finite(rc[[zc]]) & !is.na(rc$age_decade)
  if (sum(ok) < 100) next

  p_box <- ggplot(rc[ok, ], aes(x = age_decade, y = .data[[zc]])) +
    geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    theme_minimal(base_size = 13) +
    labs(title = paste0("Recentered Z-Score by Age: ", short_label),
         subtitle = "After sex-adjusted residualization",
         x = "Age Group", y = "Z-score (0 = population average)")
  print(p_box)

  ggsave(file.path(OUTPUT_DIR, paste0("step4_zscore_by_age_", short_label, ".png")),
         p_box, width = 8, height = 5, dpi = 150)
}


# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  PIPELINE COMPLETE                                      ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

cat("  STEP 1 — NHANES Training\n")
if (exists("nhanes_table") && !is.null(nhanes_table)) {
  for (i in seq_len(nrow(nhanes_table))) {
    ct <- nhanes_table[i, ]
    cat(sprintf("    %-35s r = %.4f  RMSE = %.2f  (n = %s)\n",
                ct$clock, ct$r, ct$RMSE, format(ct$n, big.mark = ",")))
  }
}

cat("\n  STEP 2 — Argentine Scoring\n")
cat(sprintf("    Cohort: %s patients (%.0f%% F, age %.0f ± %.0f)\n",
            format(nrow(results), big.mark = ","),
            100 * mean(results$sex == "F", na.rm = TRUE),
            mean(results$age, na.rm = TRUE),
            sd(results$age, na.rm = TRUE)))
if (exists("arg_table") && !is.null(arg_table)) {
  for (i in seq_len(nrow(arg_table))) {
    ct <- arg_table[i, ]
    cat(sprintf("    %-35s r = %.4f  RMSE = %.2f  n = %s\n",
                ct$clock, ct$r, ct$RMSE, format(ct$n, big.mark = ",")))
  }
}

cat("\n  STEP 3 — Recentering\n")
cat(sprintf("    Clocks recentered: %d\n", length(viable_adv)))
cat("    Methods: A (age-residual), B (age+sex), C (age-stratified)\n")

cat("\n  STEP 4 — Fast/Slow Ager Profiling\n")
cat(sprintf("    Subjects profiled: %s\n", format(n_valid, big.mark = ",")))
cat(sprintf("    Primary Z-score: %s\n", primary_z))
if (!is.null(profile_table)) {
  n_sig <- sum(!is.na(profile_table$d_P5_vs_P1) &
               abs(profile_table$d_P5_vs_P1) > 0.15 &
               !is.na(profile_table$kw_padj) &
               profile_table$kw_padj < 0.05)
  cat(sprintf("    Significant biomarkers (|d| > 0.15): %d / %d\n",
              n_sig, nrow(profile_table)))
}

cat("\n  Output directory: ", OUTPUT_DIR, "\n")
cat("  Files saved:\n")
saved <- list.files(OUTPUT_DIR, full.names = FALSE)
for (f in saved) cat("    - ", f, "\n")

cat("\n  Key data files:\n")
cat("    - data_derived/bioage_deployment_bundle.rds  (trained models)\n")
cat("    - data_derived/clinical_scored_lite.csv      (raw + recentered scores)\n")
cat("    - final_output/step4_aging_signature_heatmap.png\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║  Done. Results ready for manuscript.                    ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")

# Store everything
final_results <- list(
  nhanes_perf     = if (exists("nhanes_table")) nhanes_table else NULL,
  argentina_perf  = if (exists("arg_table")) arg_table else NULL,
  scored_data     = results,
  recentered_data = rc,
  profile_table   = profile_table,
  bundle          = bundle,
  primary_z       = primary_z
)
message(">> All results in 'final_results'")
