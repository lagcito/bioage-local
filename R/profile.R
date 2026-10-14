# R/profile.R
# Biomarker profiling: fast vs slow ager classification and Cohen's d analysis.
#
# compute_biomarker_profile(data, z_col, biomarkers, min_n)
#   - Classifies subjects into 5 ager groups via empirical quintiles of z_col
#   - For each biomarker: Spearman rho, Kruskal-Wallis, Cohen's d (P5 vs P1) with 95% CI
#   - Returns data.frame sorted by |d|, with BH-adjusted p-values
#
# build_aging_heatmap(profile_table, data, z_col, d_thresh, fdr_thresh)
#   - Filters to significant biomarkers and builds a ggplot2 tile heatmap
#   - Returns list: $plot, $data (long), $wide (pivot_wider)

compute_biomarker_profile <- function(data, z_col, biomarkers, min_n = 30) {
  stopifnot(z_col %in% names(data))

  ok   <- !is.na(data[[z_col]]) & is.finite(data[[z_col]])
  qts  <- quantile(data[[z_col]][ok], probs = c(0.10, 0.25, 0.75, 0.90))

  data$ager_group <- NA_character_
  data$ager_group[ok] <- as.character(cut(
    data[[z_col]][ok],
    breaks = c(-Inf, qts[1], qts[2], qts[3], qts[4], Inf),
    labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
  ))

  results <- lapply(biomarkers, function(bm) {
    if (!bm %in% names(data)) return(NULL)
    x  <- data[[bm]]
    g  <- data$ager_group
    ok_bm <- !is.na(x) & is.finite(x) & !is.na(g)
    if (sum(ok_bm) < min_n) return(NULL)

    p1 <- x[ok_bm & g == "P1_very_slow"]
    p5 <- x[ok_bm & g == "P5_very_fast"]
    n1 <- length(p1); n5 <- length(p5)

    cohens_d <- d_ci_lo <- d_ci_hi <- NA_real_
    if (n1 >= 5 && n5 >= 5) {
      ps <- sqrt((var(p1) * (n1 - 1) + var(p5) * (n5 - 1)) / (n1 + n5 - 2))
      if (ps > 0) {
        cohens_d <- (mean(p5) - mean(p1)) / ps
        se_d     <- sqrt((n1 + n5) / (n1 * n5) + cohens_d^2 / (2 * (n1 + n5)))
        d_ci_lo  <- cohens_d - 1.96 * se_d
        d_ci_hi  <- cohens_d + 1.96 * se_d
      }
    }

    kw <- tryCatch(kruskal.test(x[ok_bm] ~ g[ok_bm]), error = function(e) NULL)
    ct <- tryCatch(
      cor.test(data[[z_col]][ok_bm], x[ok_bm], method = "spearman", exact = FALSE),
      error = function(e) NULL
    )

    data.frame(
      biomarker  = bm,
      n_total    = sum(ok_bm),
      n_P1       = n1,
      n_P5       = n5,
      rho_z      = if (!is.null(ct)) round(ct$estimate, 4) else NA_real_,
      rho_p      = if (!is.null(ct)) ct$p.value else NA_real_,
      kw_chi2    = if (!is.null(kw)) round(kw$statistic, 2) else NA_real_,
      kw_p       = if (!is.null(kw)) kw$p.value else NA_real_,
      d_P5_vs_P1 = round(cohens_d, 3),
      d_ci_lo    = round(d_ci_lo, 3),
      d_ci_hi    = round(d_ci_hi, 3),
      mean_P1    = if (n1 >= 5) round(mean(p1), 3) else NA_real_,
      mean_P5    = if (n5 >= 5) round(mean(p5), 3) else NA_real_,
      direction  = ifelse(!is.na(cohens_d),
                   ifelse(cohens_d >  0.1, "HIGHER_in_fast",
                   ifelse(cohens_d < -0.1, "LOWER_in_fast", "~equal")),
                   NA_character_),
      stringsAsFactors = FALSE
    )
  })

  tbl <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)

  tbl$rho_padj <- p.adjust(tbl$rho_p, method = "BH")
  tbl$kw_padj  <- p.adjust(tbl$kw_p,  method = "BH")
  rownames(tbl) <- NULL
  tbl[order(-abs(tbl$d_P5_vs_P1)), ]
}


build_aging_heatmap <- function(profile_table, data, z_col,
                                d_thresh = 0.15, fdr_thresh = 0.05) {
  sig_bm <- profile_table$biomarker[
    !is.na(profile_table$d_P5_vs_P1) &
    abs(profile_table$d_P5_vs_P1) > d_thresh &
    !is.na(profile_table$kw_padj) &
    profile_table$kw_padj < fdr_thresh
  ]
  if (length(sig_bm) < 3) return(NULL)

  # Recompute ager_group from z_col
  ok  <- !is.na(data[[z_col]]) & is.finite(data[[z_col]])
  qts <- quantile(data[[z_col]][ok], probs = c(0.10, 0.25, 0.75, 0.90))
  data$ager_group <- NA_character_
  data$ager_group[ok] <- as.character(cut(
    data[[z_col]][ok],
    breaks = c(-Inf, qts[1], qts[2], qts[3], qts[4], Inf),
    labels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast")
  ))

  radar_rows <- list()
  for (bm in sig_bm) {
    if (!bm %in% names(data)) next
    x     <- data[[bm]]
    ok_bm <- !is.na(x) & is.finite(x) & !is.na(data$ager_group)
    if (sum(ok_bm) < 30) next
    x_z   <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    gm    <- tapply(x_z[ok_bm], data$ager_group[ok_bm], mean)
    for (g in names(gm)) {
      radar_rows[[length(radar_rows) + 1]] <- data.frame(
        biomarker = bm, group = g, mean_z = round(gm[g], 3),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(radar_rows) == 0) return(NULL)

  long_df <- do.call(rbind, radar_rows)
  long_df$group <- factor(long_df$group,
    levels = c("P1_very_slow", "P2_slow", "P3_average", "P4_fast", "P5_very_fast"))

  bm_order <- long_df[long_df$group %in% c("P1_very_slow", "P5_very_fast"), ]
  bm_order <- tidyr::pivot_wider(bm_order, names_from = group, values_from = mean_z)
  bm_order <- bm_order[complete.cases(bm_order), ]
  bm_order$diff <- bm_order$P5_very_fast - bm_order$P1_very_slow
  bm_order <- bm_order[order(bm_order$diff), ]
  long_df$biomarker <- factor(long_df$biomarker, levels = bm_order$biomarker)

  wide_df <- tidyr::pivot_wider(long_df, names_from = group, values_from = mean_z)

  p <- ggplot2::ggplot(long_df,
         ggplot2::aes(x = group, y = biomarker, fill = mean_z)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_z)), size = 3) +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                                  midpoint = 0, name = "Mean Z") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "Ager Group", y = "Biomarker")

  list(plot = p, data = long_df, wide = as.data.frame(wide_df))
}
