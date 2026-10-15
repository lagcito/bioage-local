###############################################################################
# validate.R - external construct validity of the local score
#
# Without local outcomes, the question is whether the score is associated, in
# the expected direction, with biomarkers the clock never used. P5 (highest
# decile) is compared with P1 (lowest decile): Cohen d with a 95% CI, Spearman
# correlation over the whole range, and a Kruskal-Wallis test across the five
# bands. P-values are Benjamini-Hochberg adjusted over the markers analysed
# together. Because P1 and P5 are extreme groups, d is an upper bound of the
# effect, not its typical size; the Spearman correlation does not select
# extremes.
###############################################################################

profile_by_band <- function(data, z_col, biomarkers, min_n = 30) {
  g <- assign_bands(data[[z_col]])
  rows <- lapply(biomarkers, function(bm) {
    if (!bm %in% names(data)) return(NULL)
    x  <- data[[bm]]
    ok <- is.finite(x) & !is.na(g)
    if (sum(ok) < min_n) return(NULL)
    p1 <- x[ok & g == "P1"]; p5 <- x[ok & g == "P5"]
    n1 <- length(p1); n5 <- length(p5)

    d <- lo <- hi <- NA_real_
    if (n1 >= 5 && n5 >= 5) {
      sp <- sqrt((var(p1) * (n1 - 1) + var(p5) * (n5 - 1)) / (n1 + n5 - 2))
      if (sp > 0) {
        d  <- (mean(p5) - mean(p1)) / sp
        se <- sqrt((n1 + n5) / (n1 * n5) + d^2 / (2 * (n1 + n5)))
        lo <- d - 1.96 * se; hi <- d + 1.96 * se
      }
    }
    kw <- tryCatch(kruskal.test(x[ok] ~ g[ok]), error = function(e) NULL)
    ct <- tryCatch(cor.test(data[[z_col]][ok], x[ok], method = "spearman", exact = FALSE),
                   error = function(e) NULL)
    data.frame(
      biomarker = bm, n_total = sum(ok), n_P1 = n1, n_P5 = n5,
      rho_z = if (!is.null(ct)) round(unname(ct$estimate), 4) else NA_real_,
      rho_p = if (!is.null(ct)) ct$p.value else NA_real_,
      kw_chi2 = if (!is.null(kw)) round(unname(kw$statistic), 2) else NA_real_,
      kw_p = if (!is.null(kw)) kw$p.value else NA_real_,
      d_P5_vs_P1 = round(d, 3), d_ci_lo = round(lo, 3), d_ci_hi = round(hi, 3),
      mean_P1 = if (n1 >= 5) round(mean(p1), 3) else NA_real_,
      mean_P5 = if (n5 >= 5) round(mean(p5), 3) else NA_real_,
      direction = ifelse(is.na(d), NA_character_,
                  ifelse(d > EFFECT_THRESHOLD, "HIGHER_in_P5",
                  ifelse(d < -EFFECT_THRESHOLD, "LOWER_in_P5", "~equal"))),
      stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  tbl$rho_padj <- p.adjust(tbl$rho_p, method = "BH")
  tbl$kw_padj  <- p.adjust(tbl$kw_p,  method = "BH")
  rownames(tbl) <- NULL
  tbl[order(-abs(tbl$d_P5_vs_P1)), ]
}

meets_criterion <- function(tbl) {
  !is.na(tbl$d_P5_vs_P1) & abs(tbl$d_P5_vs_P1) > EFFECT_THRESHOLD &
    !is.na(tbl$kw_padj) & tbl$kw_padj < FDR_THRESHOLD
}

# Tests ordered because a clinician is concerned can inflate an effect. For
# each marker: coverage among patients with a score, how many were measured in
# P5 and in P1, and a flag when coverage is below 50% and the P5/P1 ratio above 1.5.
ascertainment_table <- function(tbl, n_scored) {
  ratio <- ifelse(!is.na(tbl$n_P1) & tbl$n_P1 > 0, round(tbl$n_P5 / tbl$n_P1, 2), NA_real_)
  cover <- round(100 * tbl$n_total / n_scored, 1)
  data.frame(biomarker = tbl$biomarker, coverage = cover,
             n_P1 = tbl$n_P1, n_P5 = tbl$n_P5, ratio_P5_P1 = ratio,
             d_P5_vs_P1 = tbl$d_P5_vs_P1,
             ascertainment_flag = !is.na(cover) & !is.na(ratio) & cover < 50 & ratio > 1.5)
}

# frame: one row per patient with columns z, sex, flag_hematol and the candidate markers.
external_validity <- function(frame) {
  usable <- EXTERNAL_CANDIDATES[vapply(EXTERNAL_CANDIDATES, function(b)
    b %in% names(frame) && sum(is.finite(frame[[b]])) >= 50, logical(1))]
  external <- setdiff(usable, CLOCK_INPUT_COLUMNS)
  pure     <- setdiff(external, DERIVED_COLUMNS)

  main <- profile_by_band(frame, "z", external)
  main$derived_from_clock <- main$biomarker %in% DERIVED_COLUMNS
  main$clock <- CLOCK_NAME

  clean <- frame[!frame$flag_hematol | is.na(frame$flag_hematol), ]
  sens  <- profile_by_band(clean, "z", external)
  sens$derived_from_clock <- sens$biomarker %in% DERIVED_COLUMNS

  by_sex <- do.call(rbind, lapply(c("F", "M"), function(sx) {
    sub <- frame[is.finite(frame$z) & !is.na(frame$sex) & frame$sex == sx, ]
    if (nrow(sub) < 100) return(NULL)
    t <- profile_by_band(sub, "z", pure); t$sex <- sx; t
  }))

  list(external = main, sensitivity = sens, by_sex = by_sex,
       ascertainment = ascertainment_table(main, sum(is.finite(frame$z))),
       pure = pure)
}
