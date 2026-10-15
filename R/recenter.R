###############################################################################
# recenter.R - local reference for the clock (Method B)
#
# Fit lm(BA ~ age + sex) in the local cohort. The residual is the distance of a
# patient from the expectation for their age and sex; dividing by the residual
# SD gives a Z score. A positive Z means a higher biological age than expected
# for age and sex. It is a position among local peers, not a measured rate of
# aging, and it is not a recalibration of mortality risk.
###############################################################################

recenter_method_b <- function(df, ba) {
  ok  <- is.finite(ba) & !is.na(df$age) & df$sex %in% c("F", "M")
  fit <- lm(ba[ok] ~ df$age[ok] + df$sex[ok])
  cf  <- unname(coef(fit))
  res <- residuals(fit)

  z <- rep(NA_real_, nrow(df))
  z[ok] <- as.numeric(scale(res))

  # Extreme standardized residuals are set to missing; they are logged
  # (patient level, kept out of version control) so they can be reviewed.
  extreme <- !is.na(z) & abs(z) > Z_CENSOR
  fh <- df$flag_hematol[extreme]
  z_clean <- z
  z_clean[extreme] <- NA_real_

  zc   <- z_clean[!is.na(z_clean)]
  qts  <- quantile(zc, probs = c(0.10, 0.25, 0.75, 0.90))
  list(
    z = z_clean, extreme = extreme,
    coefs = list(intercept = cf[1], slope_age = cf[2], slope_male = cf[3],
                 resid_sd = sd(res), n_fit = sum(ok)),
    cutpoints = list(p10 = unname(qts[1]), p25 = unname(qts[2]), p50 = median(zc),
                     p75 = unname(qts[3]), p90 = unname(qts[4])),
    summary = data.frame(
      clock = CLOCK_NAME, n = length(zc),
      z_min = min(zc), z_max = max(zc), z_iqr = IQR(zc),
      z_p01 = unname(quantile(zc, 0.01)), z_p99 = unname(quantile(zc, 0.99)),
      n_gt3 = sum(abs(zc) > 3), n_gt4 = sum(abs(zc) > 4),
      n_censored = sum(extreme),
      n_censored_explained_hematol = sum(fh, na.rm = TRUE),
      p10_empirical = round(unname(qts[1]), 4), p25_empirical = round(unname(qts[2]), 4),
      p50_empirical = round(median(zc), 4),
      p75_empirical = round(unname(qts[3]), 4), p90_empirical = round(unname(qts[4]), 4)
    )
  )
}

# Five bands from the empirical percentiles of Z: P1 lowest 10%, P2 next 15%,
# P3 middle 50%, P4 next 15%, P5 highest 10%. Percentiles are taken from the
# observed distribution; normality is not assumed.
BAND_LABELS <- c("P1", "P2", "P3", "P4", "P5")

assign_bands <- function(z) {
  ok <- is.finite(z)
  q  <- quantile(z[ok], probs = c(0.10, 0.25, 0.75, 0.90))
  g  <- rep(NA_character_, length(z))
  g[ok] <- as.character(cut(z[ok], breaks = c(-Inf, q, Inf), labels = BAND_LABELS))
  g
}
