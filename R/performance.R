###############################################################################
# performance.R - how the clock behaves in the local cohort
###############################################################################

AGE_BREAKS <- c(0, 30, 40, 50, 60, 70, 80, 120)
AGE_LABELS <- c("<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

rmse_by_age_group <- function(ba, age) {
  grp <- cut(age, breaks = AGE_BREAKS, labels = AGE_LABELS, right = FALSE)
  rows <- lapply(AGE_LABELS, function(g) {
    ok <- is.finite(ba) & is.finite(age) & !is.na(grp) & grp == g
    if (sum(ok) < 10) return(NULL)
    d <- ba[ok] - age[ok]
    data.frame(clock = CLOCK_NAME, age_group = g, n = sum(ok),
               RMSE = round(sqrt(mean(d^2)), 2), MAE = round(mean(abs(d)), 2),
               bias = round(mean(d), 2), r = round(cor(age[ok], ba[ok]), 4))
  })
  do.call(rbind, rows)
}

# Bland-Altman agreement between BA and CA: parametric and empirical limits.
bland_altman <- function(ba, age) {
  ok <- is.finite(ba) & is.finite(age)
  d  <- ba[ok] - age[ok]
  data.frame(
    clock = CLOCK_NAME, n = sum(ok),
    bias = round(mean(d), 3), sd_diff = round(sd(d), 3),
    loa_lo = round(mean(d) - 1.96 * sd(d), 3), loa_hi = round(mean(d) + 1.96 * sd(d), 3),
    loa_lo_empirical = round(as.numeric(quantile(d, 0.025)), 3),
    loa_hi_empirical = round(as.numeric(quantile(d, 0.975)), 3),
    pct_within_5yr = round(100 * mean(abs(d) <= 5), 1)
  )
}

# Are the scored patients different from the unscored ones in age and sex?
# Descriptive only (standardized mean difference), no significance tests.
scored_vs_unscored <- function(df, scored) {
  s <- df[scored, ]; u <- df[!scored, ]
  smd <- (mean(s$age) - mean(u$age)) / sqrt((var(s$age) + var(u$age)) / 2)
  data.frame(
    group      = c("scored", "unscored"),
    n          = c(nrow(s), nrow(u)),
    pct        = round(100 * c(nrow(s), nrow(u)) / nrow(df), 1),
    mean_age   = round(c(mean(s$age), mean(u$age)), 1),
    sd_age     = round(c(sd(s$age), sd(u$age)), 1),
    pct_female = round(c(100 * mean(s$sex == "F", na.rm = TRUE),
                         100 * mean(u$sex == "F", na.rm = TRUE)), 1),
    smd_age    = c(round(smd, 3), NA)
  )
}
