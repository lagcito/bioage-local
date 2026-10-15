###############################################################################
# train.R - train the local clock and a nine-marker comparator in NHANES III
#           with the BioAge package, and evaluate them in NHANES IV
###############################################################################

train_clock <- function(panel) {
  if (!requireNamespace("BioAge", quietly = TRUE)) {
    stop("BioAge is required: devtools::install_github('dayoonkwon/BioAge')")
  }
  BioAge::phenoage_nhanes(panel)
}

# Agreement of a biological age with chronological age.
age_metrics <- function(ba, age) {
  ok  <- is.finite(ba) & is.finite(age)
  ba  <- ba[ok]; age <- age[ok]
  fit <- lm(ba ~ age)
  s   <- summary(fit)
  data.frame(
    n             = length(ba),
    r             = round(cor(age, ba), 4),
    r_sq          = round(s$r.squared, 4),
    slope         = round(unname(coef(fit)[2]), 4),
    intercept     = round(unname(coef(fit)[1]), 2),
    RMSE          = round(sqrt(mean(s$residuals^2)), 2),   # around the BA ~ CA line
    RMSE_identity = round(sqrt(mean((ba - age)^2)), 2),    # around BA = CA
    mean_adv      = round(mean(ba - age), 2),
    sd_adv        = round(sd(ba - age), 2),
    row.names = NULL
  )
}

# NHANES IV performance of both clocks: complete samples, and the participants
# scored by both (the fair comparison).
evaluate_nhanes <- function(local, reference) {
  d_loc <- local$data[, c("sampleID", "age", "phenoage")]
  d_ref <- reference$data[, c("sampleID", "age", "phenoage")]
  full <- rbind(cbind(clock = CLOCK_NAME, sample = "complete cases",
                      age_metrics(d_loc$phenoage, d_loc$age)),
                cbind(clock = "levine_nine_markers", sample = "complete cases",
                      age_metrics(d_ref$phenoage, d_ref$age)))
  both <- merge(d_loc, d_ref, by = "sampleID", suffixes = c("_loc", "_ref"))
  both <- both[is.finite(both$phenoage_loc) & is.finite(both$phenoage_ref), ]
  shared <- rbind(cbind(clock = CLOCK_NAME, sample = "shared participants",
                        age_metrics(both$phenoage_loc, both$age_loc)),
                  cbind(clock = "levine_nine_markers", sample = "shared participants",
                        age_metrics(both$phenoage_ref, both$age_loc)))
  rbind(full, shared)
}
