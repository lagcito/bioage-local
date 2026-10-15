###############################################################################
# registry_entry.R - draft of the registry entry for this laboratory
#
# Writes model parameters and aggregate statistics only (no patient data).
# Review it, fill in the laboratory fields, validate it with
#   python scripts/build_registry.py --check
# and contribute it as registry/clocks/<id>.json through a pull request.
###############################################################################

build_registry_entry <- function(fit, rec, perf, df) {
  cf <- fit$coef
  betas <- as.list(setNames(cf[[1]], rownames(cf)))
  betas <- betas[!names(betas) %in% c("shape", "rate")]
  markers <- lapply(PANEL_LOCAL, function(k) {
    m <- BIOMARKER_META[[k]]
    list(key = k, label = m$label, unit = m$unit, scale = m$scale, log = m$log)
  })
  list(
    id = "your-lab-clock-name", example = FALSE, status = "self-reported",
    lab = list(name = "Your laboratory", country = "XX", type = "clinical laboratory"),
    clock = list(
      name = CLOCK_NAME, family = "PhenoAge-derived (Gompertz)",
      description = "Describe the markers, what the model was trained on and evaluated on.",
      trained_on = "NHANES III", evaluated_on = "NHANES IV",
      biomarkers = markers,
      model = list(type = "gompertz_phenoage", rate = cf["rate", 1], betas = betas,
                   m_n = fit$m_n, m_d = fit$m_d, ba_n = fit$BA_n, ba_d = fit$BA_d, ba_i = fit$BA_i)),
    recentering = list(
      method = "B: lm(BA ~ age + sex), standardized residual",
      reference_n = rec$coefs$n_fit, intercept = rec$coefs$intercept,
      slope_age = rec$coefs$slope_age, slope_male = rec$coefs$slope_male,
      resid_sd = rec$coefs$resid_sd, z_cutpoints = rec$cutpoints,
      cutpoints_n = rec$summary$n),
    performance = list(
      n_scored = perf$n, r_with_age = perf$r, rmse_regression_years = perf$RMSE,
      rmse_identity_years = perf$RMSE_identity, mean_ba_minus_ca_years = perf$mean_adv),
    cohort = list(n_eligible = nrow(df), pct_female = round(100 * mean(df$sex == "F", na.rm = TRUE), 1),
                  mean_age = round(mean(df$age), 1), sd_age = round(sd(df$age), 1)),
    license = "CC-BY-4.0", submitted = format(Sys.Date())
  )
}
