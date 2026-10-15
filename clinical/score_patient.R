###############################################################################
# score_patient.R - score one patient (or a table of patients) with a registry entry
#
#   source("clinical/score_patient.R")
#   entry <- load_entry("final_output/registry_entry_draft.json")
#   score_patient(entry, age = 55, sex = "M", glucose = 95, creatinine = 0.9,
#                 lymph = 30, mcv = 88, rdw = 13.2, wbc = 6.5)
#   score_patients(entry, read.csv("clinical/example_patients_template.csv"))
#
# Inputs are in the units shown in the entry (glucose mg/dL, creatinine mg/dL,
# lymphocytes %, MCV fL, RDW %, WBC 10^3/uL). The maths is the same as
# docs/clock.js, which runs the calculator on the website.
#
# Research use only. The result is a position among local peers for a given age
# and sex. It does not predict mortality, measure the pace of aging or
# provide clinical advice.
###############################################################################

load_entry <- function(path = "final_output/registry_entry_draft.json") {
  jsonlite::read_json(path, simplifyVector = FALSE)
}

# Value on the model scale: value * scale, then ln() when the marker is logged.
to_model_scale <- function(marker, value) {
  v <- value * marker$scale
  if (isTRUE(marker$log)) log(v) else v
}

gompertz_biological_age <- function(model, x, age) {
  xb <- model$rate + sum(vapply(names(model$betas), function(k)
    model$betas[[k]] * (if (k == "age") age else x[[k]]), numeric(1)))
  M <- 1 - exp(model$m_n * exp(xb) / model$m_d)
  log(model$ba_n * log(1 - M)) / model$ba_d + model$ba_i
}

score_patient <- function(entry, age, sex, ..., quiet = FALSE) {
  values <- list(...)
  markers <- entry$clock$biomarkers
  x <- list()
  for (m in markers) {
    key <- m$key
    # accept the marker key or its common name (lncreat = creatinine, lymph = lymphocyte)
    alias <- c(lncreat = "creatinine", lymph = "lymphocyte")[key]
    v <- values[[key]] %||% (if (!is.na(alias)) values[[alias]])
    if (is.null(v)) stop("Missing value for ", m$label, " (", m$unit, ")")
    x[[key]] <- to_model_scale(m, v)
  }
  ba  <- gompertz_biological_age(entry$clock$model, x, age)
  rc  <- entry$recentering
  expected <- rc$intercept + rc$slope_age * age + rc$slope_male * (sex == "M")
  z   <- (ba - expected) / rc$resid_sd
  cut <- unlist(rc$z_cutpoints)
  pct <- approx(c(-5, cut, 5), c(0, 10, 25, 50, 75, 90, 100), xout = z, rule = 2)$y
  band <- as.character(cut(z, c(-Inf, cut["p10"], cut["p25"], cut["p75"], cut["p90"], Inf),
                           labels = c("P1", "P2", "P3", "P4", "P5")))
  out <- list(biological_age = ba, difference_to_age = ba - age, expected_biological_age = expected,
              z = z, percentile_approx = pct, band = band)
  if (!quiet) {
    cat(sprintf(paste0("Biological age %.1f y (%+.1f y vs chronological age)\n",
                       "Expected for this age and sex: %.1f y\n",
                       "Local Z score %+.2f | approximate percentile %.0f | group %s\n",
                       "Position among local peers; research use only, not clinical advice.\n"),
                ba, ba - age, expected, z, pct, band))
  }
  invisible(out)
}

score_patients <- function(entry, data) {
  keys <- vapply(entry$clock$biomarkers, function(m) m$key, "")
  rows <- lapply(seq_len(nrow(data)), function(i) {
    vals <- as.list(data[i, , drop = FALSE])
    names(vals) <- names(data)
    r <- do.call(score_patient, c(list(entry = entry, age = vals$age, sex = vals$sex, quiet = TRUE),
                                  vals[!names(vals) %in% c("age", "sex", "patient_id")]))
    as.data.frame(r)
  })
  cbind(data, do.call(rbind, rows))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
