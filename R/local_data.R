###############################################################################
# local_data.R - read the local laboratory export and put it on the model scale
###############################################################################

read_local_data <- function() {
  files <- sort(Sys.glob(file.path(PROJECT_ROOT, LOCAL_FILES)))
  if (length(files) == 0) {
    stop("No local data found for the pattern '", LOCAL_FILES, "'. ",
         "Place your export in data_raw/ (never committed) and edit LOCAL_FILES in R/config.R.")
  }
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

num <- function(x) suppressWarnings(as.numeric(x))

# Columns every laboratory must provide; the rest are optional.
REQUIRED_COLUMNS <- c("age", "sex", "glucose", "creatinine", "lymphocyte", "mcv", "rdw", "wbc")

# Returns list(df = analysis frame on the model scale, raw = raw export, row-aligned).
prepare_local_data <- function(raw) {
  missing_cols <- setdiff(c(REQUIRED_COLUMNS, LOCAL_ID_COLUMN), names(raw))
  if (length(missing_cols) > 0) {
    stop("The local export is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  age  <- num(raw$age)
  keep <- !is.na(age) & age >= AGE_RANGE[1] & age <= AGE_RANGE[2]
  raw  <- raw[keep, , drop = FALSE]
  # optional columns become NA when the laboratory does not have them
  col  <- function(name) if (name %in% names(raw)) num(raw[[name]]) else rep(NA_real_, nrow(raw))

  glucose_mg <- col("glucose")                   # mg/dL
  crp        <- col("crp") / 10                  # mg/L -> mg/dL
  creat      <- col("creatinine") * 88.4         # mg/dL -> umol/L
  wbc        <- col("wbc")
  wbc        <- ifelse(wbc > 300, wbc / 1000, wbc)   # cells/uL -> 10^3/uL

  df <- data.frame(
    patient_id = raw[[LOCAL_ID_COLUMN]],
    age        = age[keep],
    sex        = raw$sex,
    albumin    = col("albumin") * 10,            # g/dL -> g/L
    alp        = col("alp"),
    glucose    = glucose_mg / 18,                # mg/dL -> mmol/L
    lncrp      = ifelse(is.na(crp)   | crp   <= 0, NA, log(crp)),
    lncreat    = ifelse(is.na(creat) | creat <= 0, NA, log(creat)),
    lymph      = col("lymphocyte"),
    mcv        = col("mcv"),
    rdw        = col("rdw"),
    wbc        = wbc,
    stringsAsFactors = FALSE
  )

  # Indices computed from raw columns
  neu <- col("neutrophils"); lym <- col("lymphocyte"); ins <- col("insulin")
  df$nlr     <- ifelse(!is.na(neu) & !is.na(lym) & lym > 0, neu / lym, NA)
  df$homa_ir <- ifelse(!is.na(glucose_mg) & !is.na(ins) & glucose_mg > 0 & ins > 0,
                       glucose_mg * ins / 405, NA)

  # Possible hematological disease: marked, not excluded
  lym_abs <- col("lymphocyte_abs")
  df$flag_hematol <- (!is.na(wbc) & wbc > HEMATOL_WBC_MAX) |
                     (!is.na(lym) & lym > HEMATOL_LYMPH_PCT) |
                     (!is.na(lym_abs) & lym_abs > HEMATOL_LYMPH_ABSMAX)

  list(df = df, raw = raw)
}

# Score every patient with a complete panel; others stay NA.
score_clock <- function(df, panel, fit) {
  cols <- c(panel, "age")
  ok   <- complete.cases(df[, cols, drop = FALSE])
  ba   <- rep(NA_real_, nrow(df))
  if (any(ok)) {
    calc <- BioAge::phenoage_calc(data = df[ok, cols, drop = FALSE], biomarkers = panel, fit = fit)
    ba[ok] <- calc$data$phenoage
  }
  ba
}
