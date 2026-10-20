###############################################################################
# config.R - everything a laboratory edits to adapt the pipeline to its data
###############################################################################

# -- Clock definition -------------------------------------------------------
CLOCK_NAME <- "six-bioage"

# Markers of the local clock. Names follow the BioAge package conventions
# (glucose in mmol/L, lncreat = ln of creatinine in umol/L). The survival model
# is trained in NHANES III, so only markers that exist there can be used.
PANEL_LOCAL <- c("glucose", "lncreat", "lymph", "mcv", "rdw", "wbc")

# Comparator: the nine markers of Levine et al. (2018)
PANEL_REFERENCE <- c("albumin", "alp", "glucose", "lncrp", "lncreat",
                     "lymph", "mcv", "rdw", "wbc")

# How each local marker is shown to a user and converted to the model scale:
# model value = value * scale, then ln() if log = TRUE.
BIOMARKER_META <- list(
  glucose = list(label = "Glucose",           unit = "mg/dL",   scale = 1 / 18, log = FALSE),
  lncreat = list(label = "Creatinine",        unit = "mg/dL",   scale = 88.4,   log = TRUE),
  lymph   = list(label = "Lymphocytes",       unit = "%",       scale = 1,      log = FALSE),
  mcv     = list(label = "MCV",               unit = "fL",      scale = 1,      log = FALSE),
  rdw     = list(label = "RDW",               unit = "%",       scale = 1,      log = FALSE),
  wbc     = list(label = "White blood cells", unit = "10^3/uL", scale = 1,      log = FALSE)
)

# -- Local export -----------------------------------------------------------
# Files (a glob relative to the repository root) and the identifier column.
# The expected columns and units are documented in the README.
LOCAL_FILES     <- "data_raw/phenoage_master_PARTE*.csv"
LOCAL_ID_COLUMN <- "Protocolo"
AGE_RANGE       <- c(18, 120)

# -- Recentering ------------------------------------------------------------
Z_CENSOR <- 5   # standardized residuals with |Z| above this are set to missing

# -- Hematology flag (marked, never excluded from scoring) -------------------
HEMATOL_WBC_MAX      <- 30   # 10^3/uL
HEMATOL_LYMPH_PCT    <- 70   # %
HEMATOL_LYMPH_ABSMAX <- 10   # 10^3/uL

# -- External validation ----------------------------------------------------
# Raw columns of the local export that are candidates for external validation.
EXTERNAL_CANDIDATES <- c(
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
# Raw columns that are inputs of the clock: excluded to avoid circularity.
CLOCK_INPUT_COLUMNS <- c("rdw", "wbc", "mcv", "glucose", "creatinine", "lymphocyte")
# Computed from clock inputs: reported but flagged as derived.
DERIVED_COLUMNS <- c("nlr", "eGFR", "homa_ir")
EFFECT_THRESHOLD <- 0.10   # |Cohen d| criterion
FDR_THRESHOLD    <- 0.05   # Benjamini-Hochberg criterion
