###############################################################################
# score_single_patient_v3.R
# Calculadora de Edad Biológica — Referencia Poblacional Argentina
# PhenoAge (pheno_max_v2) | Versión 3
###############################################################################
#
# DESCRIPCIÓN
#   Calcula la edad biológica de un paciente con el reloj pheno_max_v2
#   (6 biomarcadores de Levine) y la compara contra la cohorte clínica
#   argentina (n = 71.421). Reporta:
#     · Edad biológica estimada y avance (BA − CA)
#     · Z-score argentino ajustado por edad y sexo (Método B)
#     · Percentil empírico exacto (ECDF de la cohorte argentina)
#     · Grupo de envejecimiento P1–P5 (cortes empíricos argentinos)
#     · Comparación aproximada vs referencia NHANES IV
#     · Flag de QC hematológico (leucemia / linfocitosis extrema)
#
# REQUIERE  (generados por pipeline/run_all.R)
#   data_derived/bioage_deployment_bundle.rds      — modelo pheno_max_v2 entrenado
#   data_derived/argentine_method_b_coefs.rds      — coeficientes Método B
#   data_derived/argentine_ecdf_reference.rds      — distribución de referencia ECDF
#
# ─────────────────────────────────────────────────────────────────────────
# USO RÁPIDO
#
#   source("clinical/score_patient.R")
#
#   # Un paciente
#   score_patient(age = 55, sex = "M",
#                 glucose = 95, creatinine = 0.9,
#                 lymph = 30, mcv = 88, rdw = 13.2, wbc = 6.5)
#
#   # Con ID y linfocitos absolutos (para QC)
#   score_patient(age = 62, sex = "F", patient_id = "PAC-001",
#                 glucose = 110, creatinine = 1.1,
#                 lymph = 25, mcv = 92, rdw = 14.5, wbc = 7.8,
#                 lymph_abs = 2.1)
#
#   # Lote desde data.frame
#   score_batch(mis_pacientes)
#
# UNIDADES DE ENTRADA
#   glucose     mg/dL       creatinine  mg/dL
#   lymph       %           mcv         fL
#   rdw         %           wbc         1000 cel/µL  o  cel/µL (autodetectado)
#   lymph_abs   1000 cel/µL (opcional, solo para QC hematológico)
###############################################################################

cat("══════════════════════════════════════════════════════════════\n")
cat("  CALCULADORA PHENOAGE v3 — Referencia Argentina + NHANES    \n")
cat("══════════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════════════════════
# 0. DIRECTORIO Y CONFIG
# ══════════════════════════════════════════════════════════════════════════

.SCRIPT_ROOT <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")),
  error = function(e) normalizePath("..")
)
.DATA_DERIVED <- file.path(.SCRIPT_ROOT, "data_derived")
.FINAL_OUTPUT <- file.path(.SCRIPT_ROOT, "final_output")
source(file.path(.SCRIPT_ROOT, "R", "config.R"), local = FALSE)

# ══════════════════════════════════════════════════════════════════════════
# 1. CARGAR MODELOS
# ══════════════════════════════════════════════════════════════════════════

# ── 1a. Modelo pheno_max_v2 (bundle de deployment) ──
.pmv2_fit <- NULL
.bundle_path <- file.path(.DATA_DERIVED, "bioage_deployment_bundle.rds")

if (file.exists(.bundle_path)) {
  .b <- readRDS(.bundle_path)
  if (!is.null(.b$subclock_models$pheno_max_v2)) {
    .pmv2_fit <- .b$subclock_models$pheno_max_v2
    message(">> Modelo pheno_max_v2 cargado")
  } else {
    message(">> ADVERTENCIA: pheno_max_v2 no encontrado en el bundle.")
  }
  rm(.b)
} else {
  message(">> ADVERTENCIA: bundle no encontrado — ", .bundle_path)
  message("   La puntuación de edad biológica no estará disponible.")
}

# ── 1b. Coeficientes Método B argentinos ──
.arg_coefs <- NULL
.arg_coefs_path <- file.path(.DATA_DERIVED, "argentine_method_b_coefs.rds")

if (file.exists(.arg_coefs_path)) {
  .arg_coefs <- readRDS(.arg_coefs_path)
  message(">> Coeficientes Método B cargados | relojes: ",
          paste(names(.arg_coefs), collapse = ", "))
} else {
  message(">> ADVERTENCIA: argentine_method_b_coefs.rds no encontrado.")
  message("   Ejecutar el pipeline completo para generarlo.")
}

# ── 1c. Distribución de referencia ECDF ──
.ecdf_fns <- list()
.ecdf_path <- file.path(.DATA_DERIVED, "argentine_ecdf_reference.rds")

if (file.exists(.ecdf_path)) {
  .ecdf_raw <- readRDS(.ecdf_path)
  for (.clk in names(.ecdf_raw)) {
    vals <- .ecdf_raw[[.clk]]
    if (length(vals) > 10) .ecdf_fns[[.clk]] <- ecdf(vals)
  }
  message(">> ECDF argentina cargada | n por reloj: ",
          paste(paste0(names(.ecdf_raw), "=",
                       format(sapply(.ecdf_raw, length), big.mark = ",")),
                collapse = ", "))
  rm(.ecdf_raw)
} else {
  message(">> INFO: argentine_ecdf_reference.rds no encontrado.")
  message("   Se usará interpolación de percentiles (menos exacta).")
}

# ══════════════════════════════════════════════════════════════════════════
# 2. CONSTANTES DE REFERENCIA
# ══════════════════════════════════════════════════════════════════════════

# Clave del reloj en los RDS de coeficientes
.CLOCK_KEY <- "pheno_pheno_max_v2"

# Cortes empíricos del Z-score argentino — se leen de
# final_output/step3_validation_summary.csv (regenerado por
# pipeline/03_recenter.R) para que esta calculadora nunca quede desalineada
# del pipeline. Los valores de abajo son sólo el fallback si ese CSV no
# está disponible (p. ej. al distribuir esta calculadora sola, sin el resto
# del repo) — último corrido conocido: cohorte n = 71.212, post-QC.
.ARG_CUTS <- c(p10 = -1.062, p25 = -0.637, p50 = -0.091,
               p75 =  0.455, p90 =  1.132)

.arg_cuts_path <- file.path(.FINAL_OUTPUT, "step3_validation_summary.csv")
if (file.exists(.arg_cuts_path)) {
  .vl <- tryCatch(read.csv(.arg_cuts_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(.vl) && nrow(.vl) >= 1 &&
      all(c("p10_empirical","p25_empirical","p75_empirical","p90_empirical") %in% names(.vl))) {
    .row <- .vl[1, ]
    .ARG_CUTS <- c(p10 = .row$p10_empirical, p25 = .row$p25_empirical,
                   p50 = if ("p50_empirical" %in% names(.vl)) .row$p50_empirical else .ARG_CUTS["p50"],
                   p75 = .row$p75_empirical, p90 = .row$p90_empirical)
    message(">> Cortes empíricos cargados desde step3_validation_summary.csv (n = ",
            format(.row$n, big.mark = ","), ")")
  } else {
    message(">> ADVERTENCIA: step3_validation_summary.csv sin las columnas esperadas — usando fallback hardcodeado.")
  }
  rm(.vl)
} else {
  message(">> ADVERTENCIA: no se encontró final_output/step3_validation_summary.csv — usando cortes de referencia hardcodeados (pueden estar desactualizados). Correr pipeline/03_recenter.R para regenerarlos.")
}

# Referencia NHANES IV para pheno_max_v2 (cohorte entrenamiento, n = 47.922)
.NHANES <- list(mean_adv = -1.65, sd_adv = 5.25, n = 47922,
                r = 0.9639, rmse = 5.20)

# Sesgo poblacional Argentina vs NHANES (mean_adv Argentina - mean_adv NHANES).
# Se recalcula desde step2_argentina_performance.csv si está disponible;
# -3.33 es el fallback del último corrido conocido.
.ARG_BIAS_VS_NHANES <- -3.33
.arg_perf_path <- file.path(.FINAL_OUTPUT, "step2_argentina_performance.csv")
if (file.exists(.arg_perf_path)) {
  .ap <- tryCatch(read.csv(.arg_perf_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(.ap) && "mean_adv" %in% names(.ap) && nrow(.ap) >= 1) {
    .ARG_BIAS_VS_NHANES <- .ap$mean_adv[1] - .NHANES$mean_adv
  }
  rm(.ap)
}

# ══════════════════════════════════════════════════════════════════════════
# 3. FUNCIONES AUXILIARES (internas, prefijo punto)
# ══════════════════════════════════════════════════════════════════════════

# ── 3a. Conversión de unidades clínicas → SI ──
.to_si <- function(glucose_mgdL, creatinine_mgdL, wbc_val) {
  list(
    glucose = glucose_mgdL / 18.0,            # mg/dL  → mmol/L
    lncreat = log(creatinine_mgdL * 88.4),    # mg/dL  → ln(µmol/L)
    wbc     = if (wbc_val > 300) wbc_val / 1000 else wbc_val  # → 10³/µL
  )
}

# ── 3b. Puntuación Gompertz con el fit object ──
.gompertz_score <- function(fit, bm, age) {
  if (is.null(fit) || is.null(fit$coef)) return(NA_real_)
  coefs <- if (is.data.frame(fit$coef)) {
    setNames(fit$coef[, 1], rownames(fit$coef))
  } else fit$coef
  rate  <- if ("rate" %in% names(coefs)) coefs["rate"] else 0
  betas <- coefs[!names(coefs) %in% c("shape", "rate")]
  x     <- c(bm, age = age)
  miss  <- setdiff(names(betas), names(x))
  if (length(miss) > 0) {
    warning("Biomarcadores faltantes: ", paste(miss, collapse = ", "))
    return(NA_real_)
  }
  xb  <- sum(betas[names(betas)] * x[names(betas)]) + rate
  M   <- 1 - exp(fit$m_n * exp(xb) / fit$m_d)
  as.numeric(log(fit$BA_n * log(1 - M)) / fit$BA_d + fit$BA_i)
}

# ── 3c. Percentil empírico (ECDF o interpolación por piezas) ──
.pct_empirical <- function(resid, clock_key = .CLOCK_KEY) {
  if (is.na(resid)) return(NA_real_)
  if (!is.null(.ecdf_fns[[clock_key]])) {
    return(round(.ecdf_fns[[clock_key]](resid) * 100, 1))
  }
  # Fallback: interpolación lineal por tramos desde cuantiles conocidos
  resid_sd <- if (!is.null(.arg_coefs[[clock_key]]))
    as.numeric(.arg_coefs[[clock_key]]$resid_sd) else 5.3551
  z <- resid / resid_sd
  z_pts  <- c(-5, .ARG_CUTS["p10"], .ARG_CUTS["p25"], .ARG_CUTS["p50"],
               .ARG_CUTS["p75"], .ARG_CUTS["p90"], 5)
  p_pts  <- c(0.5, 10, 25, 50, 75, 90, 99.5)
  round(approx(z_pts, p_pts, xout = z, rule = 2)$y, 1)
}

# ── 3d. Grupo de envejecimiento P1–P5 ──
.ager_group <- function(z) {
  if (is.na(z)) return(list(id = NA, label = "No calculable"))
  if (z <  .ARG_CUTS["p10"]) return(list(id = "P1", label = "Envejecimiento muy lento  (< P10)"))
  if (z <  .ARG_CUTS["p25"]) return(list(id = "P2", label = "Envejecimiento lento       (P10-P25)"))
  if (z <= .ARG_CUTS["p75"]) return(list(id = "P3", label = "Envejecimiento normal      (P25-P75)"))
  if (z <= .ARG_CUTS["p90"]) return(list(id = "P4", label = "Envejecimiento acelerado   (P75-P90)"))
  return(list(id = "P5",     label = "Envejecimiento muy acelerado (> P90)"))
}

# ── 3e. Flag QC hematológico ──
.hematol_flag <- function(wbc_k, lymph_pct, lymph_abs = NA) {
  flags <- character(0)
  if (!is.na(wbc_k)     && wbc_k     > 30) flags <- c(flags, sprintf("WBC = %.1f ×10³/µL (umbral >30)", wbc_k))
  if (!is.na(lymph_pct) && lymph_pct > 70) flags <- c(flags, sprintf("Linfocitos = %.1f%% (umbral >70%%)", lymph_pct))
  if (!is.na(lymph_abs) && lymph_abs > 10) flags <- c(flags, sprintf("Linfocitos abs = %.1f ×10³/µL (umbral >10)", lymph_abs))
  list(flagged = length(flags) > 0, details = flags)
}

# ── 3f. Texto interpretativo según grupo ──
.interpret <- function(group_id, pct, adv) {
  adv_str <- if (!is.na(adv)) {
    if (adv >= 0) sprintf(" (%+.1f años)", adv) else sprintf(" (%.1f años)", adv)
  } else ""
  switch(group_id,
    P1 = paste0("Envejecimiento MUY LENTO. Este paciente está en el percentil ", pct,
                " de la cohorte argentina", adv_str, ". ",
                "Su perfil de biomarcadores es consistente con un ritmo de envejecimiento ",
                "significativamente más lento que el promedio. Bajo riesgo biológico."),
    P2 = paste0("Envejecimiento LENTO. Percentil ", pct, adv_str, ". ",
                "Los biomarcadores sugieren un envejecimiento más lento que la media. ",
                "Buen perfil de riesgo biológico."),
    P3 = paste0("Envejecimiento NORMAL. Percentil ", pct, adv_str, ". ",
                "El ritmo de envejecimiento es comparable al esperado para edad y sexo. ",
                "Perfil de riesgo dentro del rango normal."),
    P4 = paste0("Envejecimiento ACELERADO. Percentil ", pct, adv_str, ". ",
                "Los biomarcadores sugieren un envejecimiento más rápido que el promedio. ",
                "Se recomienda revisar factores de riesgo modificables (glucemia, inflamación, función renal)."),
    P5 = paste0("Envejecimiento MUY ACELERADO. Percentil ", pct, adv_str, ". ",
                "El perfil de biomarcadores indica un envejecimiento acelerado en forma marcada. ",
                "Se recomienda evaluación clínica integral y manejo activo de factores de riesgo."),
    "No disponible."
  )
}

# ── 3g. Informe impreso ──
.print_v3 <- function(r) {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║       INFORME DE EDAD BIOLÓGICA — PhenoAge Argentina         ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  cat(sprintf("  Paciente : %-20s   Sexo : %s\n",
              r$patient_id,
              if (is.na(r$sex) || r$sex == "—") "no informado" else r$sex))
  cat(sprintf("  Reloj    : pheno_max_v2 (6/9 biomarcadores Levine)\n\n"))

  cat("  ── Biomarcadores ingresados ─────────────────────────────────\n")
  cat(sprintf("    Glucosa       : %6.1f mg/dL\n", r$glucose_mgdL))
  cat(sprintf("    Creatinina    : %6.2f mg/dL\n", r$creatinine_mgdL))
  cat(sprintf("    Linfocitos    : %6.1f %%\n",    r$lymph_pct))
  cat(sprintf("    VCM           : %6.1f fL\n",    r$mcv_fL))
  cat(sprintf("    ADE (RDW)     : %6.1f %%\n",    r$rdw_pct))
  cat(sprintf("    Leucocitos    : %6.2f ×10³/µL\n", r$wbc_k))

  # QC flag
  if (r$flag_hematol) {
    cat("\n")
    cat("  ⚠  ALERTA DE CALIDAD HEMATOLÓGICA\n")
    for (d in r$flag_details) cat("     →", d, "\n")
    cat("     Valores extremos pueden corresponder a neoplasia hematológica.\n")
    cat("     Interpretar el Z-score con precaución.\n")
  }

  cat("\n")
  cat("  ── Resultado principal ──────────────────────────────────────\n")
  if (is.na(r$phenoage)) {
    cat("    ** Edad biológica no disponible (bundle no cargado) **\n")
  } else {
    cat(sprintf("    Edad cronológica  : %5.1f años\n", r$age))
    cat(sprintf("    Edad biológica    : %5.1f años\n", r$phenoage))
    cat(sprintf("    Avance biológico  : %+5.1f años  ",  r$advancement))
    if      (r$advancement >  5) cat("(biológicamente más VIEJO)\n")
    else if (r$advancement < -5) cat("(biológicamente más JOVEN)\n")
    else                          cat("(dentro del rango esperado)\n")
  }

  cat("\n")
  cat("  ── Posición en la población argentina ──────────────────────\n")
  cat(sprintf("    Z-score argentino (aj. edad+sexo) : %+.2f SD\n",
              ifelse(is.na(r$z_argentina), NA, r$z_argentina)))
  cat(sprintf("    Percentil empírico                : %.1f°\n",
              ifelse(is.na(r$pct_empirical), NA, r$pct_empirical)))
  cat(sprintf("    Grupo de envejecimiento           : %s — %s\n",
              ifelse(is.na(r$ager_group), "N/A", r$ager_group),
              r$ager_label))

  if (!is.null(r$sex_note)) {
    cat(sprintf("\n    NOTA: %s\n", r$sex_note))
  }

  cat("\n")
  cat("  ── Comparación vs NHANES IV (referencia internacional) ─────\n")
  cat(sprintf("    Avance medio NHANES IV   : %+.2f ± %.2f años  (n = %s)\n",
              .NHANES$mean_adv, .NHANES$sd_adv,
              format(.NHANES$n, big.mark = ".")))
  if (!is.na(r$advancement)) {
    diff_vs_nhanes <- r$advancement - .NHANES$mean_adv
    cat(sprintf("    Avance del paciente      : %+.2f años\n", r$advancement))
    cat(sprintf("    Diferencia vs NHANES     : %+.2f años\n", diff_vs_nhanes))
    cat(sprintf("    Z NHANES (aprox.)        : %+.2f SD\n",   r$z_nhanes_approx))
    cat(sprintf("\n    Sesgo poblacional argentino vs NHANES: %.2f años\n",
                .ARG_BIAS_VS_NHANES))
    cat("    (El reloj subestima la edad en la población argentina\n")
    cat("     respecto a la referencia NHANES; el Z argentino corrige este sesgo.)\n")
  }

  cat("\n")
  cat("  ── Interpretación clínica ───────────────────────────────────\n")
  if (!is.na(r$ager_group)) {
    wrapped <- strwrap(r$interpretation, width = 62, prefix = "    ", initial = "    ")
    cat(paste(wrapped, collapse = "\n"), "\n")
  } else {
    cat("    No disponible.\n")
  }

  cat("\n")
  cat("  ── Referencia de grupos (cortes empíricos, n = 71.421) ─────\n")
  cat("    P1  Z < −1.062  ( < P10 )   Envejecimiento muy lento\n")
  cat("    P2  −1.062 a −0.637 (P10-25) Envejecimiento lento\n")
  cat("    P3  −0.637 a +0.455 (P25-75) Envejecimiento normal\n")
  cat("    P4  +0.455 a +1.132 (P75-90) Envejecimiento acelerado\n")
  cat("    P5  Z >  +1.132  ( > P90 )  Envejecimiento muy acelerado\n")
  cat("────────────────────────────────────────────────────────────────\n\n")
}

# ══════════════════════════════════════════════════════════════════════════
# 4. FUNCIÓN PRINCIPAL: score_patient()
# ══════════════════════════════════════════════════════════════════════════

#' Calcular y reportar la edad biológica de un paciente
#'
#' @param age          Edad cronológica (años)
#' @param sex          "M" o "F" [recomendado para Z-score preciso]
#' @param glucose      Glucosa sérica (mg/dL)
#' @param creatinine   Creatinina sérica (mg/dL)
#' @param lymph        Linfocitos (%)
#' @param mcv          Volumen corpuscular medio (fL)
#' @param rdw          ADE / RDW (%)
#' @param wbc          Leucocitos totales (1000 cel/µL o cel/µL, autodetectado)
#' @param lymph_abs    Linfocitos absolutos (1000 cel/µL) [opcional, solo QC]
#' @param patient_id   Identificador [opcional]
#' @param print_report Imprimir informe en consola (defecto: TRUE)
#' @return Lista con todos los resultados (invisible)
score_patient <- function(age,
                          sex          = NULL,
                          glucose,
                          creatinine,
                          lymph,
                          mcv,
                          rdw,
                          wbc,
                          lymph_abs    = NA,
                          patient_id   = NULL,
                          print_report = TRUE) {

  # ── Validaciones de entrada ──
  if (is.na(age) || age < 10 || age > 120)
    stop("Edad inválida: ", age, ". Rango esperado: 10–120 años.")
  required <- c(glucose = glucose, creatinine = creatinine,
                lymph = lymph, mcv = mcv, rdw = rdw, wbc = wbc)
  if (any(is.na(required)))
    stop("Biomarcador(es) faltante(s): ",
         paste(names(required)[is.na(required)], collapse = ", "))
  if (!is.null(sex) && !sex %in% c("M", "F"))
    stop("Sexo inválido: '", sex, "'. Usar 'M' o 'F'.")

  # ── Conversión de unidades ──
  si   <- .to_si(glucose, creatinine, wbc)
  wbc_k <- si$wbc

  # ── QC hematológico ──
  qc <- .hematol_flag(wbc_k, lymph, lymph_abs)

  # ── Edad biológica ──
  bm <- c(glucose = si$glucose, lncreat = si$lncreat,
          lymph = lymph, mcv = mcv, rdw = rdw, wbc = wbc_k)
  phenoage    <- .gompertz_score(.pmv2_fit, bm, age)
  advancement <- if (!is.na(phenoage)) phenoage - age else NA_real_

  # ── Z-score argentino (Método B: lm(BA ~ age + sex)) ──
  z_arg <- NA_real_; resid_arg <- NA_real_
  sex_note <- NULL

  if (!is.null(.arg_coefs) && .CLOCK_KEY %in% names(.arg_coefs)) {
    cf <- .arg_coefs[[.CLOCK_KEY]]

    # Codificación de sexo: F = 0, M = 1
    # Sin sexo: promedio ponderado (35% M / 65% F según cohorte argentina)
    sex_num <- if (!is.null(sex)) {
      if (sex == "M") 1 else 0
    } else {
      sex_note <- "Sexo no informado. Z-score calculado con proporción media (35% M / 65% F)."
      0.35
    }

    predicted_BA <- as.numeric(cf$intercept) +
                    as.numeric(cf$slope_age) * age +
                    as.numeric(cf$slope_sexM) * sex_num
    resid_arg <- phenoage - predicted_BA
    z_arg     <- resid_arg / as.numeric(cf$resid_sd)
  }

  # ── Percentil empírico ──
  pct <- .pct_empirical(resid_arg)

  # ── Grupo de envejecimiento ──
  ag <- .ager_group(z_arg)

  # ── Z-score NHANES (aproximado, sin ajuste de edad) ──
  z_nhanes_approx <- if (!is.na(advancement)) {
    (advancement - .NHANES$mean_adv) / .NHANES$sd_adv
  } else NA_real_

  # ── Texto de interpretación ──
  interp <- if (!is.na(ag$id)) {
    .interpret(ag$id, pct, advancement)
  } else "Edad biológica no calculable con los datos proporcionados."

  # ── Resultado ──
  result <- list(
    patient_id      = if (!is.null(patient_id)) as.character(patient_id) else "—",
    age             = age,
    sex             = if (!is.null(sex)) sex else "—",
    # Biomarcadores (unidades clínicas)
    glucose_mgdL    = glucose,
    creatinine_mgdL = creatinine,
    lymph_pct       = lymph,
    mcv_fL          = mcv,
    rdw_pct         = rdw,
    wbc_k           = wbc_k,
    # QC
    flag_hematol    = qc$flagged,
    flag_details    = qc$details,
    # Edad biológica
    phenoage        = if (!is.na(phenoage))    round(phenoage, 2)    else NA_real_,
    advancement     = if (!is.na(advancement)) round(advancement, 2) else NA_real_,
    # Referencia argentina
    residual_arg    = if (!is.na(resid_arg)) round(resid_arg, 2) else NA_real_,
    z_argentina     = if (!is.na(z_arg))     round(z_arg, 2)     else NA_real_,
    pct_empirical   = pct,
    ager_group      = ag$id,
    ager_label      = ag$label,
    # Referencia NHANES
    z_nhanes_approx = if (!is.na(z_nhanes_approx)) round(z_nhanes_approx, 2) else NA_real_,
    # Interpretación
    interpretation  = interp,
    sex_note        = sex_note
  )

  if (print_report) .print_v3(result)
  invisible(result)
}

# ══════════════════════════════════════════════════════════════════════════
# 5. FUNCIÓN BATCH: score_batch()
# ══════════════════════════════════════════════════════════════════════════

#' Puntuar un lote de pacientes desde un data.frame
#'
#' @param patients data.frame con columnas: age, glucose, creatinine, lymph,
#'                 mcv, rdw, wbc. Opcionales: sex, patient_id, lymph_abs.
#' @param verbose  TRUE para imprimir informe individual de cada paciente.
#' @return data.frame con resultados de todos los pacientes
score_batch <- function(patients, verbose = FALSE) {

  req_cols <- c("age", "glucose", "creatinine", "lymph", "mcv", "rdw", "wbc")
  miss_col <- setdiff(req_cols, names(patients))
  if (length(miss_col) > 0)
    stop("Columnas faltantes en el data.frame: ", paste(miss_col, collapse = ", "))

  if (!"patient_id" %in% names(patients))
    patients$patient_id <- paste0("PAC-", sprintf("%03d", seq_len(nrow(patients))))

  rows <- vector("list", nrow(patients))

  for (i in seq_len(nrow(patients))) {
    p <- patients[i, ]
    res <- tryCatch(
      score_patient(
        age          = p$age,
        sex          = if ("sex"       %in% names(p)) p$sex       else NULL,
        glucose      = p$glucose,
        creatinine   = p$creatinine,
        lymph        = p$lymph,
        mcv          = p$mcv,
        rdw          = p$rdw,
        wbc          = p$wbc,
        lymph_abs    = if ("lymph_abs" %in% names(p)) p$lymph_abs else NA,
        patient_id   = p$patient_id,
        print_report = verbose
      ),
      error = function(e) {
        message("  [ERROR] Paciente ", p$patient_id, ": ", e$message)
        list(patient_id = p$patient_id, age = p$age, sex = NA,
             phenoage = NA, advancement = NA,
             z_argentina = NA, z_nhanes_approx = NA,
             pct_empirical = NA, ager_group = NA, ager_label = "ERROR",
             flag_hematol = NA, interpretation = paste("ERROR:", e$message))
      }
    )
    rows[[i]] <- data.frame(
      patient_id      = res$patient_id,
      age             = res$age,
      sex             = if (!is.null(res$sex)) res$sex else NA_character_,
      phenoage        = res$phenoage,
      advancement     = res$advancement,
      z_argentina     = res$z_argentina,
      pct_empirical   = res$pct_empirical,
      ager_group      = res$ager_group,
      ager_label      = res$ager_label,
      z_nhanes_approx = res$z_nhanes_approx,
      flag_hematol    = res$flag_hematol,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  cat("\n")
  cat("══════════════════════════════════════════════════════════════\n")
  cat(sprintf("  BATCH COMPLETADO: %d pacientes procesados\n", nrow(out)))
  cat("══════════════════════════════════════════════════════════════\n\n")

  ok  <- sum(!is.na(out$z_argentina))
  err <- sum(is.na(out$z_argentina))
  cat(sprintf("  Puntuados con éxito : %d\n", ok))
  cat(sprintf("  Errores / incompletos: %d\n\n", err))

  if (ok > 0) {
    cat("  Distribución de grupos de envejecimiento:\n")
    grp_tab <- table(out$ager_group[!is.na(out$ager_group)])
    for (g in c("P1", "P2", "P3", "P4", "P5")) {
      n_g <- if (g %in% names(grp_tab)) grp_tab[g] else 0
      bar <- paste(rep("█", round(n_g / ok * 20)), collapse = "")
      cat(sprintf("    %s : %3d  %s\n", g, n_g, bar))
    }
    n_flag <- sum(out$flag_hematol, na.rm = TRUE)
    if (n_flag > 0)
      cat(sprintf("\n  ⚠  %d paciente(s) con flag hematológico\n", n_flag))
    cat("\n")
  }

  print(out[, c("patient_id", "age", "sex", "phenoage", "advancement",
                "z_argentina", "pct_empirical", "ager_group",
                "z_nhanes_approx", "flag_hematol")],
        row.names = FALSE)
  cat("\n")
  invisible(out)
}

# ══════════════════════════════════════════════════════════════════════════
# 6. CARGA DESDE ARCHIVO CSV / EXCEL
# ══════════════════════════════════════════════════════════════════════════

#' Cargar pacientes desde CSV y calcular edad biológica
#'
#' El CSV debe tener las columnas (con exactamente estos nombres):
#'   patient_id, age, sex, glucose, creatinine, lymph, mcv, rdw, wbc
#' Columna opcional: lymph_abs
#'
#' La plantilla lista para completar está en: plantilla_pacientes.csv
#'
#' @param archivo   Ruta al archivo CSV (defecto: "plantilla_pacientes.csv")
#' @param guardar   TRUE para exportar resultados a CSV (defecto: TRUE)
#' @param verbose   TRUE para imprimir informe individual de cada paciente
#' @return data.frame con resultados (también guardado en CSV si guardar=TRUE)
cargar_y_puntuar <- function(archivo = "plantilla_pacientes.csv",
                              guardar = TRUE,
                              verbose = FALSE) {

  # ── Leer archivo ──
  if (!file.exists(archivo))
    stop("Archivo no encontrado: ", archivo,
         "\nCompletá la plantilla en: plantilla_pacientes.csv")

  ext <- tolower(tools::file_ext(archivo))

  if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Instalar el paquete readxl para leer Excel:\n  install.packages('readxl')")
    pacientes <- as.data.frame(readxl::read_excel(archivo))
    message(">> Excel cargado: ", nrow(pacientes), " filas")
  } else {
    pacientes <- read.csv(archivo, stringsAsFactors = FALSE,
                          na.strings = c("", "NA", "N/A", "-"))
    message(">> CSV cargado: ", nrow(pacientes), " filas")
  }

  # Eliminar filas completamente vacías
  pacientes <- pacientes[rowSums(!is.na(pacientes)) > 1, ]

  if (nrow(pacientes) == 0)
    stop("El archivo está vacío o solo tiene encabezado.")

  cat(sprintf("   %d paciente(s) a procesar\n\n", nrow(pacientes)))

  # ── Puntuar ──
  resultados <- score_batch(pacientes, verbose = verbose)

  # ── Guardar resultados ──
  if (guardar) {
    base_name  <- tools::file_path_sans_ext(archivo)
    out_path   <- paste0(base_name, "_resultados.csv")
    write.csv(resultados, out_path, row.names = FALSE)
    cat(sprintf("\n>> Resultados guardados en: %s\n\n", out_path))
  }

  invisible(resultados)
}

#' Crear (o mostrar) la plantilla CSV lista para completar en Excel
plantilla <- function() {
  path <- file.path(.SCRIPT_ROOT, "plantilla_pacientes.csv")
  if (!file.exists(path)) {
    header <- paste(
      "patient_id,age,sex,glucose,creatinine,lymph,mcv,rdw,wbc,lymph_abs",
      "PAC-001,55,M,95,0.9,30,88,13.2,6.5,",
      "PAC-002,62,F,110,1.1,25,92,14.5,7.8,",
      sep = "\n"
    )
    writeLines(header, path)
  }
  cat("\n")
  cat("  Plantilla CSV disponible en:\n")
  cat("  ", normalizePath(path), "\n\n")
  cat("  Columnas requeridas:\n")
  cat("    patient_id  — identificador (texto libre)\n")
  cat("    age         — edad en años\n")
  cat("    sex         — M o F\n")
  cat("    glucose     — glucosa (mg/dL)\n")
  cat("    creatinine  — creatinina (mg/dL)\n")
  cat("    lymph       — linfocitos (%)\n")
  cat("    mcv         — VCM (fL)\n")
  cat("    rdw         — ADE / RDW (%)\n")
  cat("    wbc         — leucocitos (1000 cel/µL o cel/µL)\n")
  cat("  Columna opcional:\n")
  cat("    lymph_abs   — linfocitos absolutos (1000 cel/µL) [solo para QC]\n\n")
  cat("  Flujo de trabajo:\n")
  cat("    1. Abrí plantilla_pacientes.csv en Excel\n")
  cat("    2. Completá tus pacientes (podés agregar filas)\n")
  cat("    3. Guardá como CSV (mismo nombre o uno nuevo)\n")
  cat("    4. En R: cargar_y_puntuar('plantilla_pacientes.csv')\n\n")
  invisible(path)
}

# ══════════════════════════════════════════════════════════════════════════
# LISTO
# ══════════════════════════════════════════════════════════════════════════

cat(">> Calculadora lista. Funciones disponibles:\n\n")
cat("   score_patient()       — 1 paciente (parámetros directos)\n")
cat("   score_batch()         — lote desde data.frame\n")
cat("   cargar_y_puntuar()    — leer CSV/Excel y puntuar\n")
cat("   plantilla()           — ver instrucciones de la plantilla CSV\n\n")
cat("   Flujo más simple:\n")
cat("     1. plantilla()                            # ver instrucciones\n")
cat("     2. [completar plantilla_pacientes.csv]\n")
cat("     3. cargar_y_puntuar()                     # leer, calcular y guardar\n\n")
cat("   Ejemplo paciente individual:\n")
cat('   score_patient(age=55, sex="M",\n')
cat('                 glucose=95, creatinine=0.9,\n')
cat('                 lymph=30, mcv=88, rdw=13.2, wbc=6.5)\n\n')

if (is.null(.pmv2_fit))
  cat("   ** Modelo no cargado — falta bioage_deployment_bundle.rds\n\n")
if (is.null(.arg_coefs))
  cat("   ** Coeficientes argentinos no cargados — ejecutar el pipeline primero\n\n")
if (length(.ecdf_fns) == 0)
  cat("   ** ECDF no cargada — percentiles aproximados (ejecutar pipeline para exactos)\n\n")
