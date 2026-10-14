###############################################################################
# pipeline/run_all.R — Orquestador único del pipeline completo
#
# Corre, en orden, los cinco pasos que antes había que ejecutar a mano en la
# misma sesión de R (run_train.R → run_score_argentina.R → run_recentering.R
# → run_validation.R → extract_results_for_paper.R). Cada paso sigue
# dependiendo del entorno global compartido con el anterior (bundle, cohort,
# results, rc, ...) — eso no cambió — pero ahora hay un único punto de
# entrada que fija PROJECT_ROOT una vez y llama a los pasos en el orden
# correcto, en vez de depender de que alguien los pegue a mano en orden.
#
# USO
#   Rscript pipeline/run_all.R
#   # o interactivamente:
#   source("pipeline/run_all.R")
#
# REQUIERE
#   data_raw/master_PARTE1.csv, data_raw/master_PARTE2.csv
#   data_raw/phenoage_master_PARTE1.csv, data_raw/phenoage_master_PARTE2.csv
#   (ver README.md — estos datos NO se versionan, están en .gitignore)
###############################################################################

PROJECT_ROOT <- tryCatch({
  normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))
}, error = function(e) normalizePath("."))
setwd(PROJECT_ROOT)
message(">> pipeline/run_all.R — PROJECT_ROOT = ", PROJECT_ROOT)

.required_raw <- file.path(PROJECT_ROOT, "data_raw",
  c("master_PARTE1.csv", "master_PARTE2.csv",
    "phenoage_master_PARTE1.csv", "phenoage_master_PARTE2.csv"))
.missing_raw <- .required_raw[!file.exists(.required_raw)]
if (length(.missing_raw) > 0) {
  stop("Faltan datos crudos en data_raw/ (no versionados, ver .gitignore):\n  ",
       paste(.missing_raw, collapse = "\n  "),
       "\nEste pipeline no puede correr sin ellos.")
}

steps <- c(
  "01_train.R",
  "02_score_argentina.R",
  "03_recenter.R",
  "04_validate.R",
  "05_extract_manuscript.R"
)

for (s in steps) {
  message("\n", strrep("=", 70))
  message(">> PASO: ", s)
  message(strrep("=", 70))
  source(file.path(PROJECT_ROOT, "pipeline", s), local = FALSE)
}

message("\n>> pipeline/run_all.R COMPLETO — ver final_output/ y data_derived/")
