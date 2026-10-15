# Step 4. External construct validity: is the score associated, in the expected
# direction, with biomarkers the clock never used?
for (f in c("config", "local_data", "recenter", "validate")) source(file.path(PROJECT_ROOT, "R", paste0(f, ".R")))

rc  <- readRDS(file.path(PROJECT_ROOT, "data_derived", "local_recentered.rds"))
df  <- rc$df; raw <- rc$raw

frame <- data.frame(z = df$z, sex = df$sex, flag_hematol = df$flag_hematol)
for (b in EXTERNAL_CANDIDATES) {
  if (b %in% c("nlr", "homa_ir")) frame[[b]] <- df[[b]]
  else if (b %in% names(raw))     frame[[b]] <- num(raw[[b]])
}

res <- external_validity(frame)
ext <- res$external
pure <- ext[!ext$derived_from_clock, ]
n_pass <- sum(meets_criterion(pure))
message(sprintf(">> %d of %d external markers meet |d| > %.2f and FDR < %.2f",
                n_pass, nrow(pure), EFFECT_THRESHOLD, FDR_THRESHOLD))
for (sx in c("F", "M")) {
  t <- res$by_sex[res$by_sex$sex == sx, ]
  message(sprintf("   %s: %d of %d", sx, sum(meets_criterion(t)), nrow(t)))
}

out <- function(x, f) write.csv(x, file.path(PROJECT_ROOT, "final_output", f), row.names = FALSE)
out(ext,               "external_validity.csv")
out(res$sensitivity,   "external_validity_without_hematology_flag.csv")
out(res$by_sex,        "external_validity_by_sex.csv")
out(res$ascertainment, "ascertainment.csv")
