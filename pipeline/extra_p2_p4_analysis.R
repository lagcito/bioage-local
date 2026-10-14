# Extra analysis: how do P2 and P4 (intermediate bands) behave for external biomarkers?
# Recomputes Method B Z (lm BA ~ age + sex, |Z|>5 censored), the same empirical
# bands as R/profile.R, then compares EVERY band against P3 (reference) per marker.
# Output: final_output/step4_bands_vs_P3.csv
suppressMessages(library(data.table))
setwd(if (file.exists("data_derived")) "." else "..")
lite <- fread("data_derived/clinical_scored_lite.csv")
raw  <- rbind(fread("data_raw/phenoage_master_PARTE1.csv", fill = TRUE),
              fread("data_raw/phenoage_master_PARTE2.csv", fill = TRUE), fill = TRUE)
raw[, patient_id := Protocolo]
raw <- raw[!duplicated(patient_id)]

ba <- "pheno_pheno_max_v2"
d  <- lite[!is.na(get(ba)) & is.finite(get(ba)) & !is.na(age) & sex %in% c("F", "M")]
fit <- lm(d[[ba]] ~ d$age + d$sex)
d[, z := as.numeric(scale(residuals(fit)))]
d[abs(z) > 5, z := NA_real_]
d <- d[!is.na(z)]
q <- quantile(d$z, c(.10, .25, .75, .90))
d[, band := cut(z, c(-Inf, q, Inf), labels = c("P1", "P2", "P3", "P4", "P5"))]
cat("valid Z:", nrow(d), "\n"); print(table(d$band))

bm_all <- setdiff(names(raw), c("Protocolo", "source_file", "age", "sex", "date_admission",
                                "date_result", "patient_id"))
dd <- merge(d[, .(patient_id, z, band)], raw, by = "patient_id", all.x = TRUE)
# nlr / homa_ir come from lite
dd <- merge(dd, lite[, .(patient_id, nlr_l = nlr, homa_ir_l = homa_ir)], by = "patient_id", all.x = TRUE)
dd[, nlr := nlr_l][, homa_ir := homa_ir_l]

clock <- c("rdw", "wbc", "mcv", "glucose", "creatinine", "lymphocyte")
derived <- c("nlr", "eGFR", "homa_ir")
ext <- read.csv("final_output/step4_external_construct_validity.csv")$biomarker
ext <- intersect(ext, names(dd))

cohen <- function(a, b) {  # mean(a) - mean(b) in pooled SD
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2))
  dd <- (mean(a) - mean(b)) / sp
  se <- sqrt((na + nb) / (na * nb) + dd^2 / (2 * (na + nb)))
  c(d = dd, lo = dd - 1.96 * se, hi = dd + 1.96 * se)
}
res <- list()
for (bm in ext) {
  x <- suppressWarnings(as.numeric(dd[[bm]])); ok <- is.finite(x) & !is.na(dd$band)
  x <- x[ok]; g <- dd$band[ok]
  if (sum(ok) < 50 || any(table(g) < 5)) next
  ref <- x[g == "P3"]
  row <- list(biomarker = bm, n = sum(ok), derived = bm %in% derived,
              rho = unname(cor(dd$z[ok], x, method = "spearman")))
  for (b in c("P1", "P2", "P4", "P5")) {
    xb <- x[g == b]; cd <- cohen(xb, ref)
    row[[paste0("n_", b)]] <- length(xb)
    row[[paste0("d_", b, "_vs_P3")]] <- cd["d"]
    row[[paste0("lo_", b)]] <- cd["lo"]; row[[paste0("hi_", b)]] <- cd["hi"]
    row[[paste0("p_", b)]] <- suppressWarnings(wilcox.test(xb, ref)$p.value)
  }
  m <- tapply(x, g, mean)
  row$monotone <- all(diff(m) > 0) || all(diff(m) < 0)
  res[[bm]] <- as.data.frame(row)
}
out <- do.call(rbind, res)
for (b in c("P1", "P2", "P4", "P5"))
  out[[paste0("fdr_", b)]] <- p.adjust(out[[paste0("p_", b)]], "BH")
write.csv(out, "final_output/step4_bands_vs_P3.csv", row.names = FALSE)
cat("markers:", nrow(out), "\n")
