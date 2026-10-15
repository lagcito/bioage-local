# Step 3. Local reference: Method B recentering, Z score and percentile cut points.
for (f in c("config", "recenter")) source(file.path(PROJECT_ROOT, "R", paste0(f, ".R")))

sc  <- readRDS(file.path(PROJECT_ROOT, "data_derived", "local_scored.rds"))
df  <- sc$df
rec <- recenter_method_b(df, df$ba)
df$z    <- rec$z
df$band <- assign_bands(rec$z)

print(rec$summary, row.names = FALSE)
write.csv(cbind(rec$summary, as.data.frame(rec$coefs)),
          file.path(PROJECT_ROOT, "final_output", "recentering_summary.csv"), row.names = FALSE)
print(table(df$band))

# Patient-level log of the censored cases: review it locally, never commit it.
write.csv(df[rec$extreme, c("patient_id", "age", "sex", "ba", "flag_hematol")],
          file.path(PROJECT_ROOT, "data_derived", "censored_extreme_z.csv"), row.names = FALSE)
saveRDS(list(df = df, raw = sc$raw, rec = rec), file.path(PROJECT_ROOT, "data_derived", "local_recentered.rds"))
