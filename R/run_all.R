# Runs every criterion in order (about 35 minutes on 7 cores, most of it in
# 03_invariance_dif) and writes the
# CSV files in results/. Run from the repository root: Rscript R/run_all.R
for (s in c("R/01_item_parameters.R", "R/02_factor_structure.R", "R/04_onset_age.R",
            "R/05_allocation.R", "R/03_invariance_dif.R", "R/06_pilot_precision.R", "R/07_design_grid.R")) {
  cat("\n==========", s, format(Sys.time()), "==========\n")
  t0 <- Sys.time(); source(s, echo = FALSE, local = new.env())
  cat("-- elapsed:", format(round(Sys.time() - t0, 1)), "\n")
}
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
