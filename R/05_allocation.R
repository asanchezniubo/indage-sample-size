# Allocation of the sample, at a fixed size: (a) how the non-DS sample is spread
# across the eight age bands, and (b) how the DS minimum is spread across the
# three decades. Neither question changes the sample size; both change precision.
# Output: results/05_allocation_bands.csv and results/05_allocation_ds.csv
source("R/00_assumptions.R")
source("R/onset_functions.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(505)

# (a) Age-band profiles for the non-DS sample (n = N_NON), DS sample fixed.
bands <- do.call(rbind, lapply(names(ALLOC_PROFILES), function(nm) {
  non_b5 <- largest_remainder(sum(N_NON_BAND), ALLOC_PROFILES[[nm]])
  do.call(rbind, lapply(names(ONSET_SCENARIOS), function(sn) {
    sc <- ONSET_SCENARIOS[[sn]]
    t_ds  <- true_onset(sc$A0_ds,  sc$B_ds,  c(35, 64.9))
    t_non <- true_onset(sc$A0_non, sc$B_non, c(35, 79.9))
    r <- do.call(rbind, parallel::mclapply(seq_len(ONSET_REPS), function(i) {
      a_ds  <- sample_ages(c(largest_remainder(sum(N_DS_BAND), rep(DS_PER_DECADE/2, each = 2)), 0, 0))
      a_non <- sample_ages(non_b5)
      y_ds  <- mu(a_ds,  sc$A0_ds,  sc$B_ds)  + rnorm(length(a_ds))  + rnorm(length(a_ds),  0, sqrt(ERR_VAR))
      y_non <- mu(a_non, sc$A0_non, sc$B_non) + rnorm(length(a_non)) + rnorm(length(a_non), 0, sqrt(ERR_VAR))
      c(ds = onset_age(a_ds, y_ds, c(35, 64.9)), non = onset_age(a_non, y_non, c(35, 79.9)))
    }, mc.cores = N_CORES))
    summarise_onset(r, t_ds, t_non,
                    data.frame(profile = nm, scenario = sn,
                               band_65plus = sum(non_b5[7:8]),
                               pct_65plus = round(100 * sum(non_b5[7:8]) / sum(non_b5), 1)))
  }))
}))
write.csv(bands, file.path(RESULTS, "05_allocation_bands.csv"), row.names = FALSE)
print(bands[, c("profile", "scenario", "pct_65plus", "halfwidth_non_DS", "halfwidth_DS")], digits = 3)

# (b) Splits of the DS minimum across decades, at the adopted total.
ds <- do.call(rbind, lapply(names(DS_SPLITS), function(nm) {
  sp <- DS_SPLITS[[nm]]; stopifnot(sum(sp) == sum(N_DS_BAND))
  ds_b5 <- c(largest_remainder(sum(sp), rep(sp / 2, each = 2)), 0, 0)
  do.call(rbind, lapply(names(ONSET_SCENARIOS), function(sn) {
    sc <- ONSET_SCENARIOS[[sn]]
    t_ds  <- true_onset(sc$A0_ds,  sc$B_ds,  c(35, 64.9))
    t_non <- true_onset(sc$A0_non, sc$B_non, c(35, 79.9))
    r <- do.call(rbind, parallel::mclapply(seq_len(ONSET_REPS), function(i) {
      a_ds  <- sample_ages(ds_b5); a_non <- sample_ages(N_NON_BAND)
      y_ds  <- mu(a_ds,  sc$A0_ds,  sc$B_ds)  + rnorm(length(a_ds))  + rnorm(length(a_ds),  0, sqrt(ERR_VAR))
      y_non <- mu(a_non, sc$A0_non, sc$B_non) + rnorm(length(a_non)) + rnorm(length(a_non), 0, sqrt(ERR_VAR))
      c(ds = onset_age(a_ds, y_ds, c(35, 64.9)), non = onset_age(a_non, y_non, c(35, 79.9)))
    }, mc.cores = N_CORES))
    summarise_onset(r, t_ds, t_non, data.frame(ds_split = nm, scenario = sn))
  }))
}))
write.csv(ds, file.path(RESULTS, "05_allocation_ds.csv"), row.names = FALSE)
print(ds[, c("ds_split", "scenario", "halfwidth_DS", "undetected_DS", "power_difference")], digits = 3)
