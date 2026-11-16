# Criterion 4 (binding) — precision of the onset age of decline in each stratum
# (DS / non-DS) and power to detect that it is earlier in DS, as a function of
# the size of the DS sample. The non-DS sample is fixed at N_NON_BAND, because
# the two strata are modelled separately and no longer share a fixed total.
# Onset age = youngest age at which the fitted mean falls ONSET_SD_DROP SDs
# below the mean of the 35-44 group of the same stratum.
# Output: results/04_onset_age.csv
source("R/00_assumptions.R")
source("R/onset_functions.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(404)

# A DS sample of size n keeps the adopted shape across decades (115/115/70).
ds_bands <- function(n) {
  w <- rep(DS_PER_DECADE / 2, each = 2)
  c(largest_remainder(n, w), 0, 0)
}

onset_rep <- function(sc, n_ds, non_b5 = N_NON_BAND) {
  a_ds  <- sample_ages(ds_bands(n_ds)); a_non <- sample_ages(non_b5)
  y_ds  <- mu(a_ds,  sc$A0_ds,  sc$B_ds)  + rnorm(length(a_ds))  + rnorm(length(a_ds),  0, sqrt(ERR_VAR))
  y_non <- mu(a_non, sc$A0_non, sc$B_non) + rnorm(length(a_non)) + rnorm(length(a_non), 0, sqrt(ERR_VAR))
  c(ds = onset_age(a_ds, y_ds, c(35, 64.9)), non = onset_age(a_non, y_non, c(35, 79.9)))
}

res <- do.call(rbind, lapply(names(ONSET_SCENARIOS), function(nm) {
  sc <- ONSET_SCENARIOS[[nm]]
  t_ds  <- true_onset(sc$A0_ds,  sc$B_ds,  c(35, 64.9))
  t_non <- true_onset(sc$A0_non, sc$B_non, c(35, 79.9))
  do.call(rbind, lapply(ONSET_N_DS, function(n_ds) {
    r <- do.call(rbind, parallel::mclapply(seq_len(ONSET_REPS),
                                           function(i) onset_rep(sc, n_ds),
                                           mc.cores = N_CORES))
    summarise_onset(r, t_ds, t_non,
                    data.frame(scenario = nm, n_DS = n_ds, n_non_DS = sum(N_NON_BAND),
                               total = n_ds + sum(N_NON_BAND)))
  }))
}))
write.csv(res, file.path(RESULTS, "04_onset_age.csv"), row.names = FALSE)
print(res, digits = 3)
