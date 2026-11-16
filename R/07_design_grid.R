# Choice of the sample size, with the two strata free: a grid over the size of
# the non-DS sample and the size of the DS sample. The adopted design is the
# smallest cell that meets both acceptance criteria under the conservative
# scenario (ONSET_MAX_HALFWIDTH and ONSET_MIN_POWER).
# Output: results/07_design_grid.csv
source("R/00_assumptions.R")
source("R/onset_functions.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(707)

res <- do.call(rbind, lapply(ONSET_N_NON, function(n_non)
  do.call(rbind, lapply(ONSET_N_DS, function(n_ds)
    do.call(rbind, lapply(names(ONSET_SCENARIOS), function(nm) {
      sc <- ONSET_SCENARIOS[[nm]]
      t_ds  <- true_onset(sc$A0_ds,  sc$B_ds,  c(35, 64.9))
      t_non <- true_onset(sc$A0_non, sc$B_non, c(35, 79.9))
      non_b5 <- largest_remainder(n_non, ALLOC_PROFILES[["adopted"]])
      ds_b5  <- c(largest_remainder(n_ds, rep(DS_PER_DECADE / 2, each = 2)), 0, 0)
      r <- do.call(rbind, parallel::mclapply(seq_len(ONSET_REPS), function(i) {
        a_ds <- sample_ages(ds_b5); a_non <- sample_ages(non_b5)
        y_ds  <- mu(a_ds,  sc$A0_ds,  sc$B_ds)  + rnorm(length(a_ds))  + rnorm(length(a_ds),  0, sqrt(ERR_VAR))
        y_non <- mu(a_non, sc$A0_non, sc$B_non) + rnorm(length(a_non)) + rnorm(length(a_non), 0, sqrt(ERR_VAR))
        c(ds = onset_age(a_ds, y_ds, c(35, 64.9)), non = onset_age(a_non, y_non, c(35, 79.9)))
      }, mc.cores = N_CORES))
      out <- summarise_onset(r, t_ds, t_non,
                             data.frame(scenario = nm, n_non_DS = n_non, n_DS = n_ds,
                                        total = n_non + n_ds))
      out$meets_criteria <- out$halfwidth_DS <= ONSET_MAX_HALFWIDTH &
                            out$power_difference >= ONSET_MIN_POWER
      out
    }))))))
write.csv(res, file.path(RESULTS, "07_design_grid.csv"), row.names = FALSE)
print(res[res$scenario == "conservative",
          c("n_non_DS", "n_DS", "total", "halfwidth_DS", "halfwidth_non_DS",
            "power_difference", "meets_criteria")], digits = 3)
