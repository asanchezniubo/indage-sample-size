# Precision of the pilot decision rules (Section 3.2.4) by pilot size.
# Each NeuronUP test is completed by half of the sample (planned missing design).
# Proportions: Wald 95% interval; correlation: Fisher z; ICC: 95% interval width
# from Bonett (2002) with k = 2 measurements.
# Decision: pilot of 150 with moderate ID oversampled to half; retest n = 50.
# Output: results/06_pilot_precision.csv, results/06_pilot_retest.csv
source("R/00_assumptions.R")
ci_p <- function(p, n) 1.96 * sqrt(p * (1 - p) / n)
ci_r <- function(r, n) { z <- atanh(r); s <- 1 / sqrt(n - 3); tanh(c(z - 1.96 * s, z + 1.96 * s)) }
bonett_w <- function(rho, n, k = 2)
  sqrt(8 * 1.96^2 * (1 - rho)^2 * (1 + (k - 1) * rho)^2 / (k * (k - 1) * (n - 1)))
res <- do.call(rbind, lapply(c(100, 150, 200, 250), function(N) {
  nt <- N / 2; nm <- nt / 2                       # moderate ID = half of the sample
  data.frame(n_pilot = N, n_per_test = nt,
             halfwidth_floor15 = ci_p(.15, nt), halfwidth_completion80 = ci_p(.80, nt),
             r50_lower = ci_r(.5, nt)[1], r50_upper = ci_r(.5, nt)[2],
             n_moderate_per_test = nm, halfwidth_floor15_moderate = ci_p(.15, nm),
             halfwidth_icc60_informant_moderate = bonett_w(.6, N / 2) / 2)
}))
retest <- data.frame(n_retest = c(30, 40, 50, 60),
                     halfwidth_icc80 = bonett_w(.8, c(30, 40, 50, 60)) / 2)
write.csv(res, file.path(RESULTS, "06_pilot_precision.csv"), row.names = FALSE)
write.csv(retest, file.path(RESULTS, "06_pilot_retest.csv"), row.names = FALSE)
print(res, digits = 2); print(retest, digits = 2)
