# Criterion 1 — precision of graded response model (GRM) item parameters per
# domain, by N. Output: results/01_item_parameters.csv
suppressPackageStartupMessages(library(mirt))
source("R/00_assumptions.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(101)

a_true <- matrix(GRM_A)
d_true <- -GRM_A * GRM_B                    # mirt intercepts: d = -a * b

one_rep <- function(N) {
  dat <- simdata(a = a_true, d = d_true, N = N, itemtype = "graded")
  if (any(apply(dat, 2, function(x) length(unique(x))) < 4)) return(c(NA, NA, 1))
  fit <- tryCatch(mirt(dat, 1, "graded", verbose = FALSE, TOL = 1e-4),
                  error = function(e) NULL)
  if (is.null(fit) || !extract.mirt(fit, "converged")) return(c(NA, NA, 1))
  est <- coef(fit, simplify = TRUE, IRTpars = TRUE)$items
  if (mean(est[, "a"]) < 0) est[, "a"] <- -est[, "a"]   # arbitrary sign
  c(sqrt(mean((est[, "a"] - GRM_A)^2)),
    sqrt(mean((est[, c("b1", "b2", "b3")] - GRM_B)^2)), 0)
}

res <- do.call(rbind, lapply(GRM_N, function(N) {
  r <- do.call(rbind, parallel::mclapply(seq_len(GRM_REPS), function(i) one_rep(N),
                                         mc.cores = N_CORES))
  data.frame(N = N,
             rmse_discrimination = mean(r[, 1], na.rm = TRUE),
             rmse_thresholds     = mean(r[, 2], na.rm = TRUE),
             prop_rmse_discr_le_025 = mean(r[, 1] <= GRM_MAX_RMSE_A, na.rm = TRUE),
             prop_rmse_thr_le_025   = mean(r[, 2] <= GRM_MAX_RMSE_B, na.rm = TRUE),
             failures = mean(r[, 3]))
}))
write.csv(res, file.path(RESULTS, "01_item_parameters.csv"), row.names = FALSE)
print(res, digits = 3)
