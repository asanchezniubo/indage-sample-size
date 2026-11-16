# Criterion 2 — stability of the factor solution of each part (logic of SENECA,
# Lorenzo-Seva & Ferrando, 2024): CRMSR between the population reproduced common
# variance (Lambda Phi Lambda') and the sample one (A A', rotation-invariant).
# Ordinal items (4 categories) -> polychoric correlations -> 5-factor EFA.
# Output: results/02_factor_structure.csv
suppressPackageStartupMessages(library(lavaan))
source("R/00_assumptions.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(202)

population_model <- function(spec) {
  m <- sum(spec$n_items); k <- length(spec$n_items)
  L <- matrix(0, m, k); row <- 0
  for (j in seq_len(k)) {
    idx <- row + seq_len(spec$n_items[j]); L[idx, j] <- spec$loadings[[j]]; row <- max(idx)
  }
  Phi <- matrix(spec$phi, k, k); diag(Phi) <- 1
  common <- L %*% Phi %*% t(L)
  Sigma <- common; diag(Sigma) <- 1
  list(L = L, Phi = Phi, common = common, Sigma = Sigma, m = m, k = k)
}

# Asymmetric thresholds (moderate ceiling): 15%, 25%, 35%, 25% per category.
TAU <- qnorm(c(.15, .40, .75))

crmsr_rep <- function(pop, N) {
  Z <- MASS::mvrnorm(N, rep(0, pop$m), pop$Sigma)
  X <- as.data.frame(apply(Z, 2, function(z) findInterval(z, TAU) + 1))
  if (any(sapply(X, function(x) length(unique(x))) < SEN_CATS)) return(NA)
  R <- tryCatch(lavCor(X, ordered = names(X)), error = function(e) NULL)
  if (is.null(R)) return(NA)
  R <- unclass(R)
  fa <- tryCatch(factanal(covmat = R, factors = pop$k, n.obs = N, rotation = "none"),
                 error = function(e) NULL)
  if (is.null(fa)) return(NA)
  A <- unclass(fa$loadings)
  D <- (A %*% t(A)) - pop$common
  sqrt(mean(D[lower.tri(D)]^2))
}

res <- do.call(rbind, lapply(c("IC", "FA"), function(part) {
  pop <- population_model(if (part == "IC") SEN_IC else SEN_FA)
  do.call(rbind, lapply(SEN_N, function(N) {
    v <- unlist(parallel::mclapply(seq_len(SEN_REPS), function(i) crmsr_rep(pop, N),
                                   mc.cores = N_CORES))
    data.frame(part = part, N = N, m = pop$m,
               mean_communality = round(mean(diag(pop$common)), 2),
               crmsr_median = median(v, na.rm = TRUE),
               prop_crmsr_le_05 = mean(v <= .05, na.rm = TRUE),
               prop_crmsr_le_04 = mean(v <= .04, na.rm = TRUE),
               prop_crmsr_le_03 = mean(v <= .03, na.rm = TRUE),
               failures = mean(is.na(v)))
  }))
}))
write.csv(res, file.path(RESULTS, "02_factor_structure.csv"), row.names = FALSE)
print(res, digits = 3)
