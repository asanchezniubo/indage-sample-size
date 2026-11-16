# Criterion 3 — power to detect lack of invariance (uniform DIF) in one item,
# with the other 4 items of the domain as anchors. The focal group also has a
# latent mean 0.5 SD lower (impact), as expected at older ages and in DS.
# Output: results/03_invariance_dif.csv
suppressPackageStartupMessages(library(mirt))
source("R/00_assumptions.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(303)

a <- matrix(GRM_A); b <- GRM_B
d_of <- function(a, b) -as.vector(a) * b

dif_rep <- function(n_ref, n_foc, shift) {
  th_r <- matrix(rnorm(n_ref)); th_f <- matrix(rnorm(n_foc, -0.5))
  b_f <- b; b_f[5, ] <- b_f[5, ] + shift
  d1 <- simdata(a = a, d = d_of(a, b), Theta = th_r, itemtype = "graded")
  d2 <- simdata(a = a, d = d_of(a, b_f), Theta = th_f, itemtype = "graded")
  dat <- rbind(d1, d2); g <- rep(c("ref", "foc"), c(n_ref, n_foc))
  if (any(sapply(split(as.data.frame(dat), g), function(D)
    any(sapply(D, function(x) length(unique(x))) < 4)))) return(NA)
  it <- colnames(dat)
  m0 <- tryCatch(multipleGroup(dat, 1, g, verbose = FALSE, quadpts = 31, TOL = 1e-3,
                   invariance = c(it, "free_means", "free_var")), error = function(e) NULL)
  m1 <- tryCatch(multipleGroup(dat, 1, g, verbose = FALSE, quadpts = 31, TOL = 1e-3,
                   invariance = c(it[1:4], "free_means", "free_var")), error = function(e) NULL)
  if (is.null(m0) || is.null(m1)) return(NA)
  lr <- 2 * (extract.mirt(m1, "logLik") - extract.mirt(m0, "logLik"))
  pchisq(lr, df = 4, lower.tail = FALSE) < DIF_ALPHA        # a + 3 thresholds
}

comparisons <- DIF_COMPARISONS
for (k in DIF_K_DS) comparisons[[sprintf("DS (%d) vs non-DS (%d)", 3 * k, 810 - 3 * k)]] <-
  c(810 - 3 * k, 3 * k)

res <- do.call(rbind, lapply(names(comparisons), function(nm) {
  n <- comparisons[[nm]]
  do.call(rbind, lapply(c(0, DIF_SHIFT), function(sh) {
    v <- unlist(parallel::mclapply(seq_len(DIF_REPS), function(i) dif_rep(n[1], n[2], sh),
                                   mc.cores = N_CORES))
    data.frame(comparison = nm, n_reference = n[1], n_focal = n[2], dif_shift = sh,
               rejection_rate = mean(v, na.rm = TRUE), failures = mean(is.na(v)))
  }))
}))
write.csv(res, file.path(RESULTS, "03_invariance_dif.csv"), row.names = FALSE)
print(res, digits = 3)
