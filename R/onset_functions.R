# Functions shared by 04_onset_age.R, 05_allocation.R and 07_design_grid.R.
suppressPackageStartupMessages(library(splines))

ERR_VAR <- 1 / ONSET_REL - 1                       # measurement error (latent var = 1)
mu <- function(age, A0, B) ifelse(age < A0, 0, -B * (age - A0) / 10)

# Onset age = youngest age at which the fitted mean (natural spline) falls
# ONSET_SD_DROP SDs below the mean of the 35-44 group; searched from age 45.
onset_age <- function(age, y, range) {
  fit <- lm(y ~ ns(age, df = ONSET_DF_SPLINE))
  ref <- age >= 35 & age < 45
  threshold <- mean(fitted(fit)[ref]) - ONSET_SD_DROP * sd(y[ref])
  grid <- seq(range[1], range[2], by = 0.1)
  pr <- predict(fit, newdata = data.frame(age = grid))
  i <- which(pr < threshold & grid >= 45)[1]
  if (is.na(i)) NA else grid[i]
}

# Same definition applied to the population curve.
true_onset <- function(A0, B, range) {
  age <- seq(35, range[2], by = 0.01); m <- mu(age, A0, B)
  threshold <- mean(m[age < 45]) - ONSET_SD_DROP * sqrt(1 + ERR_VAR)
  age[which(m < threshold & age >= 45)[1]]
}

sample_ages <- function(n_per_band5, last_width = 10) {
  unlist(mapply(function(lo, n, w) runif(n, lo, lo + w), BANDS, n_per_band5,
                c(rep(5, length(BANDS) - 1), last_width), SIMPLIFY = FALSE))
}

# Largest-remainder allocation of n units across relative weights w.
largest_remainder <- function(n, w) {
  raw <- n * w / sum(w); base <- floor(raw); k <- n - sum(base)
  if (k > 0) { o <- order(-(raw - base))[seq_len(k)]; base[o] <- base[o] + 1 }
  stopifnot(sum(base) == n)
  base
}

# Summary of one set of replications of the two onset ages.
summarise_onset <- function(r, true_ds, true_non, meta) {
  d <- r[, "non"] - r[, "ds"]; se <- sd(d, na.rm = TRUE)
  q_ds  <- quantile(r[, "ds"],  c(.025, .975), na.rm = TRUE)
  q_non <- quantile(r[, "non"], c(.025, .975), na.rm = TRUE)
  cbind(meta, data.frame(
    true_DS = round(true_ds, 1), mean_DS = mean(r[, "ds"], na.rm = TRUE),
    halfwidth_DS = unname(diff(q_ds) / 2), undetected_DS = mean(is.na(r[, "ds"])),
    true_non_DS = round(true_non, 1), mean_non_DS = mean(r[, "non"], na.rm = TRUE),
    halfwidth_non_DS = unname(diff(q_non) / 2), undetected_non_DS = mean(is.na(r[, "non"])),
    mean_difference = mean(d, na.rm = TRUE),
    power_difference = mean(!is.na(d) & d / se > qnorm(0.975))))
}
