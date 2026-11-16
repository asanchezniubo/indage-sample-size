# =============================================================================
# INDAGE — sample size simulation (AEI 2026 proposal, Section 3.2.3)
# 00_assumptions.R: every planning assumption in one place.
# Run all scripts from the repository root, e.g. `Rscript R/run_all.R`.
# =============================================================================
# Values are planning assumptions (no pilot data yet). Each block states its
# origin and what should be revised once pilot data are available.

set.seed(20261111)
N_CORES <- max(1, parallel::detectCores() - 1)
RESULTS <- "results"
dir.create(RESULTS, showWarnings = FALSE)

# ---- Target sample (Table 1 of Section 3.2.3) -------------------------------
# Counts by five-year band, 35-39 ... 70+ (70+ simulated as 70-79).
# The design has two samples that are modelled separately (Section 3.2.6), so
# their sizes are set separately: adults with Down syndrome (DS) are recruited
# between 35 and 64 only, because survival makes them scarce above 65.
BANDS      <- c(35, 40, 45, 50, 55, 60, 65, 70)          # lower limits
N_BAND     <- c(175, 170, 151, 144, 115, 106,  71,  68)  # whole sample, Table 1
N_DS_BAND  <- c( 57,  58,  57,  58,  35,  35,   0,   0)  # national DS minimum
N_NON_BAND <- N_BAND - N_DS_BAND                         # adults of other aetiologies
N_TOTAL    <- sum(N_BAND); N_DS <- sum(N_DS_BAND); N_NON <- sum(N_NON_BAND)
stopifnot(N_TOTAL == 1000, N_DS == 300, N_NON == 700, all(N_NON_BAND > 0))
# DS minimum per decade (35-44, 45-54, 55-64), as stated in Section 3.2.3.
DS_PER_DECADE <- c(115, 115, 70)
stopifnot(sum(DS_PER_DECADE) == N_DS)

# ---- Criterion 1: item parameter estimation (graded response model) ---------
# One domain = 5 ordinal items with 4 categories; discriminations 1.0-2.2.
# Thresholds shifted towards the "no difficulty" end (mild-moderate ID, 35+):
# a moderate ceiling, as expected in self-reported functional ability.
GRM_N    <- c(200, 300, 500, 700, 1000)
GRM_REPS <- 300
GRM_A    <- c(1.0, 1.3, 1.6, 1.9, 2.2)
GRM_B    <- rbind(c(-2.2, -1.2, 0.0),   # b1 < b2 < b3 (theta scale)
                  c(-2.0, -1.0, 0.3),
                  c(-1.8, -0.8, 0.5),
                  c(-1.6, -0.5, 0.8),
                  c(-1.4, -0.3, 1.0))
GRM_MAX_RMSE_A <- 0.25                  # planning acceptance criteria
GRM_MAX_RMSE_B <- 0.25

# ---- Criterion 2: factor structure (logic of SENECA; Lorenzo-Seva & Ferrando, 2024)
# One part of the scale = 5 factors (domains). Data are generated from an
# assumed population model, because the original procedure starts from a real
# initial sample (here: the first 250 nationwide participants). CRMSR = discrepancy
# between population (L Phi L') and sample (A A') reproduced common variance.
SEN_N      <- c(300, 400, 500, 600, 700, 800, 1000, 1200)
SEN_REPS   <- 150
SEN_CATS   <- 4
# Intrinsic capacity: 5 domains with 6/3/4/6/4 indicators (cognition, locomotor,
# sensory, psychological, vitality); moderate loadings.
SEN_IC <- list(n_items = c(6, 3, 4, 6, 4),
               loadings = list(c(.70,.65,.60,.60,.55,.50), c(.65,.60,.55),
                               c(.60,.55,.55,.50), c(.70,.65,.60,.60,.55,.50),
                               c(.60,.55,.50,.45)),
               phi = 0.40)                  # correlation between domains
# Functional ability: 5 domains x 5 items; somewhat higher loadings (ADL/IADL).
SEN_FA <- list(n_items = rep(5, 5),
               loadings = rep(list(c(.80,.75,.70,.65,.60)), 5),
               phi = 0.50)

# ---- Criterion 3: invariance (power to detect DIF) --------------------------
# 5-item GRM domain; one item with uniform DIF (shifted thresholds) in the focal
# group; likelihood-ratio test in mirt::multipleGroup with the other 4 as anchors.
DIF_REPS   <- 150
DIF_SHIFT  <- c(0.5)        # threshold shift (theta); 0 = type I error
DIF_ALPHA  <- 0.05
DIF_K_DS   <- c(60, 100, 150) # DS per decade contrasted here
DIF_COMPARISONS <- list(
  "65+ (139) vs 35-44 (345)"   = c(sum(N_BAND[1:2]), sum(N_BAND[7:8])),
  "55-64 (221) vs 35-44 (345)" = c(sum(N_BAND[1:2]), sum(N_BAND[5:6])))

# ---- Criterion 4: precision of the onset age (binding) ----------------------
# Domain score (T-score relative to the 35-44 group of each stratum): flat up
# to A0, then linear decline of B SD per decade. Planning assumptions:
#   non-DS: A0 = 50, B = 0.35  -> 0.5 SD drop around age 64
#   DS    : A0 = 42, B = 0.70  -> 0.5 SD drop around age 49
# Conservative scenario (slower decline in DS): B_DS = 0.50.
ONSET_SD_DROP   <- 0.5      # definition of onset age (Section 3.2.6)
ONSET_REL       <- 0.85     # reliability of domain scores
ONSET_DF_SPLINE <- 3
ONSET_REPS      <- 16000   # The half-width is a 2.5/97.5 quantile of the estimate
                            # distribution, so few replications underestimate it:
                            # 1,000 give +-5.0 years in the binding cell and 4,000
                            # give 4.95 and 5.15 on two seeds. 16,000 pin it down.
ONSET_SCENARIOS <- list(
  base         = list(A0_non = 50, B_non = 0.35, A0_ds = 42, B_ds = 0.70),
  conservative = list(A0_non = 50, B_non = 0.35, A0_ds = 42, B_ds = 0.50))
ONSET_N_DS      <- c(180, 240, 300, 360)   # DS sample sizes contrasted
ONSET_N_NON     <- c(600, 700, 800)        # non-DS sample sizes contrasted
ONSET_MAX_HALFWIDTH <- 5    # target: 95% interval of DS onset age within +-5 years
ONSET_MIN_POWER     <- 0.80

# ---- Allocation profiles compared in 05_allocation.R ------------------------
# Relative weights across the eight five-year bands. The adopted profile keeps
# the 35-44 reference large, because every onset age is defined relative to its
# mean, and supports the upper boundary, because natural splines are linear
# beyond the boundary knots.
ALLOC_PROFILES <- list(
  "2025 proposal"        = c(153, 153, 153, 132, 119, 100,  55,  55),
  "proportional to frame"= c(124, 124, 124, 124, 124, 126,  90,  84),
  "uniform"              = rep(115, 8),
  "adopted"              = c(150, 150, 125, 115, 105,  95,  90,  90),
  "U-shaped"             = c(160, 160,  90,  90,  90,  90, 120, 120))
# Splits of the DS minimum across the three decades, at a fixed total.
DS_SPLITS <- list("115/115/70 (adopted)" = c(115, 115, 70),
                  "100/100/100"          = c(100, 100, 100),
                  "120/120/60"           = c(120, 120, 60),
                  "140/110/50"           = c(140, 110, 50),
                  "100/120/80"           = c(100, 120, 80))

# Quick test: SMOKE=1 Rscript R/0X_*.R
if (Sys.getenv("SMOKE") == "1") {
  GRM_REPS <- 4; SEN_REPS <- 4; DIF_REPS <- 4; ONSET_REPS <- 20
  GRM_N <- c(300, 1000); SEN_N <- c(400, 1000)
  ONSET_N_DS <- c(180, 300); ONSET_N_NON <- c(700)
  ALLOC_PROFILES <- ALLOC_PROFILES[c(1, 4)]; DS_SPLITS <- DS_SPLITS[c(1, 3)]
  DIF_COMPARISONS <- DIF_COMPARISONS[1]; DIF_K_DS <- 60
}
