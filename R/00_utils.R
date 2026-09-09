###############################################################################
## 00_utils.R
##
## Shared setup, publication theme, and bootstrap/confidence-band helpers
## used by both Model A (R/01_model_A.R) and Model B (R/02_model_B.R).
##
## These functions were previously duplicated verbatim in both model
## scripts; they are defined once here and sourced by everything else.
###############################################################################

# --------------------------------------------------------------------------
# Packages
# --------------------------------------------------------------------------
required_packages <- c(
  "dplyr", "tibble", "ggplot2", "future", "future.apply",
  "parallelly", "scales", "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(future)
  library(future.apply)
  library(parallelly)
  library(scales)
  library(patchwork)
})

if (packageVersion("ggplot2") < package_version("3.5.0")) {
  stop("ggplot2 >= 3.5.0 is required. Run install.packages('ggplot2').")
}

# --------------------------------------------------------------------------
# Shared study constants (Table 1 of the manuscript)
# --------------------------------------------------------------------------
START_DATE      <- as.Date("2026-05-15")  # day 1 of the observation period
POPULATION      <- 48656601               # DRC population (UN WPP)
Q_DETECTION     <- 1.0                    # probability of detection after import
ANNUAL_TRAVELERS <- 2160000                # DRC -> Uganda, mobility dataset
DAILY_TRAVELERS  <- ANNUAL_TRAVELERS / 365 # ~5,918 travelers/day

GI_MEAN <- 13.7   # generation-interval mean (days), Zaire ebolavirus proxy
GI_SD   <- 4.5    # generation-interval SD (days)

N_BOOTSTRAP               <- 1000L   # parametric-bootstrap replications
N_CONFIDENCE_TRAJECTORIES <- 10000L  # Poisson trajectories for the CI band
SIMULTANEOUS_CI_LEVEL     <- 0.95

# GI parameter ranges used in the Figure 4 sensitivity heatmaps
GI_MEAN_RANGE <- seq(11.1, 19.4, length.out = 10)
GI_SD_RANGE   <- seq(2.3, 9.3, length.out = 8)

# Published R0 ranges shown as reference bands in Figure 3
PUBLISHED_R0_ALL_LOWER    <- 1.10
PUBLISHED_R0_ALL_UPPER    <- 10.00
PUBLISHED_R0_UGANDA_LOWER <- 1.34
PUBLISHED_R0_UGANDA_UPPER <- 2.70

set.seed(123)

# --------------------------------------------------------------------------
# Parallel backend
# --------------------------------------------------------------------------
n_workers <- max(1L, as.integer(parallelly::availableCores()[1]) - 1L, na.rm = TRUE)
if (is.na(n_workers) || n_workers < 1L) n_workers <- 1L

# --------------------------------------------------------------------------
# Publication ggplot theme, shared by every figure script
# --------------------------------------------------------------------------
theme_publication <- theme_classic(base_size = 20, base_family = "Arial") +
  theme(
    text          = element_text(color = "black"),
    plot.title    = element_blank(),
    axis.title    = element_text(size = 18, face = "bold"),
    axis.text     = element_text(size = 15, color = "black"),
    axis.ticks    = element_line(color = "black"),
    axis.line     = element_line(color = "black", linewidth = 1),
    legend.text   = element_text(size = 14),
    legend.title  = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    strip.text    = element_text(color = "black", face = "bold", size = 18)
  )

# --------------------------------------------------------------------------
# Load the observed imported-case time series (single source of truth)
# --------------------------------------------------------------------------
load_observed_cases <- function(path = "data/imported_cases_uganda.csv") {
  observed <- read.csv(path, stringsAsFactors = FALSE)
  observed$date <- as.Date(observed$date)
  observed
}

format_calendar_date <- function(x) {
  x <- as.Date(x, origin = "1970-01-01")
  paste(as.integer(format(x, "%d")), month.abb[as.integer(format(x, "%m"))])
}

make_weekly_breaks <- function(start_date, end_date, by = "7 days") {
  seq.Date(as.Date(start_date), as.Date(end_date), by = by)
}

# --------------------------------------------------------------------------
# R0 from the exponential growth rate r via the Euler-Lotka equation,
# assuming a gamma-distributed generation interval (Eq. 6 in the manuscript)
# --------------------------------------------------------------------------
calculate_R0 <- function(r, GI_mean = GI_MEAN, GI_sd = GI_SD) {
  gamma_shape <- GI_mean^2 / GI_sd^2
  gamma_scale <- GI_sd^2 / GI_mean

  base <- 1 + gamma_scale * r
  result <- rep(NA_real_, length(r))
  valid <- is.finite(base) & base > 0
  result[valid] <- base[valid]^gamma_shape
  result
}

# --------------------------------------------------------------------------
# Parametric-bootstrap confidence-band construction
#
# Shared by Model A and Model B. Both models:
#   1. Simulate Y_t^(b) ~ Poisson(mu_hat_t) from the fitted mean curve
#   2. Refit the model to each bootstrap replicate, keeping the joint
#      (i0, r) pair
#   3. Reconstruct one fitted mean curve per successful bootstrap replicate
#   4. Sample those curves and draw a final Poisson count trajectory from
#      each, to add measurement/sampling variation on top of parameter
#      uncertainty
#   5. Build a *simultaneous* 95% band from the standardized maximum
#      absolute deviation over the fitted period
# --------------------------------------------------------------------------

#' Generate one Poisson count trajectory per sampled bootstrap mean curve.
generate_poisson_confidence_trajectories <- function(bootstrap_mean_matrix,
                                                       n_simulation, seed) {
  bootstrap_mean_matrix <- as.matrix(bootstrap_mean_matrix)
  storage.mode(bootstrap_mean_matrix) <- "double"

  if (nrow(bootstrap_mean_matrix) < 1L || ncol(bootstrap_mean_matrix) < 2L) {
    stop("At least two valid joint-bootstrap mean curves are required.")
  }
  if (any(!is.finite(bootstrap_mean_matrix)) || any(bootstrap_mean_matrix < 0)) {
    stop("bootstrap_mean_matrix contains invalid Poisson means.")
  }

  set.seed(seed)
  sampled_curve_index <- sample.int(ncol(bootstrap_mean_matrix), n_simulation, replace = TRUE)

  count_trajectory_matrix <- vapply(
    seq_len(n_simulation),
    function(s) {
      lambda_s <- bootstrap_mean_matrix[, sampled_curve_index[s]]
      stats::rpois(length(lambda_s), lambda = lambda_s)
    },
    numeric(nrow(bootstrap_mean_matrix))
  )

  list(count_trajectory_matrix = count_trajectory_matrix,
       sampled_curve_index = sampled_curve_index)
}

#' Simultaneous 95% confidence band from the standardized maximum absolute
#' deviation of bootstrap trajectories around the fitted mean curve.
calculate_simultaneous_confidence_band <- function(fitted_mean, bootstrap_mean_matrix,
                                                     level = SIMULTANEOUS_CI_LEVEL) {
  fitted_mean <- as.numeric(fitted_mean)
  bootstrap_mean_matrix <- as.matrix(bootstrap_mean_matrix)
  storage.mode(bootstrap_mean_matrix) <- "double"

  if (length(fitted_mean) != nrow(bootstrap_mean_matrix)) {
    stop(sprintf(
      "Length mismatch: fitted_mean has %d elements, bootstrap_mean_matrix has %d rows.",
      length(fitted_mean), nrow(bootstrap_mean_matrix)
    ))
  }
  if (any(!is.finite(fitted_mean))) stop("fitted_mean contains non-finite values.")
  if (ncol(bootstrap_mean_matrix) < 2L) stop("At least two valid bootstrap mean curves are required.")

  pointwise_sd <- apply(bootstrap_mean_matrix, 1, sd, na.rm = TRUE)
  sd_floor <- max(sqrt(.Machine$double.eps), 1e-10 * max(abs(fitted_mean), na.rm = TRUE))
  pointwise_sd_safe <- pmax(pointwise_sd, sd_floor)

  standardized <- sweep(bootstrap_mean_matrix, 1, fitted_mean, "-")
  standardized <- sweep(standardized, 1, pointwise_sd_safe, "/")
  max_abs_deviation <- apply(abs(standardized), 2, max, na.rm = TRUE)

  critical_value <- as.numeric(
    quantile(max_abs_deviation, probs = level, na.rm = TRUE, names = FALSE, type = 8)
  )

  list(
    lower = pmax(0, fitted_mean - critical_value * pointwise_sd_safe),
    upper = fitted_mean + critical_value * pointwise_sd_safe,
    pointwise_sd = pointwise_sd,
    critical_value = critical_value,
    maximum_absolute_deviation = max_abs_deviation
  )
}
