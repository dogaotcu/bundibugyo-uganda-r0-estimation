###############################################################################
## 01_model_A.R
##
## Model A: exponential-growth traveler model fitted to the pre-closure
## period (Equation 1 in the manuscript):
##
##   E[c_t] = (N / (m*q)) * i0 * exp(r*t)
##   c_t ~ Poisson(E[c_t])
##
## For any candidate r, the conditional MLE of i0 is available in closed
## form, so only r is optimized numerically (profile likelihood). 95% CIs
## and simultaneous confidence bands are obtained by parametric bootstrap.
##
## Fitted for three alternative analysis end days (22, 23, 24 -> 5, 6,
## 7 June 2026) to produce the Figure 3A end-date sensitivity panel; the
## day-23 fit is the model's headline result and feeds Figure 2A.
###############################################################################

source("R/00_utils.R")

END_DAYS <- 22:24              # analysis end days under evaluation
MAIN_END_DAY <- 23L            # headline result reported in the abstract

observed <- load_observed_cases()
if (max(END_DAYS) > nrow(observed)) stop("END_DAYS exceeds the available data.")

# --------------------------------------------------------------------------
# Build the analysis data frame truncated at a given end day
# --------------------------------------------------------------------------
make_analysis_data <- function(end_day) {
  d <- observed[seq_len(end_day), ]
  data.frame(
    day = d$day,
    time = d$day - 1L,
    daily_cases = d$daily_cases,
    cumulative_cases = d$cumulative_cases,
    N_t = rep(DAILY_TRAVELERS, end_day)
  )
}

# --------------------------------------------------------------------------
# Profile likelihood: for a given r, i0 has a closed-form conditional MLE
# --------------------------------------------------------------------------
calculate_profile_values <- function(r, data) {
  if (length(r) != 1L || !is.finite(r)) return(NULL)

  total_cases <- sum(data$daily_cases)
  if (!is.finite(total_cases) || total_cases <= 0) return(NULL)

  log_exposure <- log(data$N_t / (POPULATION * Q_DETECTION))
  log_h <- log_exposure + r * data$time
  if (any(!is.finite(log_h))) return(NULL)

  # log-sum-exp for numerical stability
  log_h_max <- max(log_h)
  log_sum_h <- log_h_max + log(sum(exp(log_h - log_h_max)))
  log_i0 <- log(total_cases) - log_sum_h

  log_mu <- log_i0 + log_h
  if (any(!is.finite(log_mu)) || any(log_mu > 700)) return(NULL)

  mu <- exp(log_mu)
  if (any(!is.finite(mu)) || any(mu <= 0)) return(NULL)

  list(i0 = exp(log_i0), log_i0 = log_i0, r = r, mu_t = mu, cumulative_mu = cumsum(mu))
}

profile_poisson_nll <- function(par, data) {
  model <- calculate_profile_values(r = as.numeric(par[1]), data = data)
  if (is.null(model)) return(1e100)
  nll <- -sum(dpois(data$daily_cases, lambda = model$mu_t, log = TRUE))
  if (!is.finite(nll)) return(1e100)
  nll
}

# --------------------------------------------------------------------------
# Fit r by multi-start BFGS, then recover i0 analytically
# --------------------------------------------------------------------------
fit_traveler_poisson <- function(data, starting_r = NULL) {
  if (any(!is.finite(data$daily_cases)) || any(data$daily_cases < 0) ||
      sum(data$daily_cases) <= 0) return(NULL)

  r_starts <- if (is.null(starting_r) || !is.finite(starting_r)) {
    c(-1, -0.5, -0.3, -0.2, -0.1, -0.05, 0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1)
  } else {
    unique(c(starting_r, starting_r + c(-0.2, -0.1, -0.05, 0.05, 0.1, 0.2), 0))
  }

  fits <- lapply(r_starts, function(r0) {
    tryCatch(
      optim(c(r = r0), profile_poisson_nll, data = data, method = "BFGS",
            control = list(maxit = 20000L, reltol = 1e-12)),
      error = function(e) NULL
    )
  })

  ok <- vapply(fits, function(x) !is.null(x) && is.finite(x$value) && all(is.finite(x$par)), logical(1))
  if (!any(ok)) return(NULL)
  fits <- fits[ok]
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]

  r_hat <- as.numeric(best$par[1])
  model <- calculate_profile_values(r = r_hat, data = data)
  if (is.null(model)) return(NULL)

  list(
    i0 = model$i0, log_i0 = model$log_i0, r = r_hat,
    R0 = calculate_R0(r_hat),
    mu_t = model$mu_t, cumulative_mu = model$cumulative_mu,
    nll = best$value, loglik = -best$value,
    convergence = best$convergence, message = best$message
  )
}

calculate_curve_from_parameters <- function(i0, r, data) {
  if (!is.finite(i0) || i0 <= 0 || !is.finite(r)) return(NULL)
  log_mu <- log(i0) + log(data$N_t / (POPULATION * Q_DETECTION)) + r * data$time
  if (any(!is.finite(log_mu)) || any(log_mu > 700)) return(NULL)
  mu <- exp(log_mu)
  if (any(!is.finite(mu)) || any(mu <= 0)) return(NULL)
  list(daily = mu, cumulative = cumsum(mu))
}

run_one_bootstrap <- function(bootstrap_id, fit0, data) {
  y_boot <- rpois(nrow(data), lambda = fit0$mu_t)
  boot_data <- data
  boot_data$daily_cases <- y_boot
  boot_data$cumulative_cases <- cumsum(y_boot)

  fit_boot <- fit_traveler_poisson(boot_data, starting_r = fit0$r)
  if (is.null(fit_boot) || fit_boot$convergence != 0L ||
      !is.finite(fit_boot$i0) || !is.finite(fit_boot$r)) {
    return(data.frame(bootstrap_id, i0 = NA_real_, log_i0 = NA_real_,
                       r = NA_real_, R0 = NA_real_, convergence = FALSE))
  }
  data.frame(bootstrap_id, i0 = fit_boot$i0, log_i0 = fit_boot$log_i0,
             r = fit_boot$r, R0 = fit_boot$R0, convergence = TRUE)
}

generate_joint_bootstrap_mean_curves <- function(valid_bootstrap, data) {
  mat <- vapply(seq_len(nrow(valid_bootstrap)), function(b) {
    curve <- calculate_curve_from_parameters(valid_bootstrap$i0[b], valid_bootstrap$r[b], data)
    if (is.null(curve)) return(rep(NA_real_, nrow(data)))
    curve$daily
  }, numeric(nrow(data)))
  ok <- apply(mat, 2, function(x) all(is.finite(x)) && all(x > 0))
  mat[, ok, drop = FALSE]
}

# --------------------------------------------------------------------------
# Fit + bootstrap for every analysis end day
# --------------------------------------------------------------------------
run_model_A <- function() {
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  summary_list <- list()
  fitted_list  <- list()
  bootstrap_list <- list()

  for (end_day in END_DAYS) {
    cat(sprintf("\n[Model A] fitting through day %d (%s)\n",
                 end_day, format(START_DATE + end_day - 1, "%d %b %Y")))

    data_k <- make_analysis_data(end_day)
    fit_k <- fit_traveler_poisson(data_k)
    if (is.null(fit_k)) { warning("Fit failed for end_day ", end_day); next }

    boot_k <- future.apply::future_lapply(
      seq_len(N_BOOTSTRAP),
      function(b) run_one_bootstrap(b, fit_k, data_k),
      future.seed = 123L + end_day, future.scheduling = 2
    ) |> dplyr::bind_rows()
    boot_k$end_day <- end_day
    bootstrap_list[[as.character(end_day)]] <- boot_k

    valid_boot <- boot_k |>
      dplyr::filter(convergence, is.finite(i0), i0 > 0, is.finite(r), is.finite(R0))
    if (nrow(valid_boot) < 2L) { warning("Too few valid bootstrap fits for end_day ", end_day); next }

    daily_matrix <- generate_joint_bootstrap_mean_curves(valid_boot, data_k)
    traj <- generate_poisson_confidence_trajectories(
      daily_matrix, N_CONFIDENCE_TRAJECTORIES, seed = 123L + end_day
    )
    daily_band <- calculate_simultaneous_confidence_band(fit_k$mu_t, traj$count_trajectory_matrix)

    cumulative_matrix <- apply(traj$count_trajectory_matrix, 2, cumsum)
    cumulative_band <- calculate_simultaneous_confidence_band(fit_k$cumulative_mu, cumulative_matrix)

    i0_CI <- quantile(valid_boot$i0, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    r_CI  <- quantile(valid_boot$r,  c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    R0_CI <- quantile(valid_boot$R0, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

    cat(sprintf("  i0 = %.2f | r = %.4f (95%% CI %.4f-%.4f) | R0 = %.2f (95%% CI %.2f-%.2f)\n",
                fit_k$i0, fit_k$r, r_CI[1], r_CI[2], fit_k$R0, R0_CI[1], R0_CI[2]))

    summary_list[[as.character(end_day)]] <- data.frame(
      end_day = end_day, total_cases = sum(data_k$daily_cases),
      i0_hat = fit_k$i0, i0_lower = i0_CI[1], i0_upper = i0_CI[2],
      r_hat = fit_k$r,   r_lower = r_CI[1],   r_upper = r_CI[2],
      R0_hat = fit_k$R0, R0_lower = R0_CI[1], R0_upper = R0_CI[2],
      logLik = fit_k$loglik, bootstrap_success = nrow(valid_boot), bootstrap_total = N_BOOTSTRAP
    )

    fitted_list[[as.character(end_day)]] <- data.frame(
      end_day = end_day, day = data_k$day,
      observed_daily = data_k$daily_cases,
      fitted_daily = fit_k$mu_t,
      fitted_daily_lower = daily_band$lower,
      fitted_daily_upper = daily_band$upper,
      observed_cumulative = data_k$cumulative_cases,
      fitted_cumulative = fit_k$cumulative_mu,
      fitted_cumulative_lower = cumulative_band$lower,
      fitted_cumulative_upper = cumulative_band$upper
    )
  }

  list(
    summary  = dplyr::bind_rows(summary_list),
    fitted   = dplyr::bind_rows(fitted_list),
    bootstrap = dplyr::bind_rows(bootstrap_list),
    main_end_day = MAIN_END_DAY
  )
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  model_A_results <- run_model_A()
  dir.create("output", showWarnings = FALSE)
  saveRDS(model_A_results, "output/model_A_results.rds")
  cat("\n[Model A] saved output/model_A_results.rds\n")
  print(model_A_results$summary)
}
