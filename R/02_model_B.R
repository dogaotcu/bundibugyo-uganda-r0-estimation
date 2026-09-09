###############################################################################
## 02_model_B.R
##
## Model B: whole-period / border-closure model
## (Equations 2-3 in the manuscript). Before closure (t < t1), imported
## incidence follows the same exponential-growth traveler model as Model A.
## From closure day t1 onward, expected incidence is a mixture of:
##   (i)  pre-closure travelers whose symptoms onset after the incubation
##        period ("residual" pre-closure share), and
##   (ii) travelers who still cross at a reduced rate (1 - zeta) after
##        closure.
##
## log(i0) and r are estimated jointly by maximum-likelihood optimization
## (no closed-form profiling, unlike Model A). 95% CIs and simultaneous
## confidence bands are obtained by parametric bootstrap.
##
## Fitted for three alternative assumed closure days (12, 13, 14 -> 26,
## 27, 28 May 2026) for the Figure 3B sensitivity panel; the day-13 fit
## (27 May, the reported closure date) is the headline result and feeds
## Figure 2B.
###############################################################################

source("R/00_utils.R")

CLOSURE_DAYS <- 12:14           # assumed closure days under evaluation
MAIN_CLOSURE_DAY <- 13L         # reported closure date: 27 May 2026
ZETA <- 0.95                    # assumed border-closure effectiveness

INCUBATION_SHAPE <- 2.4
INCUBATION_MEAN  <- 10.5
INCUBATION_SCALE <- INCUBATION_MEAN / INCUBATION_SHAPE

LOG_I0_LOWER <- log(1e-8)
LOG_I0_UPPER <- log(1e12)
R_LOWER <- -1 / INCUBATION_SCALE + 1e-7
R_UPPER <- 1.0

observed <- load_observed_cases()
N_DAYS <- nrow(observed)
OBS_TIME <- observed$day - 1L
TRAVELER_FRACTION <- DAILY_TRAVELERS / (POPULATION * Q_DETECTION)


calculate_hybrid_basis <- function(r, closure_day) {
  if (length(r) != 1L || !is.finite(r) || r <= R_LOWER || r > R_UPPER) return(NULL)
  if (length(closure_day) != 1L || !is.finite(closure_day) ||
      closure_day < 2L || closure_day > N_DAYS) return(NULL)

  t1 <- closure_day - 1L  # first displayed day affected by closure
  beta_r <- 1 + INCUBATION_SCALE * r
  if (!is.finite(beta_r) || beta_r <= 0) return(NULL)
  tilted_scale <- INCUBATION_SCALE / beta_r

  G_t1 <- pgamma(t1, shape = INCUBATION_SHAPE, scale = tilted_scale)
  if (!is.finite(G_t1) || G_t1 <= 0) return(NULL)

  log_A <- log(TRAVELER_FRACTION)
  log_h_total <- preclosure_share <- postclosure_share <- rep(NA_real_, N_DAYS)

  for (j in seq_len(N_DAYS)) {
    T_now <- OBS_TIME[j]

    if (T_now < t1) {
      # Before closure: simple exponential traveler model, no convolution
      log_h_total[j] <- log_A + r * T_now
      preclosure_share[j] <- 1
      postclosure_share[j] <- 0
      next
    }

    # After closure: normalized incubation-period convolution
    dt <- T_now - t1
    survival_after <- pgamma(dt, shape = INCUBATION_SHAPE, scale = tilted_scale, lower.tail = FALSE)
    survival_T <- pgamma(T_now, shape = INCUBATION_SHAPE, scale = tilted_scale, lower.tail = FALSE)
    residual_pre <- max(survival_after - survival_T, 0)  # residual pre-closure onsets

    G_after <- pgamma(dt, shape = INCUBATION_SHAPE, scale = tilted_scale)
    post_passage <- (1 - ZETA) * G_after                 # reduced post-closure passage

    total <- residual_pre + post_passage
    if (!is.finite(total) || total <= 0) return(NULL)

    ratio <- total / G_t1
    if (!is.finite(ratio) || ratio <= 0) return(NULL)

    log_h_total[j] <- log_A + r * T_now + log(ratio)
    preclosure_share[j] <- residual_pre / total
    postclosure_share[j] <- post_passage / total
  }

  if (any(!is.finite(log_h_total)) || any(!is.finite(preclosure_share)) ||
      any(!is.finite(postclosure_share))) return(NULL)

  list(log_h_total = log_h_total, preclosure_share = preclosure_share,
       postclosure_share = postclosure_share, t1 = t1, tilted_scale = tilted_scale, G_t1 = G_t1)
}

# --------------------------------------------------------------------------
# Joint Poisson likelihood for (log_i0, r)
# --------------------------------------------------------------------------
joint_poisson_nll <- function(par, closure_day, y) {
  log_i0 <- as.numeric(par[1]); r <- as.numeric(par[2])
  if (!is.finite(log_i0) || !is.finite(r) ||
      log_i0 < LOG_I0_LOWER || log_i0 > LOG_I0_UPPER ||
      r <= R_LOWER || r > R_UPPER) return(1e100)

  basis <- calculate_hybrid_basis(r, closure_day)
  if (is.null(basis)) return(1e100)

  log_mu <- log_i0 + basis$log_h_total
  if (any(!is.finite(log_mu)) || any(log_mu > 700)) return(1e100)
  mu <- exp(log_mu)
  if (any(!is.finite(mu)) || any(mu <= 0)) return(1e100)

  nll <- -sum(dpois(y, lambda = mu, log = TRUE))
  if (!is.finite(nll)) return(1e100)
  nll
}

calculate_joint_model <- function(par, closure_day, y) {
  log_i0 <- as.numeric(par[1]); r <- as.numeric(par[2])
  basis <- calculate_hybrid_basis(r, closure_day)
  if (is.null(basis)) return(NULL)

  log_mu <- log_i0 + basis$log_h_total
  if (any(!is.finite(log_mu)) || any(log_mu > 700)) return(NULL)
  mu <- exp(log_mu)
  if (any(!is.finite(mu)) || any(mu <= 0)) return(NULL)

  loglik <- sum(dpois(y, lambda = mu, log = TRUE))
  if (!is.finite(loglik)) return(NULL)

  list(par = c(log_i0 = log_i0, r = r), log_i0 = log_i0, i0 = exp(log_i0), r = r,
       R0 = calculate_R0(r), mu_total = mu, cumulative_mu = cumsum(mu),
       mu_preclosure = mu * basis$preclosure_share, mu_postclosure = mu * basis$postclosure_share,
       loglik = loglik, nll = -loglik)
}

# --------------------------------------------------------------------------
# Optimization helpers: multi-start L-BFGS-B with box constraints
# --------------------------------------------------------------------------
clamp_parameter_vector <- function(par) {
  c(log_i0 = min(max(as.numeric(par[1]), LOG_I0_LOWER + 1e-8), LOG_I0_UPPER - 1e-8),
    r      = min(max(as.numeric(par[2]), R_LOWER + 1e-8),      R_UPPER - 1e-8))
}

run_joint_optim_once <- function(start_par, closure_day, y, maxit = 20000L) {
  start_par <- clamp_parameter_vector(start_par)
  tryCatch(
    optim(start_par, joint_poisson_nll, closure_day = closure_day, y = y,
          method = "L-BFGS-B",
          lower = c(log_i0 = LOG_I0_LOWER, r = R_LOWER),
          upper = c(log_i0 = LOG_I0_UPPER, r = R_UPPER),
          control = list(maxit = maxit, factr = 1e7, pgtol = 1e-10, lmm = 20)),
    error = function(e) NULL
  )
}

is_valid_optim_result <- function(result) {
  !is.null(result) && is.list(result) &&
    !is.null(result$convergence) && result$convergence == 0L &&
    !is.null(result$value) && is.finite(result$value) &&
    !is.null(result$par) && length(result$par) == 2L && all(is.finite(result$par))
}

make_start_from_r <- function(r_start, closure_day, y) {
  basis <- calculate_hybrid_basis(r_start, closure_day)
  if (is.null(basis)) return(NULL)
  h_start <- exp(basis$log_h_total)
  i0_start <- sum(y) / sum(h_start)  # crude scale estimate for a stable start
  if (!is.finite(i0_start) || i0_start <= 0) return(NULL)
  clamp_parameter_vector(c(log_i0 = log(i0_start), r = r_start))
}

ORIGINAL_R_STARTS <- c(-0.20, -0.15, -0.10, -0.05, 0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50)

fit_original_joint_model <- function(y, closure_day) {
  set.seed(123L + closure_day)
  r_starts <- ORIGINAL_R_STARTS[ORIGINAL_R_STARTS > R_LOWER & ORIGINAL_R_STARTS < R_UPPER]

  starts <- Filter(Negate(is.null), lapply(r_starts, make_start_from_r, closure_day = closure_day, y = y))
  if (length(starts) == 0L) return(NULL)

  fits <- lapply(starts, run_joint_optim_once, closure_day = closure_day, y = y)
  ok <- vapply(fits, is_valid_optim_result, logical(1))
  if (!any(ok)) return(NULL)
  fits <- fits[ok]
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]

  model_hat <- calculate_joint_model(best$par, closure_day, y)
  if (is.null(model_hat)) return(NULL)

  c(model_hat, list(convergence = best$convergence, message = best$message))
}

# --------------------------------------------------------------------------
# Parametric bootstrap:
# --------------------------------------------------------------------------
BOOTSTRAP_RETRY_STARTS <- 3L

run_one_parametric_bootstrap <- function(bootstrap_id, fit0, closure_day) {
  y_boot <- rpois(N_DAYS, lambda = fit0$mu_total)
  fail <- data.frame(bootstrap_id, closure_day, log_i0 = NA_real_, i0 = NA_real_,
                      r = NA_real_, R0 = NA_real_, convergence = FALSE)
  if (sum(y_boot) <= 0) return(fail)

  attempt <- run_joint_optim_once(fit0$par, closure_day, y_boot, maxit = 5000L)

  if (!is_valid_optim_result(attempt)) {
    for (retry in seq_len(BOOTSTRAP_RETRY_STARTS)) {
      retry_start <- clamp_parameter_vector(fit0$par + c(rnorm(1, 0, 0.08), rnorm(1, 0, 0.03)))
      retry_fit <- run_joint_optim_once(retry_start, closure_day, y_boot, maxit = 5000L)
      if (is_valid_optim_result(retry_fit)) { attempt <- retry_fit; break }
    }
  }
  if (!is_valid_optim_result(attempt)) return(fail)

  model_boot <- calculate_joint_model(attempt$par, closure_day, y_boot)
  if (is.null(model_boot)) return(fail)

  data.frame(bootstrap_id, closure_day, log_i0 = model_boot$log_i0, i0 = model_boot$i0,
             r = model_boot$r, R0 = model_boot$R0, convergence = TRUE)
}

generate_joint_bootstrap_mean_curves <- function(valid_bootstrap, closure_day) {
  if (nrow(valid_bootstrap) < 2L) stop("At least two successful joint-bootstrap estimates are required.")
  mat <- vapply(seq_len(nrow(valid_bootstrap)), function(b) {
    basis_b <- calculate_hybrid_basis(valid_bootstrap$r[b], closure_day)
    if (is.null(basis_b)) return(rep(NA_real_, N_DAYS))
    mu_b <- exp(valid_bootstrap$log_i0[b] + basis_b$log_h_total)
    if (any(!is.finite(mu_b)) || any(mu_b <= 0)) return(rep(NA_real_, N_DAYS))
    mu_b
  }, numeric(N_DAYS))
  ok <- apply(mat, 2, function(x) all(is.finite(x)) && all(x > 0))
  mat[, ok, drop = FALSE]
}

# --------------------------------------------------------------------------
# Fit + bootstrap for every assumed closure day
# --------------------------------------------------------------------------
run_model_B <- function() {
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  summary_list <- fitted_list <- bootstrap_list <- list()

  for (closure_day in CLOSURE_DAYS) {
    closure_date <- START_DATE + closure_day - 1L
    cat(sprintf("\n[Model B] assumed closure day %d (%s)\n", closure_day, format(closure_date, "%d %b %Y")))

    fit_k <- fit_original_joint_model(observed$daily_cases, closure_day)
    if (is.null(fit_k)) { warning("Fit failed for closure_day ", closure_day); next }

    boot_k <- future.apply::future_lapply(
      seq_len(N_BOOTSTRAP),
      function(b) run_one_parametric_bootstrap(b, fit_k, closure_day),
      future.seed = 123L + closure_day, future.scheduling = 2
    ) |> dplyr::bind_rows()
    bootstrap_list[[as.character(closure_day)]] <- boot_k

    valid_boot <- boot_k |>
      dplyr::filter(convergence, is.finite(i0), i0 > 0, is.finite(r), is.finite(R0))
    if (nrow(valid_boot) < 2L) { warning("Too few valid bootstrap fits for closure_day ", closure_day); next }

    daily_matrix <- generate_joint_bootstrap_mean_curves(valid_boot, closure_day)
    traj <- generate_poisson_confidence_trajectories(
      daily_matrix, N_CONFIDENCE_TRAJECTORIES, seed = 24680L + closure_day
    )
    daily_band <- calculate_simultaneous_confidence_band(fit_k$mu_total, traj$count_trajectory_matrix)

    cumulative_matrix <- apply(traj$count_trajectory_matrix, 2, cumsum)
    cumulative_band <- calculate_simultaneous_confidence_band(fit_k$cumulative_mu, cumulative_matrix)

    i0_CI <- quantile(valid_boot$i0, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    r_CI  <- quantile(valid_boot$r,  c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    R0_CI <- quantile(valid_boot$R0, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

    cat(sprintf("  i0 = %.2f | r = %.4f (95%% CI %.4f-%.4f) | R0 = %.2f (95%% CI %.2f-%.2f)\n",
                fit_k$i0, fit_k$r, r_CI[1], r_CI[2], fit_k$R0, R0_CI[1], R0_CI[2]))

    summary_list[[as.character(closure_day)]] <- data.frame(
      closure_day = closure_day, total_cases = sum(observed$daily_cases),
      i0_hat = fit_k$i0, i0_lower = i0_CI[1], i0_upper = i0_CI[2],
      r_hat = fit_k$r,   r_lower = r_CI[1],   r_upper = r_CI[2],
      R0_hat = fit_k$R0, R0_lower = R0_CI[1], R0_upper = R0_CI[2],
      logLik = fit_k$loglik, bootstrap_success = nrow(valid_boot), bootstrap_total = N_BOOTSTRAP
    )

    fitted_list[[as.character(closure_day)]] <- data.frame(
      closure_day = closure_day, day = observed$day, date = observed$date,
      observed_daily = observed$daily_cases,
      fitted_daily = fit_k$mu_total,
      fitted_lower = daily_band$lower,
      fitted_upper = daily_band$upper,
      fitted_cumulative = fit_k$cumulative_mu,
      fitted_cumulative_lower = cumulative_band$lower,
      fitted_cumulative_upper = cumulative_band$upper
    )
  }

  list(
    summary = dplyr::bind_rows(summary_list),
    fitted  = dplyr::bind_rows(fitted_list),
    bootstrap = dplyr::bind_rows(bootstrap_list),
    main_closure_day = MAIN_CLOSURE_DAY
  )
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  model_B_results <- run_model_B()
  dir.create("output", showWarnings = FALSE)
  saveRDS(model_B_results, "output/model_B_results.rds")
  cat("\n[Model B] saved output/model_B_results.rds\n")
  print(model_B_results$summary)
}
