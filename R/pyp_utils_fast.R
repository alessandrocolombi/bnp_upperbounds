# Fast Pitman-Yor process utilities (R port of pyp_utils_fast.py)

new_fit_result <- function(alpha, theta, loglik, converged, meta = list()) {
  structure(
    list(
      alpha = as.numeric(alpha),
      theta = as.numeric(theta),
      loglik = as.numeric(loglik),
      converged = as.logical(converged),
      meta = meta
    ),
    class = "fit_result"
  )
}

print.fit_result <- function(x, ...) {
  cat(
    sprintf(
      "FitResult(alpha=%.10g, theta=%.10g, loglik=%.10g, converged=%s)\n",
      x$alpha,
      x$theta,
      x$loglik,
      ifelse(isTRUE(x$converged), "TRUE", "FALSE")
    )
  )
  if (length(x$meta) > 0L) {
    cat("meta:\n")
    print(x$meta)
  }
  invisible(x)
}

occupancy_counts <- function(x) {
  x <- as.vector(x)
  if (length(x) == 0L) {
    return(integer(0))
  }
  table(x)
}

sizes_from_data <- function(x, sort = TRUE) {
  counts <- occupancy_counts(x)
  sizes <- as.integer(unname(counts))
  if (isTRUE(sort)) {
    sizes <- base::sort(sizes, decreasing = TRUE)
  }
  sizes
}

mrn_over_kn <- function(sizes, rmax = 10L) {
  sizes <- as.integer(sizes)
  K <- length(sizes)
  out <- numeric(rmax)
  if (K == 0L) {
    return(out)
  }
  for (r in seq_len(rmax)) {
    out[r] <- mean(sizes == r)
  }
  out
}

.make_pyp_context <- function(sizes) {
  sizes_num <- as.numeric(sizes)
  unique_sizes <- base::sort(unique(sizes_num))
  size_counts <- tabulate(match(sizes_num, unique_sizes), nbins = length(unique_sizes))

  list(
    sizes = sizes_num,
    n = as.integer(sum(sizes_num)),
    K = as.integer(length(sizes_num)),
    sorted_sizes = base::sort(sizes_num, decreasing = TRUE),
    size_unique = unique_sizes,
    size_counts = as.numeric(size_counts)
  )
}

.ctx_term_sizes <- function(ctx, alpha) {
  if (length(ctx$size_unique) < floor(ctx$K / 2)) {
    return(sum(ctx$size_counts * lgamma(ctx$size_unique - alpha)) - ctx$K * lgamma(1.0 - alpha))
  }
  sum(lgamma(ctx$sorted_sizes - alpha)) - ctx$K * lgamma(1.0 - alpha)
}

.ctx_dalpha_term_sizes <- function(ctx, alpha) {
  if (length(ctx$size_unique) < floor(ctx$K / 2)) {
    return(-sum(ctx$size_counts * digamma(ctx$size_unique - alpha)) + ctx$K * digamma(1.0 - alpha))
  }
  -sum(digamma(ctx$sorted_sizes - alpha)) + ctx$K * digamma(1.0 - alpha)
}

.dlogl_dtheta_vectorized <- function(alpha, thetas, n, K) {
  if (K <= 1L) {
    return(digamma(thetas + 1.0) - digamma(thetas + n))
  }
  u <- thetas / alpha
  (digamma(u + K) - digamma(u + 1.0)) / alpha + digamma(thetas + 1.0) - digamma(thetas + n)
}

.logl_theta_term_vectorized <- function(alpha, thetas, K) {
  if (K <= 1L) {
    return(rep(0.0, length(thetas)))
  }
  u <- thetas / alpha
  (K - 1.0) * log(alpha) + lgamma(u + K) - lgamma(u + 1.0)
}

.compute_loglik <- function(alpha, theta, n, K, term_sizes) {
  term_theta <- 0.0
  if (K > 1L) {
    u <- theta / alpha
    term_theta <- (K - 1.0) * log(alpha) + lgamma(u + K) - lgamma(u + 1.0)
  }
  term_norm <- lgamma(theta + 1.0) - lgamma(theta + n)
  term_sizes + term_theta + term_norm
}

.validate_params <- function(alpha, theta, sizes, eps_alpha = 1e-8, eps_theta = 1e-12) {
  if (!(alpha > eps_alpha && alpha < (1.0 - eps_alpha))) {
    return(FALSE)
  }
  if (!(theta > eps_theta)) {
    return(FALSE)
  }
  if (any(sizes - alpha <= 0.0)) {
    return(FALSE)
  }
  TRUE
}

log_eppf_pyp <- function(alpha, theta, sizes, eps_alpha = 1e-8, eps_theta = 1e-12) {
  sizes <- as.numeric(sizes)
  n <- sum(sizes)
  K <- length(sizes)

  if (!.validate_params(alpha, theta, sizes, eps_alpha, eps_theta)) {
    return(-Inf)
  }

  term_sizes <- sum(lgamma(sizes - alpha)) - K * lgamma(1.0 - alpha)

  term_theta <- 0.0
  if (K > 1L) {
    u <- theta / alpha
    term_theta <- (K - 1.0) * log(alpha) + lgamma(u + K) - lgamma(u + 1.0)
  }

  term_norm <- lgamma(theta + 1.0) - lgamma(theta + n)
  as.numeric(term_sizes + term_theta + term_norm)
}

grad_log_eppf_pyp <- function(alpha, theta, sizes, eps_alpha = 1e-8, eps_theta = 1e-12) {
  sizes <- as.numeric(sizes)
  n <- sum(sizes)
  K <- length(sizes)

  if (!.validate_params(alpha, theta, sizes, eps_alpha, eps_theta)) {
    return(c(d_alpha = NaN, d_theta = NaN))
  }

  d_alpha_sizes <- -sum(digamma(sizes - alpha)) + K * digamma(1.0 - alpha)

  d_theta_theta <- 0.0
  d_alpha_theta <- 0.0
  if (K > 1L) {
    u <- theta / alpha
    psi_diff <- digamma(u + K) - digamma(u + 1.0)
    d_theta_theta <- psi_diff / alpha
    d_alpha_theta <- (K - 1.0) / alpha - (theta / (alpha^2)) * psi_diff
  }

  d_theta_norm <- digamma(theta + 1.0) - digamma(theta + n)

  c(
    d_alpha = as.numeric(d_alpha_sizes + d_alpha_theta),
    d_theta = as.numeric(d_theta_theta + d_theta_norm)
  )
}

.find_optimal_theta <- function(alpha, ctx, theta_lo, theta_hi, n_grid = 32L) {
  n <- ctx$n
  K <- ctx$K

  term_sizes <- .ctx_term_sizes(ctx, alpha)

  log_thetas <- seq(log10(theta_lo), log10(theta_hi), length.out = n_grid)
  thetas <- 10^log_thetas

  dvals <- .dlogl_dtheta_vectorized(alpha, thetas, n, K)

  sign_prod <- dvals[-length(dvals)] * dvals[-1]
  sign_changes <- which(sign_prod < 0.0)

  if (length(sign_changes) == 0L) {
    lls <- term_sizes +
      .logl_theta_term_vectorized(alpha, thetas, K) +
      lgamma(thetas + 1.0) -
      lgamma(thetas + n)
    best_idx <- which.max(lls)
    return(list(theta = as.numeric(thetas[best_idx]), loglik = as.numeric(lls[best_idx]), status = "grid"))
  }

  idx <- sign_changes[1]
  bracket_lo <- as.numeric(thetas[idx])
  bracket_hi <- as.numeric(thetas[idx + 1L])
  lower <- min(bracket_lo, bracket_hi)
  upper <- max(bracket_lo, bracket_hi)

  theta_star <- tryCatch(
    {
      uniroot(
        f = function(t) .dlogl_dtheta_vectorized(alpha, t, n, K),
        interval = c(lower, upper),
        tol = 1e-12,
        maxiter = 250L
      )$root
    },
    error = function(e) {
      NA_real_
    }
  )

  if (!is.finite(theta_star)) {
    lls <- term_sizes +
      .logl_theta_term_vectorized(alpha, thetas, K) +
      lgamma(thetas + 1.0) -
      lgamma(thetas + n)
    best_idx <- which.max(lls)
    return(list(theta = as.numeric(thetas[best_idx]), loglik = as.numeric(lls[best_idx]), status = "grid_fallback"))
  }

  ll_star <- .compute_loglik(alpha, theta_star, n, K, term_sizes)
  list(theta = as.numeric(theta_star), loglik = as.numeric(ll_star), status = "brentq")
}

fit_pyp_profile_continuous <- function(
  sizes,
  theta_lo = 1e-10,
  theta_hi = 1e12,
  a_lo = 1e-6,
  a_hi = 1 - 1e-6,
  theta_grid = 64L,
  xatol = 1e-6
) {
  sizes <- as.integer(sizes)
  ctx <- .make_pyp_context(sizes)

  obj <- function(a) {
    out <- .find_optimal_theta(as.numeric(a), ctx, theta_lo, theta_hi, n_grid = theta_grid)
    -out$loglik
  }

  opt <- optimize(obj, interval = c(a_lo, a_hi), tol = xatol)
  a_hat <- as.numeric(opt$minimum)
  theta_out <- .find_optimal_theta(a_hat, ctx, theta_lo, theta_hi, n_grid = 256L)

  new_fit_result(
    alpha = a_hat,
    theta = theta_out$theta,
    loglik = theta_out$loglik,
    converged = TRUE,
    meta = list(status = theta_out$status, alpha_opt = opt)
  )
}

fit_pyp_profile <- function(
  sizes,
  alpha_grid = NULL,
  theta_lo = 1e-10,
  theta_hi = 1e12,
  refine = TRUE,
  verbose = FALSE
) {
  sizes_np <- as.integer(sizes)
  ctx <- .make_pyp_context(sizes_np)

  if (is.null(alpha_grid)) {
    alpha_grid <- c(
      seq(0.001, 0.05, length.out = 8),
      seq(0.05, 0.2, length.out = 12),
      seq(0.2, 0.8, length.out = 20),
      seq(0.8, 0.95, length.out = 10),
      seq(0.95, 0.999, length.out = 5)
    )
  }

  best_ll <- -Inf
  best_alpha <- 0.5
  best_theta <- 1.0
  best_status <- ""

  for (a in alpha_grid) {
    out <- .find_optimal_theta(a, ctx, theta_lo, theta_hi, n_grid = 128L)
    if (isTRUE(verbose)) {
      cat(sprintf("alpha=%.4f theta*=%.3g ll=%.6g (%s)\n", a, out$theta, out$loglik, out$status))
    }
    if (out$loglik > best_ll) {
      best_ll <- out$loglik
      best_alpha <- a
      best_theta <- out$theta
      best_status <- out$status
    }
  }

  if (isTRUE(refine)) {
    lo <- max(1e-8, best_alpha - 0.04)
    hi <- min(1.0 - 1e-8, best_alpha + 0.04)
    refine_grid <- seq(lo, hi, length.out = 250)

    for (a in refine_grid) {
      out <- .find_optimal_theta(a, ctx, theta_lo, theta_hi, n_grid = 40L)
      if (isTRUE(verbose)) {
        cat(sprintf("  [refine] alpha=%.4f theta*=%.3g ll=%.6g (%s)\n", a, out$theta, out$loglik, out$status))
      }
      if (out$loglik > best_ll) {
        best_ll <- out$loglik
        best_alpha <- a
        best_theta <- out$theta
        best_status <- out$status
      }
    }
  }

  new_fit_result(
    alpha = best_alpha,
    theta = best_theta,
    loglik = best_ll,
    converged = TRUE,
    meta = list(
      method = "profile_optimized",
      status = best_status,
      n = ctx$n,
      K = ctx$K
    )
  )
}

dlog_eppf_dtheta <- function(alpha, theta, n, K) {
  term_theta <- 0.0
  if (K > 1L) {
    u <- theta / alpha
    term_theta <- (digamma(u + K) - digamma(u + 1.0)) / alpha
  }
  term_norm <- digamma(theta + 1.0) - digamma(theta + n)
  as.numeric(term_theta + term_norm)
}

d2log_eppf_dtheta2 <- function(alpha, theta, n, K) {
  term_theta <- 0.0
  if (K > 1L) {
    u <- theta / alpha
    term_theta <- -(trigamma(u + 1.0) - trigamma(u + K)) / (alpha^2)
  }
  term_norm <- trigamma(theta + 1.0) - trigamma(theta + n)
  as.numeric(term_theta + term_norm)
}

profile_theta_for_alpha <- function(alpha, sizes, theta_lo = 1e-10, theta_hi = 1e12, n_grid = 80L) {
  sizes_np <- as.integer(sizes)
  ctx <- .make_pyp_context(sizes_np)
  .find_optimal_theta(alpha, ctx, theta_lo, theta_hi, n_grid = n_grid)
}
