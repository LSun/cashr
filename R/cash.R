#' Penalty weights for the Gaussian-derivative coefficients
#'
#' \eqn{\gamma_l = 0} for odd \eqn{l} and \eqn{\lambda / \rho^{l/2}} for even \eqn{l}
#' (Sun & Stephens, Sec. 3).
#' @param L order of the Gaussian-derivative expansion.
#' @param lambda,rho penalty parameters; \code{lambda = NULL} means no penalty and
#'   \code{rho = NULL} a constant penalty \code{lambda} on the even orders.
#' @return numeric vector of length \code{L}.
#' @export
omega_penalty <- function(L, lambda = 10, rho = 0.5) {
  if (is.null(lambda)) return(rep(0, L))
  l <- seq_len(L)
  pen <- if (is.null(rho)) rep(lambda, L) else lambda / sqrt(rho^l)
  pen[l %% 2 == 1] <- 0
  pen
}

#' q-values from local false discovery rates
#' @param lfdr numeric vector.
#' @export
qval.from.lfdr <- function(lfdr) {
  q <- rep(NA_real_, length(lfdr))
  ok <- !is.na(lfdr)
  if (!any(ok)) return(q)
  o <- order(lfdr[ok])
  qq <- numeric(sum(ok))
  qq[o] <- cumsum(lfdr[ok][o]) / seq_along(o)
  q[ok] <- qq
  q
}

#' Empirical Bayes normal means with correlated noise
#'
#' Fits the model of Sun & Stephens (arXiv:1812.07488):
#' \deqn{X_j = \theta_j + s_j Z_j,\quad \theta_j \sim g = \sum_k \pi_k N(0, \sigma_k^2),\quad
#'       Z_j \sim f(z) = \sum_{l=0}^L \omega_l \varphi^{(l)}(z)/\sqrt{l!},}
#' where the \eqn{Z_j} are marginally \eqn{N(0,1)} but may be correlated. The prior \eqn{g}
#' and the realized noise distribution \eqn{f} are estimated jointly by maximising the
#' penalised likelihood (eq. 18) by exact block-coordinate ascent.
#'
#' @param x vector of observations (e.g. estimated effects).
#' @param s standard errors: a positive scalar or a vector the length of \code{x}.
#' @param deltaAt0 include a point mass at 0 in the prior.
#' @param gd.order order \eqn{L} of the Gaussian-derivative expansion.
#' @param omega.lambda,omega.rho,omega.pen penalty on \eqn{\omega}; see [omega_penalty()].
#'   \code{omega.pen} (length \code{gd.order}) overrides the other two.
#' @param mixsd.mult ratio between consecutive prior standard deviations.
#' @param null.weight Dirichlet prior weight on the null component (\eqn{\lambda_0 =} \code{null.weight - 1}).
#' @param gd.priority update the noise model \eqn{\omega} before the prior \eqn{\pi} in each
#'   iteration (default). The fit starts from "everything is null" (\eqn{\pi} = point mass at 0,
#'   \eqn{N(0,1)} noise), so with \code{gd.priority = TRUE} the noise is fitted first and signal
#'   enters only where the noise model cannot explain the data. This is the conservative choice:
#'   on all-null correlated data the \eqn{\pi}-first order can stop at a local optimum that calls
#'   most observations signal. \code{FALSE} updates \eqn{\pi} first, the cashr 0.2.0 default.
#' @param solver \code{"newton"} (default; fast, no dependencies) or \code{"mosek"} (the 2018
#'   algorithm: SQUAREM with REBayes and MOSEK; needs Rmosek and a licence).
#' @param tol relative tolerance on the objective for the newton solver.
#' @param maxiter maximal number of outer iterations (\code{"mosek"}: SQUAREM fixed-point
#'   evaluations; the 2018 default was 50).
#' @param grid points where \eqn{f(z) \ge 0} is enforced.
#' @return An object of class \code{"cash"}: \code{fitted_g} (data frame \code{pi}, \code{mean},
#'   \code{sd}), \code{omega}, \code{pi0}, \code{penloglik}, \code{loglik}, \code{lfdr},
#'   \code{qvalue}, \code{pm} (posterior mean), \code{lfsr}, \code{svalue}, \code{niter},
#'   \code{converged}.
#' @examples
#' set.seed(1)
#' theta <- c(rep(0, 800), rnorm(200, 0, 2))
#' x <- theta + rnorm(1000)
#' fit <- cash(x, 1)
#' fit
#' sum(fit$qvalue <= 0.1)
#' @export
cash <- function(x, s = 1, deltaAt0 = TRUE, gd.order = 10,
                 omega.lambda = 10, omega.rho = 0.5, omega.pen = NULL,
                 mixsd.mult = sqrt(2), null.weight = 10, gd.priority = TRUE,
                 solver = c("newton", "mosek"), tol = 1e-8, maxiter = NULL,
                 grid = seq(-10, 10, by = 0.001)) {
  solver <- match.arg(solver)
  t_start <- proc.time()[["elapsed"]]
  x <- as.numeric(x)
  if (length(s) == 1L) s <- rep(s, length(x))
  if (length(s) != length(x) || any(!is.finite(s)) || any(s <= 0) || any(!is.finite(x))) {
    stop("x must be finite and s positive (a scalar or a vector the length of x)")
  }
  L <- gd.order
  sd <- autoselect_mixsd(x, s, mixsd.mult)
  if (deltaAt0) sd <- c(0, sd)
  K <- length(sd)
  lam <- rep(0, K)
  if (deltaAt0) lam[1] <- null.weight - 1
  gam <- if (!is.null(omega.pen)) omega.pen else omega_penalty(L, omega.lambda, omega.rho)
  if (length(gam) != L) stop("omega.pen must have length gd.order")

  arr <- cash_arrays(x, s, sd, L)
  F <- arr$F

  if (solver == "newton") {
    if (is.null(maxiter)) maxiter <- 500
    B <- gauss_deriv_ratio(grid, L)[, -1, drop = FALSE]
    ## start from "everything is null": pi = point mass at 0 (if there is one), N(0, 1) noise
    pi <- if (deltaAt0) c(1, rep(0, K - 1)) else rep(1 / K, K); omega <- rep(0, L)
    objective <- function(pi, omega) {
      dens <- contract_KL(F, pi, c(1, omega))
      if (any(dens <= 0)) return(-Inf)
      sum(log(dens)) + sum(lam[lam > 0] * log(pi[lam > 0])) - sum(gam * abs(omega))
    }
    ## Each block is solved exactly, so the objective can only increase; a sub-solve that
    ## returns a worse point (numerical failure) is rejected.
    pi_update <- function(pi, omega) {
      new <- solve_pi(contract_L(F, c(1, omega)), lam, pi)$pi
      if (objective(new, omega) >= objective(pi, omega)) new else pi
    }
    omega_update <- function(pi, omega) {
      new <- solve_omega(contract_K(F, pi), B, gam, omega)$w
      if (objective(pi, new) >= objective(pi, omega)) new else omega
    }
    obj_old <- -Inf; converged <- FALSE
    for (it in seq_len(maxiter)) {
      if (gd.priority) {
        omega <- omega_update(pi, omega); pi <- pi_update(pi, omega)
      } else {
        pi <- pi_update(pi, omega); omega <- omega_update(pi, omega)
      }
      obj <- objective(pi, omega)
      if (obj - obj_old <= tol * (1 + abs(obj))) { converged <- TRUE; break }
      obj_old <- obj
    }
    w_full <- c(1, omega)
  } else {
    legacy_check()
    if (is.null(maxiter)) maxiter <- 50
    scale <- exp(arr$log_scale)
    Fu <- aperm(F * scale, c(2, 3, 1))                       # (K, L+1, n), unscaled as in 2018
    matrix_lik_z <- std.gaussDeriv(grid, L)
    res <- legacy_biopt(Fu, matrix_lik_z, lam + 1, gam, maxiter = maxiter, gd.priority = gd.priority)
    pi <- res$pi; w_full <- res$omega; it <- res$niter; converged <- res$converged
    obj <- sum(log(contract_KL(F, pi, w_full))) + sum(lam[lam > 0] * log(pi[lam > 0])) -
      sum(gam * abs(w_full[-1]))
  }

  num_k <- contract_L(F, w_full)
  dens <- as.vector(num_k %*% pi)
  lfdr <- if (deltaAt0) pmin(pmax(pi[1] * num_k[, 1] / dens, 0), 1) else rep(0, length(x))
  pm <- contract_KL(posterior_mean_array(arr), pi, w_full) / dens
  pos <- pmin(pmax(contract_KL(positive_prob_array(arr), pi, w_full) / dens, 0), 1 - lfdr)
  neg <- 1 - lfdr - pos
  lfsr <- pmin(pos, neg) + lfdr
  out <- list(fitted_g = data.frame(pi = pi, mean = 0, sd = sd), gd.order = L, omega = w_full,
              pi0 = if (deltaAt0) pi[1] else 0,
              penloglik = obj + sum(arr$log_scale), loglik = sum(log(dens)) + sum(arr$log_scale),
              niter = it, converged = converged, solver = solver,
              lfdr = lfdr, qvalue = qval.from.lfdr(lfdr), pm = pm, lfsr = lfsr,
              svalue = qval.from.lfdr(lfsr), x = x, s = s, deltaAt0 = deltaAt0,
              seconds = proc.time()[["elapsed"]] - t_start)
  class(out) <- "cash"
  out
}

#' @export
print.cash <- function(x, ...) {
  cat(sprintf("cash fit (%s): n = %d, K = %d, L = %d, pi0 = %.4f, penloglik = %.4f, %d iterations%s, %.2fs\n",
              x$solver, length(x$x), nrow(x$fitted_g), x$gd.order, x$pi0, x$penloglik, x$niter,
              if (x$converged) "" else " (not converged)", x$seconds))
  cat("omega:", format(round(x$omega, 4)), "\n")
  invisible(x)
}

#' @export
summary.cash <- function(object, ...) {
  object[c("fitted_g", "gd.order", "omega", "pi0", "penloglik", "niter", "converged")]
}

#' Local false sign rates and s-values of a cash fit (kept for 0.2.0 compatibility)
#' @param output a \code{"cash"} object.
#' @export
get_svalue <- function(output) list(lfsr = output$lfsr, svalue = output$svalue)

#' CDF of the fitted prior g
#' @param fit a \code{"cash"} object.
#' @param q numeric vector.
#' @export
cdf_g <- function(fit, q) {
  g <- fit$fitted_g
  vapply(q, function(t) sum(g$pi * ifelse(g$sd > 0, stats::pnorm(t, 0, pmax(g$sd, 1e-300)), as.numeric(t >= 0))), 0)
}
