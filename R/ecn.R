#' Fit the exchangeable correlated noise (ECN) model to correlated N(0,1) z-scores
#'
#' Estimates the realized empirical distribution of correlated, marginally N(0,1) z-scores
#' as \eqn{f(z) = \sum_{l=0}^L \omega_l \varphi^{(l)}(z)/\sqrt{l!}} by penalised maximum
#' likelihood subject to \eqn{f \ge 0}.
#'
#' @param z vector of z-scores.
#' @inheritParams cash
#' @return An object of class \code{"ecn"} with \code{omega}, \code{loglik}, \code{gd.order}, \code{z}.
#' @examples
#' z <- rnorm(1000)
#' fit <- ecn(z)
#' plot(fit)
#' @export
ecn <- function(z, gd.order = 10, omega.lambda = 10, omega.rho = 0.5, omega.pen = NULL,
                solver = c("newton", "mosek"), grid = seq(-10, 10, by = 0.001)) {
  solver <- match.arg(solver)
  L <- gd.order
  gam <- if (!is.null(omega.pen)) omega.pen else omega_penalty(L, omega.lambda, omega.rho)
  if (length(gam) != L) stop("omega.pen must have length gd.order")
  Mw <- std.gaussDeriv(z, L)
  if (solver == "newton") {
    B <- gauss_deriv_ratio(c(grid, z), L)[, -1, drop = FALSE]
    w <- c(1, solve_omega(Mw, B, gam, rep(0, L))$w)
  } else {
    legacy_check()
    w <- c(1, w_mosek(Mw, std.gaussDeriv(grid, L), gam)$w)
  }
  out <- list(gd.order = L, omega = w, loglik = sum(log(Mw %*% w)), z = z, solver = solver)
  class(out) <- "ecn"
  out
}

#' Method-of-moments estimate of the ECN coefficients
#'
#' \eqn{\hat\omega_l = \overline{(-1)^l He_l(z)/\sqrt{l!}}} (orthogonality of Hermite polynomials).
#' @param z vector of z-scores.
#' @param gd.order maximal order.
#' @export
ecn_moments <- function(z, gd.order = 10) colMeans(gauss_deriv_ratio(z, gd.order))

#' @export
print.ecn <- function(x, ...) {
  cat(sprintf("ecn fit (%s): n = %d, L = %d, loglik = %.4f\n", x$solver, length(x$z), x$gd.order, x$loglik))
  cat("omega:", format(round(x$omega, 4)), "\n")
  invisible(x)
}

#' @export
summary.ecn <- function(object, ...) object[c("gd.order", "omega", "loglik")]

#' Density of a fitted ECN model
#' @param fit an \code{"ecn"} or \code{"cash"} object.
#' @param z positions.
#' @export
density_ecn <- function(fit, z) as.vector(std.gaussDeriv(z, fit$gd.order) %*% fit$omega)

#' Histogram of z-scores with the fitted ECN density
#' @param x an \code{"ecn"} object.
#' @param breaks passed to [graphics::hist()].
#' @param std.norm also draw the N(0,1) density.
#' @param main title.
#' @param legend draw a legend.
#' @param ... further arguments to [graphics::hist()].
#' @export
plot.ecn <- function(x, breaks = 100, std.norm = TRUE, main = "", legend = TRUE, ...) {
  z <- x$z
  xp <- seq(-max(abs(z)) - 2, max(abs(z)) + 2, length = 1000)
  yp <- density_ecn(x, xp)
  h <- graphics::hist(z, breaks, plot = FALSE)
  ymax <- max(h$density, yp, if (std.norm) stats::dnorm(0))
  graphics::hist(z, breaks, prob = TRUE, ylim = c(0, ymax), main = main, xlab = "z", ...)
  graphics::lines(xp, yp, col = "red", lwd = 2)
  if (std.norm) graphics::lines(xp, stats::dnorm(xp), col = "blue", lwd = 2)
  if (legend) {
    graphics::legend("topright", lty = 1, lwd = 2, bty = "n",
                     col = c("red", if (std.norm) "blue"), c("Gaussian derivatives", if (std.norm) "N(0,1)"))
  }
  invisible(x)
}
