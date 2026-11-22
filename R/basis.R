# Gaussian-derivative bases and the likelihood / posterior arrays of the cash model.
#
#   X_j = theta_j + s_j Z_j,  theta_j ~ g = sum_k pi_k N(0, sigma_k^2),
#   Z_j ~ f(z; omega) = sum_{l=0}^L omega_l phi^(l)(z) / sqrt(l!),  omega_0 = 1.
#
# Hermite quantities use the orthonormal probabilists' recursion
#   h_0 = 1, h_1 = y, h_{l+1} = (y h_l - sqrt(l) h_{l-1}) / sqrt(l+1),
# so h_l = He_l / sqrt(l!) and phi^(l)(y) / sqrt(l!) = (-1)^l h_l(y) phi(y).

#' Orthonormal probabilists' Hermite polynomials
#'
#' @param y numeric vector.
#' @param L maximal order.
#' @return \code{length(y)} x \code{(L+1)} matrix with columns \eqn{He_l(y)/\sqrt{l!}}, \eqn{l = 0..L}.
#' @export
hermite_orthonormal <- function(y, L) {
  H <- matrix(0, length(y), L + 1)
  H[, 1] <- 1
  if (L >= 1) H[, 2] <- y
  if (L >= 2) for (l in 1:(L - 1)) H[, l + 2] <- (y * H[, l + 1] - sqrt(l) * H[, l]) / sqrt(l + 1)
  H
}

#' Standardized Gaussian derivatives
#'
#' Values of \eqn{\varphi^{(l)}(x)/\sqrt{l!}}, \eqn{l = 0, \dots, L}: the basis of the
#' exchangeable correlated noise (ECN) model.
#'
#' @param x numeric vector of positions.
#' @param gd.order maximal order \eqn{L}.
#' @return \code{length(x)} x \code{(gd.order+1)} matrix.
#' @examples
#' x <- seq(-5, 5, length = 200)
#' matplot(x, std.gaussDeriv(x, 4), type = "l")
#' @export
std.gaussDeriv <- function(x, gd.order) {
  gauss_deriv_ratio(x, gd.order) * stats::dnorm(x)
}

# phi^(l)(z) / (sqrt(l!) phi(z)) = (-1)^l h_l(z)
gauss_deriv_ratio <- function(z, L) {
  sweep(hermite_orthonormal(z, L), 2, (-1)^(0:L), "*")
}

autoselect_mixsd <- function(x, s, mult = sqrt(2)) {
  s <- s[s != 0]
  smin <- min(s) / 10
  smax <- if (all(x^2 <= s^2)) 8 * smin else 2 * sqrt(max(x^2 - s^2))
  if (mult == 0) return(c(0, smax / 2))
  npoint <- ceiling(log2(smax / smin) / log2(mult))
  mult^((-npoint):0) * smax
}

# Hermite array for an n x K matrix y: list of L+1 matrices (n x K)
hermite_list <- function(y, L) {
  H <- vector("list", L + 1)
  H[[1]] <- matrix(1, nrow(y), ncol(y))
  if (L >= 1) H[[2]] <- y
  if (L >= 2) for (l in 1:(L - 1)) H[[l + 2]] <- (y * H[[l + 1]] - sqrt(l) * H[[l]]) / sqrt(l + 1)
  H
}

# Likelihood array F[j, k, l+1] = p_jkl = s_j^l / tau^(l+1) phi^(l)(x/tau) / sqrt(l!), row-scaled:
# row j is divided by exp(log_scale[j]) = max_k phi(x_j / tau_jk) / tau_jk.
cash_arrays <- function(x, s, sd, L) {
  n <- length(x); K <- length(sd)
  tau <- sqrt(outer(s^2, sd^2, "+"))
  y <- x / tau
  logg0 <- -0.5 * y^2 - 0.5 * log(2 * pi) - log(tau)
  log_scale <- apply(logg0, 1, max)
  g0 <- exp(logg0 - log_scale)
  r <- s / tau
  H <- hermite_list(y, L + 1)
  F <- array(0, c(n, K, L + 1))
  rl <- matrix(1, n, K)
  for (l in 0:L) {
    F[, , l + 1] <- g0 * rl * (-1)^l * H[[l + 1]]
    rl <- rl * r
  }
  list(x = x, s = s, sd = sd, L = L, tau = tau, y = y, g0 = g0, r = r, H = H,
       log_scale = log_scale, F = F)
}

# Numerators of E[theta | x]: PM = -sigma^2 s^l tau^-(l+2) phi^(l+1)(y) / sqrt(l!)
posterior_mean_array <- function(a) {
  n <- length(a$x); K <- length(a$sd); L <- a$L
  sig2_tau <- matrix(a$sd^2, n, K, byrow = TRUE) / a$tau
  out <- array(0, c(n, K, L + 1))
  rl <- matrix(1, n, K)
  for (l in 0:L) {
    out[, , l + 1] <- a$g0 * rl * sig2_tau * (-1)^l * sqrt(l + 1) * a$H[[l + 2]]
    rl <- rl * a$r
  }
  out
}

# Numerators of P(theta > 0 | x): s^l / sqrt(l!) d^l/dx^l [N(x; 0, tau^2) Phi(c x)], c = sigma/(s tau),
# expanded with Leibniz' rule. Point-mass components contribute 0.
positive_prob_array <- function(a) {
  n <- length(a$x); K <- length(a$sd); L <- a$L
  out <- array(0, c(n, K, L + 1))
  pos <- which(a$sd > 0)
  if (length(pos) == 0) return(out)
  sub <- function(M) M[, pos, drop = FALSE]
  y <- sub(a$y); r <- sub(a$r); g0 <- sub(a$g0)
  q <- matrix(a$sd[pos], n, length(pos), byrow = TRUE) / sub(a$tau)
  cx <- y * outer(1 / a$s, a$sd[pos])
  Hy <- lapply(a$H[1:(L + 1)], sub)
  Hc <- hermite_list(cx, max(L - 1, 0))
  Phi_c <- stats::pnorm(cx)
  phi_c <- stats::dnorm(cx)
  for (l in 0:L) {
    term <- r^l * (-1)^l * Hy[[l + 1]] * Phi_c
    if (l >= 1) for (m in 1:l) {
      coef <- choose(l, m) * exp(0.5 * (lgamma(l - m + 1) + lgamma(m) - lgamma(l + 1)))
      term <- term + coef * (-1)^(l - 1) * r^(l - m) * q^m * Hy[[l - m + 1]] * Hc[[m]] * phi_c
    }
    out[, pos, l + 1] <- g0 * term
  }
  out
}

# contract a (n, K, L+1) array with weights over K and/or L
contract_KL <- function(A, pi, w) {          # -> n
  d <- dim(A)
  M <- matrix(A, d[1] * d[2], d[3]) %*% w     # (n*K)
  as.vector(matrix(M, d[1], d[2]) %*% pi)
}
contract_L <- function(A, w) {               # -> n x K
  d <- dim(A)
  matrix(matrix(A, d[1] * d[2], d[3]) %*% w, d[1], d[2])
}
contract_K <- function(A, pi) {              # -> n x (L+1)
  d <- dim(A)
  out <- matrix(0, d[1], d[3])
  for (l in seq_len(d[3])) out[, l] <- A[, , l] %*% pi
  out
}
