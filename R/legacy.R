# The 2018 (v0.2.0) fitting algorithm, kept for reproducing the paper's original code path:
# SQUAREM over the two block updates, with the pi-step solved by REBayes::KWDual and the
# omega-step by MOSEK. The 2018 code used Rmosek's `scopt` interface, which Rmosek >= 9
# removed; `w_mosek` solves the same problem with exponential cones (Rmosek >= 10).
# Used by cash(..., solver = "mosek") and ecn(..., solver = "mosek"). Needs Rmosek, REBayes,
# SQUAREM and a MOSEK licence.

legacy_check <- function() {
  for (p in c("Rmosek", "REBayes", "SQUAREM")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("solver = \"mosek\" needs package ", p, " (and a MOSEK licence); use solver = \"newton\"",
           if (p == "Rmosek") paste0(".\nTo install Rmosek, see https://docs.mosek.com/12.0/rmosek/install-interface.html") else "")
    }
  }
}

# max sum_j log(a_j + A_j w) - sum_l gamma_l |w_l|  s.t.  b + B w >= 0   (via MOSEK, exp cones)
w_mosek <- function(matrix_lik_w, matrix_lik_z, w_prior) {
  A <- matrix_lik_w[, -1, drop = FALSE]; a <- matrix_lik_w[, 1]
  B <- matrix_lik_z[, -1, drop = FALSE]; b <- matrix_lik_z[, 1]
  L <- ncol(A); n <- nrow(A); nB <- nrow(B)
  if (is.null(w_prior)) w_prior <- rep(0, L)
  Bs <- B / b                                  # same feasible set, better conditioned
  nv <- 2 * L + n                              # variables (w, u, t)
  I <- diag(L)
  Alin <- rbind(cbind(Bs, matrix(0, nB, L)), cbind(-I, I), cbind(I, I))
  P <- list(sense = "max")
  P$c <- c(rep(0, L), -w_prior, rep(1, n))
  P$A <- Matrix::Matrix(cbind(Alin, matrix(0, nrow(Alin), n)), sparse = TRUE)
  P$bc <- rbind(c(rep(-1, nB), rep(0, 2 * L)), rep(Inf, nB + 2 * L))
  P$bx <- rbind(rep(-Inf, nv), rep(Inf, nv))
  rowsA <- rep(3 * (0:(n - 1)) + 1, times = L)
  P$F <- Matrix::sparseMatrix(i = c(rowsA, 3 * (1:n)), j = c(rep(1:L, each = n), 2 * L + (1:n)),
                              x = c(as.vector(A), rep(1, n)), dims = c(3 * n, nv))
  P$g <- as.vector(rbind(a, 1, 0))
  P$cones <- matrix(list("PEXP", 3, NULL), nrow = 3, ncol = n)
  rownames(P$cones) <- c("type", "dim", "conepar")
  z <- Rmosek::mosek(P, opts = list(verbose = 0))
  list(w = z$sol$itr$xx[1:L], status = z$sol$itr$solsta)
}

# pi-step via the Koenker-Mizera dual (REBayes::KWDual), as in ashr's mixIP
mix_ip <- function(matrix_lik, prior, rtol = 1e-6) {
  n <- nrow(matrix_lik); k <- ncol(matrix_lik)
  A <- rbind(diag(length(prior)), matrix_lik)
  w <- c(prior - 1, rep(1, n))
  A <- A[w != 0, , drop = FALSE]; w <- w[w != 0]
  res <- REBayes::KWDual(A, rep(1, k), w / sum(w), control = list(rtol = rtol))
  res$f / sum(res$f)
}

# array_F: (K, L+1, n) unscaled likelihood array, as in the 2018 code
legacy_biopt <- function(array_F, matrix_lik_z, pi_prior, w_prior, maxiter = 50, tol = 1e-7,
                         gd.priority = FALSE) {
  K <- dim(array_F)[1]; Lw <- dim(array_F)[2]
  penloglik <- function(pinw) {
    pi <- pinw[1:K]; w <- pinw[-(1:K)]
    dens <- colSums(t(apply(pi * array_F, 2, colSums)) * w)
    sub <- pi_prior != 1
    sum(log(pmax(0, dens))) + sum((pi_prior - 1)[sub] * log(pi[sub])) - sum(abs(w[-1]) * w_prior)
  }
  ## gd.priority as in cashr 0.2.0's bifixpoint(): omega step first, then pi
  pi_step <- function(w) mix_ip(apply(aperm(array_F, c(2, 1, 3)) * w, 2, colSums), pi_prior)
  w_step <- function(pi) c(1, w_mosek(apply(array_F * pi, 2, colSums), matrix_lik_z, w_prior)$w)
  fixpt <- function(pinw) {
    pi <- pinw[1:K]; w <- pinw[-(1:K)]
    if (gd.priority) { w_new <- w_step(pi); c(pi_step(w_new), w_new) }
    else { pi_new <- pi_step(w); c(pi_new, w_step(pi_new)) }
  }
  ctrl <- list(K = 1, method = 3, square = TRUE, step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1,
               objfn.inc = 1, tol = tol, maxiter = maxiter, trace = FALSE)
  res <- SQUAREM::squarem(par = c(1, rep(0, K - 1), 1, rep(0, Lw - 1)), fixptfn = fixpt,
                          objfn = function(p) -penloglik(p), control = ctrl)
  pi <- pmax(0, res$par[1:K])
  list(pi = pi / sum(pi), omega = res$par[-(1:K)], niter = res$fpevals, converged = res$convergence)
}
