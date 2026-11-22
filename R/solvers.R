# Convex sub-problem solvers (port of the Python package's solvers.py).
#
# Both blocks of the biconvex fit have the form
#   minimise  c'x - sum_i w_i log(M_i x + m_i)   subject to  G x + g > 0
# with few variables (K ~ 20 or L ~ 10) but many rows (n observations, ~20k grid
# constraints): ideal for a log-barrier Newton method, where each step costs one pass
# over the rows plus a small Cholesky factorisation.

log_barrier_newton <- function(x0, cvec, M, m, w, G = NULL, g = NULL, tol = 1e-12, mu = 50,
                               t0 = 1, max_newton = 200, newton_tol = 1e-10) {
  # M may have fewer columns than x has entries: it then acts on the leading entries only.
  # Residuals r = M x + m and q = G x + g are updated along the search direction, so each
  # line-search trial is cheap; an accepted point is re-checked with r and q recomputed exactly.
  x <- x0
  d <- length(x)
  p <- ncol(M); ip <- seq_len(p)
  n_ineq <- if (is.null(G)) 0L else nrow(G)
  bval <- function(r, q, x, t) {
    if (any(r <= 0) || (n_ineq && any(q <= 0))) return(Inf)
    t * (sum(cvec * x) - sum(w * log(r))) - (if (n_ineq) sum(log(q)) else 0)
  }
  t <- if (n_ineq) t0 else 1
  r <- as.vector(M %*% x[ip]) + m
  gap_target <- tol * (1 + abs(sum(cvec * x) - sum(w * log(r))))
  steps <- 0L
  repeat {
    r <- as.vector(M %*% x[ip]) + m
    q <- if (n_ineq) as.vector(G %*% x) + g else NULL
    f0 <- bval(r, q, x, t)
    for (it in seq_len(max_newton)) {
      wr <- w / r
      grad <- t * cvec
      grad[ip] <- grad[ip] - t * as.vector(crossprod(M, wr))
      H <- matrix(0, d, d)
      H[ip, ip] <- t * crossprod(M * sqrt(wr / r))
      if (n_ineq) {
        iq <- 1 / q
        grad <- grad - as.vector(crossprod(G, iq))
        H <- H + crossprod(G * iq)
      }
      diag(H) <- diag(H) + 1e-14 * max(1, max(abs(diag(H))))
      dx <- tryCatch({R <- chol(H); -backsolve(R, forwardsolve(t(R), grad))},
                     error = function(e) -qr.solve(H, grad))
      dec2 <- -sum(grad * dx)
      steps <- steps + 1L
      if (dec2 <= 2 * newton_tol || dec2 <= 1e-13 * (1 + abs(f0))) break
      Mdx <- as.vector(M %*% dx[ip])
      Gdx <- if (n_ineq) as.vector(G %*% dx) else NULL
      step <- 1; slope <- sum(grad * dx); ok <- FALSE
      while (step > 1e-10) {
        r1 <- r + step * Mdx
        q1 <- if (n_ineq) q + step * Gdx else NULL
        f1 <- bval(r1, q1, x + step * dx, t)
        if (f1 <= f0 + 0.25 * step * slope) {
          ## accept only if the point is strictly feasible when r and q are recomputed exactly:
          ## a slack of ~1e-16 can be positive when updated along dx yet round to 0 at x
          xn <- x + step * dx
          r1 <- as.vector(M %*% xn[ip]) + m
          q1 <- if (n_ineq) as.vector(G %*% xn) + g else NULL
          f1 <- bval(r1, q1, xn, t)
          if (f1 <= f0 + 0.25 * step * slope) { ok <- TRUE; break }
        }
        step <- step / 2
      }
      if (!ok) break
      x <- x + step * dx; r <- r1; q <- q1
      f0 <- f1
    }
    if (n_ineq == 0 || n_ineq / t < gap_target) break
    t <- t * mu
  }
  list(x = x, steps = steps)
}

# pi-step: max sum_j log(Mpi pi)_j + sum_k lam_k log pi_k over the simplex.
# Minimising (n + sum lam) sum(pi) - sum log(Mpi pi) - sum lam log(pi) over pi >= 0 has its
# optimum exactly on sum(pi) = 1 (Koenker-Mizera), so no equality constraint is needed.
solve_pi <- function(Mpi, lam, pi0, tol = 1e-12) {
  n <- nrow(Mpi); K <- ncol(Mpi)
  pen <- lam > 0
  M <- rbind(Mpi, diag(K)[pen, , drop = FALSE])
  w <- c(rep(1, n), lam[pen])
  x0 <- pmax(pi0, 1e-3 / K); x0 <- x0 / sum(x0)
  res <- log_barrier_newton(x0, rep(sum(w), K), M, rep(0, nrow(M)), w, diag(K), rep(0, K),
                            tol = tol, t0 = 1e2)
  pi <- pmax(res$x, 0)
  list(pi = pi / sum(pi), steps = res$steps)
}

local_minima <- function(v) {
  if (diff(range(v)) <= 1e-12) return(integer(0))
  m <- length(v)
  left <- c(TRUE, v[-1] <= v[-m])
  right <- c(v[-m] <= v[-1], TRUE)
  which(left & right)
}

add_windows <- function(active, centers, half) {
  m <- length(active)
  for (i in centers) active[max(1, i - half):min(m, i + half)] <- TRUE
  active
}

# omega-step: max sum_j log(a + A w) - sum_l gamma_l |w_l|  s.t. 1 + B w >= 0 on the grid.
# Mw: n x (L+1) (column 1 is the omega_0 = 1 part). B: grid rows of f(z)/phi(z), l = 1..L,
# ordered along the grid. Cutting planes: f/phi is a degree-L polynomial with few local minima,
# so solve on a coarse subset plus windows around the minima, check the full grid and add any
# violated region. A subset solution feasible on the full grid is optimal for the full problem.
# `state` (an environment) keeps the working set across calls.
solve_omega <- function(Mw, B, gamma, w0, tol = 1e-12, coarse = 50, window = 25, near = 0.25,
                        max_rounds = 20) {
  ## Each row enters as log(a_i + A_i w), so scaling a row by a positive constant shifts the
  ## objective by a constant and leaves the optimum unchanged. Rows are scaled to max |entry| 1:
  ## an observation far in the tail of the current fit (null density ~1e-170 at the all-null
  ## start) would otherwise overflow the barrier Hessian through 1 / r^2.
  Mw <- Mw / apply(abs(Mw), 1, max)
  a <- Mw[, 1]; A <- Mw[, -1, drop = FALSE]
  n <- nrow(A); L <- ncol(A); m_grid <- nrow(B)
  P <- which(gamma > 0); nP <- length(P)
  cvec <- c(rep(0, L), gamma[P])
  E <- matrix(0, nP, L); E[cbind(seq_len(nP), P)] <- 1
  G_pen <- rbind(cbind(-E, diag(nP)), cbind(E, diag(nP)))
  feasible <- function(wv) all(a + A %*% wv > 0) && all(1 + B %*% wv > 0)
  w <- if (feasible(0.99 * w0)) 0.99 * w0 else rep(0, L)
  ## The working set is rebuilt on every call: carrying it across outer iterations piles up
  ## near-duplicate tail rows, the barrier Hessian becomes nearly singular and Newton stalls.
  active <- rep(FALSE, m_grid); active[seq(1, m_grid, by = coarse)] <- TRUE; active[m_grid] <- TRUE
  cur <- as.vector(1 + B %*% w)
  mins <- local_minima(cur)
  active <- add_windows(active, mins[cur[mins] < near], window)
  steps <- 0L
  for (round in seq_len(max_rounds)) {
    idx <- which(active)
    G <- rbind(cbind(B[idx, , drop = FALSE], matrix(0, length(idx), nP)), G_pen)
    g <- c(rep(1, length(idx)), rep(0, 2 * nP))
    x0 <- c(w, abs(w[P]) + 1e-2 * (1 + abs(w[P])))
    res <- log_barrier_newton(x0, cvec, A, a, rep(1, n), G, g, tol = tol, t0 = 1)
    steps <- steps + res$steps
    w_new <- res$x[1:L]
    full <- as.vector(1 + B %*% w_new)
    mins <- local_minima(full)
    active <- add_windows(active, mins[full[mins] < near], window)
    if (all(full > 0)) return(list(w = w_new, steps = steps))
    active[which.min(full)] <- TRUE
    full_old <- as.vector(1 + B %*% w)
    lo <- 0; hi <- 1
    for (k in 1:60) {
      mid <- (lo + hi) / 2
      if (all(full_old + mid * (full - full_old) > 0) && all(a + A %*% (w + mid * (w_new - w)) > 0)) lo <- mid else hi <- mid
    }
    w <- w + 0.99 * lo * (w_new - w)
  }
  stop("omega cutting-plane loop did not converge")
}
