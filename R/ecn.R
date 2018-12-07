#' @title Fit empirical distribution of observed correlated N(0,1) noise using Exchangeable Correlated Noise (ECN) model
#'
#' @export

ecn = function (z, gd.order = 10, omega.lambda = 10, omega.rho = 0.5, omega.pen = NULL) {
  L <- gd.order

  if (is.null(omega.pen)) {
    if (is.null(omega.lambda)) {
      omega_prior = rep(0, L)
    } else if (is.null(omega.rho)) {
      omega_prior = rep(omega.lambda, L)
      omega_prior[seq(1, L, by = 2)] = 0
    } else {
      omega_prior = omega.lambda / sqrt(omega.rho^(1 : L))
      omega_prior[seq(1, L, by = 2)] = 0
    }
  } else if (length(omega_pen) == length(L)) {
    omega_prior = omega.pen
  } else {
    stop('The length of penalty parameters of omega should be L.')
  }

  hermite = Hermite(L)
  gd0.std = dnorm(z)
  matrix_lik_w = cbind(gd0.std)
  for (i in 1 : L) {
    gd.std = (-1)^i * hermite[[i]](z) * gd0.std / sqrt(factorial(i))
    matrix_lik_w = cbind(matrix_lik_w, gd.std)
  }

  z.extra = seq(-10, 10, by = 0.001)
  gd0.std = dnorm(z.extra)
  matrix_lik_z = cbind(gd0.std)
  for (i in 1 : L) {
    gd.std = (-1)^i * hermite[[i]](z.extra) * gd0.std / sqrt(factorial(i))
    matrix_lik_z = cbind(matrix_lik_z, gd.std)
  }

  w.fit = w.mosek(matrix_lik_w, matrix_lik_z, omega_prior, w.init = NULL)
  w.hat = c(1, w.fit$w)
  w.status = w.fit$status
  loglik.hat = sum(log(matrix_lik_w %*% w.hat))

  output <- list(gd.order = L, omega = w.hat, loglik = loglik.hat, status = w.status, z = z)
  class(output) <- 'ecn'

  return(output)
}

#' @export

summary.ecn <- function (output, ...) {
  output[1 : 4]
}

#' @export

print.ecn <- function (output, ...) {
  print(summary.ecn(output))
}

gdfit.mom = function (z, gd.ord) {
  hermite.list = orthopolynom::hermite.he.polynomials(gd.ord)
  hermite.coef = orthopolynom::polynomial.coefficients(hermite.list)
  moments = c()
  for (i in 0 : gd.ord) {
    moments[i + 1] = mean(z^i)
  }
  w = c()
  for (i in 0 : gd.ord) {
    w[i + 1] = sum(moments[1 : (i + 1)] / sqrt(factorial(i)) * hermite.coef[[i + 1]]) * (-1)^i
  }
  return(list(gd.ord = gd.ord, w = w))
}

#' @export

plot.ecn = function (output, symm = TRUE, breaks = 100, std.norm = TRUE, main, legend = TRUE) {
  z <- output$z
  w <- output$omega
  gd.ord <- output$gd.order

  if (symm) {
    x.plot = seq(- max(abs(z)) - 2, max(abs(z)) + 2, length = 1000)
  } else {
    x.plot = seq(min(z) - 2, max(z) + 2, length = 1000)
  }

  hermite = Hermite(gd.ord)
  gd0.std = dnorm(x.plot)
  matrix_lik_plot = cbind(gd0.std)
  for (i in 1 : gd.ord) {
    gd.std = (-1)^i * hermite[[i]](x.plot) * gd0.std / sqrt(factorial(i))
    matrix_lik_plot = cbind(matrix_lik_plot, gd.std)
  }
  y.plot = matrix_lik_plot %*% w
  z.hist = hist(z, breaks, plot = FALSE)
  if (std.norm) {
    y.max = max(z.hist$density, y.plot, dnorm(0))
  } else {
    y.max = max(z.hist$density, y.plot)
  }
  if (!missing(main)) {
    hist(z, breaks, prob = TRUE, ylim = c(0, y.max), main = main)
  } else {
    hist(z, breaks, prob = TRUE, ylim = c(0, y.max))
  }
  lines(x.plot, y.plot, col = "blue")
  if (legend) {
    legend("topright", lty = 1, col = "blue", "Gaussian Derivatives")
  }
  if (std.norm) {
    lines(x.plot, dnorm(x.plot), col = "red")
    if (legend) {
      legend("topleft", lty = 1, col = "red", "N(0, 1)")
    }
  }
}


w.mosek = function (matrix_lik_w, matrix_lik_z, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[, 1]
  B = matrix_lik_z[, -1]
  b = matrix_lik_z[, 1]
  m = ncol(A)
  nA = nrow(A)
  nB = nrow(B)
  AB <- rbind(A, B)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    g.init <- as.vector(matrix_lik_w %*% w.init)
    v.init <- c(1 / g.init, rep(0, nB))
    v.init.list <- list(xx = v.init)
    P$sol <- list(itr = v.init.list, bas = v.init.list)
  }
  P$c <- c(a, b)
  P$A <- Matrix::Matrix(t(AB), sparse = TRUE)
  if (is.null(w_prior) | all(w_prior == 0) | missing(w_prior)) {
    P$bc <- rbind(rep(0, m), rep(0, m))
  } else {
    P$bc <- rbind(-w_prior, w_prior)
  }
  P$bx <- rbind(rep(0, nA + nB), rep(Inf, nA + nB))
  opro <- matrix(list(), nrow = 5, ncol = nA)
  rownames(opro) <- c("type", "j", "f", "g", "h")
  opro[1, ] <- as.list(rep("log", nA))
  opro[2, ] <- as.list(1 : nA)
  opro[3, ] <- as.list(rep(-1, nA))
  opro[4, ] <- as.list(rep(1, nA))
  opro[5, ] <- as.list(rep(0, nA))
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$suc - z$sol$itr$slc
  list(w = w, status = status)
}
