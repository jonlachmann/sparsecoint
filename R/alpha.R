#' Function to estimate alpha
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param zbeta estimate of short-run effects
#' @param rank cointegration rank
#' @param Omega estimate of inverse error covariance matrix
#' @param P transformation matrix P derived from Omega
#' @param beta estimate of cointegrating vector
#' @param intercept F do not include intercept, T include intercept in estimation short-run effects
#' @return A list containing:
#' alpha: estimate of adjustment coefficients
#' alphastar: estimate of transformed adjustment coefficients
nts.alpha.procrusted <- function(Y, X, Z, zbeta, rank, Omega, P, beta, intercept = F, exo = NULL) {
  # Data matrices
  if (!is.null(exo)) X <- cbind(exo, X)
  if (intercept) X <- cbind(1, X)
  Ymatrix <- Y - X %*% zbeta

  Xmatrix <- Z %*% beta
  decomp <- eigen(Omega)
  Pmin <- (decomp$vectors) %*% diag(1 / sqrt(decomp$values)) %*% solve(decomp$vectors)

  # Singular Value Decomposition to compute ALPHA
  SingDecomp <- svd(t(Xmatrix) %*% Ymatrix %*% P)
  alpha_star <- t(SingDecomp$u %*% t(SingDecomp$v))
  if (rank == 1) {
    alpha_star <- matrix(alpha_star, ncol = 1)
  }
  alpha <- Pmin %*% alpha_star

  out <- list(alpha = alpha, alpha_star = alpha_star)
}
