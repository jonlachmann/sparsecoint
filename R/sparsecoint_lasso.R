#' Main function to perform sparse cointegration
#' @param p number of lagged differences
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param r cointegration rank
#' @param alpha initial value for adjustment coefficients
#' @param beta initial value for cointegrating vector
#' @param max.iter maximum number of iterations
#' @param conv convergence parameter
#' @param lambda_gamma tuning paramter short-run effects
#' @param lambda_beta tuning paramter cointegrating vector
#' @param rho_omega tuning parameter inverse error covariance matrix
#' @param cutoff cutoff value time series cross-validation approach
#' @param tol tolerance parameter glmnet function
#' @return A list containing:
#' beta: estimate of cointegrating vectors
#' alpha: estimate of adjustment coefficients
#' gamma: estimate of short-run effects
#' omega: estimate of inverse covariance matrix
SparseCointegration_Lasso <- function(data, p, r, alpha = NULL, beta = NULL, max.iter = 10, conv = 10^-2, lambda_gamma, rho_omega, lambda_beta, cutoff = 0.8, tol = 1e-04, intercept=FALSE) {
  Y <- data$diff
  X <- data$diff_lag
  Z <- data$level

  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]

  # Obtain starting values
  init <- sparseCointegrationInit(Y, X, Z, r, p, q, alpha, beta)

  beta <- init$beta
  alpha <- init$alpha

  # Convergence parameters: initialization
  iter <- 1
  diff.obj <- 10 * conv
  Omega <- diag(1, q)
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter + 1)
  residuals <- calcResiduals(Y, X, Z, p, beta, alpha, intercept)
  value.obj[1, ] <- (1 / n) * sum(diag(residuals %*% Omega %*% t(residuals))) - log(det(Omega))


  while ((iter < max.iter) & (diff.obj > conv)) {
    fit <- sparseCointegrationFit(Y, Z, X, alpha, Omega, beta, p, r, lambda_gamma, lambda_beta, rho_omega, cutoff, intercept, tol)

    # Check convergence
    residuals <- calcResiduals(Y, X, Z, fit$gamma, fit$beta, fit$alpha, intercept)
    value.obj[1 + iter, ] <- (1 / n) * sum(diag((residuals) %*% fit$omega %*% t(residuals))) - log(det(fit$omega))
    diff.obj <- abs(value.obj[1 + iter, ] - value.obj[iter, ]) / abs(value.obj[iter, ])
    alpha <- fit$alpha
    beta <- fit$beta
    Omega <- fit$omega
    iter <- iter + 1
  }

  fit$iter <- iter

  return(fit)
}

sparseCointegrationFit <- function (Y, Z, X, alpha, omega, beta, p, rank, lambda_gamma, lambda_beta, omega_rho, cutoff, intercept, tol, fixed=FALSE) {
  # Obtain Gamma
  gamma_fit <- nts.gamma.lassoreg(Y = Y, X = X, Z = Z,
                                  alpha = alpha, beta = beta, Omega = omega,
                                  p = p, lambda = lambda_gamma, intercept = intercept, tol = tol, fixed=fixed)
  # Obtain Alpha
  alpha_fit <- nts.alpha.procrusted(Y = Y, X = X, Z = Z,
                                    zbeta = gamma_fit$gamma, Omega = omega,
                                    rank = rank, P = gamma_fit$P, beta = beta, intercept = intercept)
  # Obtain Beta and Omega
  beta_fit <- nts.beta(Y = Y, X = X, Z = Z,
                       gamma = gamma_fit$gamma, alpha = alpha_fit$alpha, alphastar = alpha_fit$alpha_star,
                       lambda = lambda_beta, rho_omega = omega_rho,
                       rank = rank, P = gamma_fit$P, cutoff = cutoff, intercept=intercept, tol = tol)

  return(list(gamma=gamma_fit$gamma, alpha=alpha_fit$alpha, beta=beta_fit$beta, omega=beta_fit$omega,
              gamma_lambda=gamma_fit$lambda, beta_lambda=beta_fit$lambda, omega_rho=beta_fit$rho))
}

sparseCointegrationInit <- function (Y, X, Z, rank, p, q, alpha.init=NULL, beta.init=NULL) {
  if (is.null(alpha.init) || is.null(beta.init)) {
    if (is.null(alpha.init) && is.null(beta.init)) {
      beta.init <- initBeta(Y, Z, rank)
      alpha.init <- initAlpha(Y, Z, X, beta.init, p)
    } else if (is.null(alpha.init)) {
      alpha.init <- initAlpha(Y, Z, X, beta.init, p)
    }
    beta.init <- initBetaFinal(Y, Z, X, alpha.init, p, q, rank)
  }
  return(list(alpha=alpha.init, beta=beta.init))
}


initBeta <- function (Y, Z, rank) {
  SigmaYY <- diag(apply(Y, 2, var))
  SigmaZZ <- diag(apply(Z, 2, var))

  SigmaYZ <- cov(Y, Z)
  SigmaZY <- cov(Z, Y)

  matrix <- solve(SigmaZZ) %*% SigmaZY %*% solve(SigmaYY) %*% SigmaYZ
  decomp <- eigen(matrix)
  beta.init <- decomp$vectors[, 1:rank]
  return(beta.init)
}

initBetaFinal <- function (Y, Z, X, alpha.init, p, q, rank) {
  beta.init <- matrix(NA, ncol = rank, nrow = q)
  Yinit <- (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)) %*% alpha.init
  for (i in 1:rank) {
    fit <- lars(x = Z, y = Yinit[, i], type = "lasso", normalize = F, intercept = F)$beta[-1, ]
    beta.init[, i] <- fit[nrow(fit), ]
  }
  return(beta.init)
}

initAlpha <- function (Y, Z, X, beta.init, p) {
  SVDinit <- svd(t(Z %*% beta.init) %*% (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)))
  alpha.init <- t(SVDinit$u %*% t(SVDinit$v))
  return(alpha.init)
}