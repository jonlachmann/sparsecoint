#' Main function to perform sparse cointegration
#' @param p: number of lagged differences
#' @param Y: Response Time Series
#' @param X: Time Series in Differences
#' @param Z: Time Series in Levels
#' @param r: cointegration rank
#' @param alpha.init: initial value for adjustment coefficients
#' @param beta.init: initial value for cointegrating vector
#' @param max.iter: maximum number of iterations
#' @param conv: convergence parameter
#' @param lambda.gamma: tuning paramter short-run effects
#' @param lambda_beta: tuning paramter cointegrating vector
#' @param rho.glasso: tuning parameter inverse error covariance matrix
#' @param cutoff: cutoff value time series cross-validation approach
#' @param glmnetthresh: tolerance parameter glmnet function
#' @return A list containing:
#' beta: estimate of cointegrating vectors
#' alpha: estimate of adjustment coefficients
#' gamma: estimate of short-run effects
#' omega: estimate of inverse covariance matrix
SparseCointegration_Lasso <- function(data, p, r, alpha.init = NULL, beta.init = NULL, max.iter = 10, conv = 10^-2, lambda.gamma = matrix(seq(from = 0.01, to = 0.001, length = 5), nrow = 1), rho.glasso = seq(from = 1, to = 0.1, length = 5), lambda_beta = matrix(seq(from = 0.1, to = 0.001, length = 100), nrow = 1), cutoff = 0.8, glmnetthresh = 1e-04) {
  Y <- data$diff
  X <- data$diff_lag
  Z <- data$level

  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]

  # Obtain starting values
  init <- sparseCointegrationInit(Y, X, Z, r, p, q, alpha.init, beta.init)

  beta.init <- init$beta
  alpha.init <- init$alpha
  Pi.init <- init$Pi

  # Convergence parameters: initialization
  it <- 1
  diff.obj <- 10 * conv
  Omega.init <- diag(1, q)
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter + 1)
  residuals <- Y - X %*% matrix(rep(diag(1, q), p - 1), ncol = q, byrow = T) - Z %*% beta.init %*% t(alpha.init)
  value.obj[1, ] <- (1 / n) * sum(diag(residuals %*% Omega.init %*% t(residuals))) - log(det(Omega.init))


  while ((it < max.iter) & (diff.obj > conv)) {
    # Obtain Gamma
    gamma.fit <- nts.gamma.lassoreg(Y = Y, X = X, Z = Z, Pi = Pi.init, p = p, Omega = Omega.init, lambda.gamma = lambda.gamma, glmnetthresh = glmnetthresh)
    # Obtain Alpha
    alpha.fit <- nts.alpha.procrusted(Y = Y, X = X, Z = Z, zbeta = gamma.fit$gamma, r = r, Omega = Omega.init, P = gamma.fit$P, beta = beta.init)
    # Obtain Beta and Omega
    beta.fit <- nts.beta(Y = Y, X = X, Z = Z, zbeta = gamma.fit$gamma, rank = r, P = gamma.fit$P, alpha = alpha.fit$ALPHA, alphastar = alpha.fit$ALPHAstar, lambda_grid = lambda_beta, rho.glasso = rho.glasso, cutoff = cutoff, glmnetthresh = glmnetthresh)

    # Check convergence
    beta.new <- beta.fit$beta
    residuals <- Y - X %*% gamma.fit$gamma - Z %*% beta.fit$beta %*% t(alpha.fit$ALPHA)
    value.obj[1 + it, ] <- (1 / n) * sum(diag((residuals) %*% beta.fit$omega %*% t(residuals))) - log(det(beta.fit$omega))
    diff.obj <- abs(value.obj[1 + it, ] - value.obj[it, ]) / abs(value.obj[it, ])
    alpha.init <- alpha.fit$ALPHA
    beta.init <- beta.new
    Pi.init <- alpha.init %*% t(beta.new)
    Omega.init <- beta.fit$omega
    it <- it + 1
  }

  out <- list(beta = beta.fit$beta,
              alpha = alpha.fit$ALPHA,
              it = it,
              gamma = gamma.fit$gamma,
              omega = beta.fit$omega,
              beta.lambda=beta.fit$lambda,
              gamma.lambda=gamma.fit$lambda)
}

sparseCointegrationInit <- function (Y, X, Z, rank, p, q, alpha.init=NULL, beta.init=NULL) {
  if (is.null(alpha.init) & is.null(beta.init)) {
    SigmaYY <- diag(apply(Y, 2, var))
    SigmaZZ <- diag(apply(Z, 2, var))

    SigmaYZ <- cov(Y, Z)
    SigmaZY <- cov(Z, Y)

    matrix <- solve(SigmaZZ) %*% SigmaZY %*% solve(SigmaYY) %*% SigmaYZ
    decomp <- eigen(matrix)
    beta.init <- decomp$vectors[, 1:rank]

    SVDinit <- svd(t(Z %*% beta.init) %*% (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)))
    alpha.init <- t(SVDinit$u %*% t(SVDinit$v))
    beta.init <- matrix(NA, ncol = rank, nrow = q)
    Yinit <- (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)) %*% alpha.init
    for (i in 1:rank) {
      FIT <- lars(x = Z, y = Yinit[, i], type = "lasso", normalize = F, intercept = F)$beta[-1, ]
      beta.init[, i] <- FIT[nrow(FIT), ]
    }
    Pi.init <- alpha.init %*% t(beta.init)
  }

  if (is.null(alpha.init) & !is.null(beta.init)) {
    SVDinit <- svd(t(Z %*% beta.init) %*% (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)))
    alpha.init <- t(SVDinit$u %*% t(SVDinit$v))
    beta.init <- matrix(NA, ncol = rank, nrow = q)
    Yinit <- (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)) %*% alpha.init
    for (i in 1:rank) {
      FIT <- lars(x = Z, y = Yinit[, i], type = "lasso", normalize = F, intercept = F)$beta[-1, ]
      beta.init[, i] <- FIT[nrow(FIT), ]
    }
    Pi.init <- alpha.init %*% t(beta.init)
  }

  if (!is.null(alpha.init) & !is.null(beta.init)) {
    Pi.init <- alpha.init %*% t(beta.init)
  }

  return(list(Pi=Pi.init, alpha=alpha.init, beta=beta.init))
}