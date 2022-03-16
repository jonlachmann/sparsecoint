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
SparseCointegration_Lasso <- function(data, p, r, alpha.init = NULL, beta.init = NULL, max.iter = 10, conv = 10^-2, lambda.gamma, rho.glasso = seq(from = 1, to = 0.1, length = 5), lambda_beta, cutoff = 0.8, glmnetthresh = 1e-04) {
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
  iter <- 1
  diff.obj <- 10 * conv
  Omega.init <- diag(1, q)
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter + 1)
  residuals <- Y - X %*% matrix(rep(diag(1, q), p - 1), ncol = q, byrow = T) - Z %*% beta.init %*% t(alpha.init)
  value.obj[1, ] <- (1 / n) * sum(diag(residuals %*% Omega.init %*% t(residuals))) - log(det(Omega.init))


  while ((iter < max.iter) & (diff.obj > conv)) {
    # Obtain Gamma
    gamma_fit <- nts.gamma.lassoreg(Y = Y, X = X, Z = Z, Pi = Pi.init, p = p, Omega = Omega.init, lambda.gamma = lambda.gamma, glmnetthresh = glmnetthresh)
    # Obtain Alpha
    alpha_fit <- nts.alpha.procrusted(Y = Y, X = X, Z = Z, zbeta = gamma_fit$gamma, r = r, Omega = Omega.init, P = gamma_fit$P, beta = beta.init)
    # Obtain Beta and Omega
    beta_fit <- nts.beta(Y = Y, X = X, Z = Z, zbeta = gamma_fit$gamma, rank = r, P = gamma_fit$P, alpha = alpha_fit$alpha, alphastar = alpha_fit$alpha_star, lambda_grid = lambda_beta, rho.glasso = rho.glasso, cutoff = cutoff, glmnetthresh = glmnetthresh)

    # Check convergence
    beta.new <- beta_fit$beta
    residuals <- Y - X %*% gamma_fit$gamma - Z %*% beta_fit$beta %*% t(alpha_fit$alpha)
    value.obj[1 + iter, ] <- (1 / n) * sum(diag((residuals) %*% beta_fit$omega %*% t(residuals))) - log(det(beta_fit$omega))
    diff.obj <- abs(value.obj[1 + iter, ] - value.obj[iter, ]) / abs(value.obj[iter, ])
    alpha.init <- alpha_fit$alpha
    beta.init <- beta.new
    Pi.init <- alpha.init %*% t(beta.new)
    Omega.init <- beta_fit$omega
    iter <- iter + 1
  }

  out <- list(beta = beta_fit$beta,
              alpha = alpha_fit$alpha,
              iter = iter,
              gamma = gamma_fit$gamma,
              omega = beta_fit$omega,
              beta.lambda=beta_fit$lambda,
              gamma.lambda=gamma_fit$lambda)
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