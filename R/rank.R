#' Function to determine the cointegration rank using the rank selection criterion
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param r.init initial value of cointegration rank
#' @param p number of lags to be included
#' @param max.iter.lasso maximum number of iterations PML
#' @param conv.lasso convergence parameter
#' @param max.iter.r maximu number of iterations to compute cointegration rank
#' @param beta.init initial value for beta
#' @param alpha.init initial value for alpha
#' @param rho.glasso tuning parameter for inverse covariance matrix
#' @param lambda.gamma tuning parameter for GAMMA
#' @param lambda.beta tuning parameter for BETA
#' @param tol tolerance parameter glmnet function
#' @return A List containing;
#' rhat: estimated cointegration rank
#' it.r: number of iterations
#' rhat_iterations: estimate of cointegration rank in each iteration
#' mu: value of mu
#' decomp: eigenvalue decomposition
determineRank <- function(data, r.init = NULL, p, max.iter.lasso = 3, conv.lasso = 10^-2, max.iter.r = 5, beta.init, alpha.init, rho.glasso = 0.1, lambda.gamma, lambda_beta, intercept=FALSE, tol = 1e-04) {
  Y <- data$diff
  X <- data$diff_lag
  Z <- data$level

  # Starting values
  if (is.null(beta.init) & is.null(alpha.init)) {
    q <- ncol(Y)
    beta.init <- matrix(1, ncol = 1, nrow = q)
    alpha.init <- initAlpha(Y, Z, X, beta.init, p)
    beta.init <- initBetaFinal(Y, Z, X, alpha.init, p, q, 1)
    beta.init <- apply(beta.init, 2, normalisationUnit)
  }

  # Preliminaries
  if (is.null(r.init)) {
    r.init <- ncol(Y)
  }
  diff.r <- 1
  iter <- 1
  rhat_iterations <- matrix(ncol = max.iter.r, nrow = 1)

  while ((iter < max.iter.r) & (diff.r > 0)) {
    # Initialization
    beta.init.fit <- matrix(0, ncol = r.init, nrow = ncol(Y))
    alpha.init.fit <- matrix(0, ncol = r.init, nrow = ncol(Y))
    if (r.init == 0) {
      diff.r <- 0
    } else {
      if (ncol(beta.init.fit) > ncol(beta.init)) {
        beta.init.fit[, seq_len(ncol(beta.init))] <- beta.init
        alpha.init.fit[, seq_len(ncol(beta.init))] <- alpha.init
      } else {
        beta.init.fit[, 1:r.init] <- beta.init[, 1:r.init]
        alpha.init.fit[, 1:r.init] <- alpha.init[, 1:r.init]
      }

      # Sparse cointegration fit
      lasso_fit <- rankSelectionCriterion(Y = Y, X = X, Z = Z, Beta = beta.init.fit, alpha = alpha.init.fit, p = p, r = r.init, max.iter = max.iter.lasso, conv = conv.lasso, lambda.gamma = lambda.gamma, lambda_beta = lambda_beta, rho.glasso = rho.glasso, intercept=intercept, tol = tol)

      # Response and predictors in penalized reduced rank regression
      khat <- calculateKhat(Y, X, Z, lasso_fit$zbeta, intercept)

      # Convergence checking
      rhat_iterations[, iter] <- khat
      diff.r <- abs(r.init - khat)
      r.init <- khat
      iter <- iter + 1
    }
  }

  out <- list(rhat = khat, iter = iter, rhat_iterations = rhat_iterations)
}

calculateKhat <- function (Y, X, Z, zbeta, intercept) {
  if (intercept) X <- cbind(1, X)
  Y.new <- Y - X %*% zbeta
  X.new <- Z
  M <- t(X.new) %*% X.new
  P <- X.new %*% ginv(M) %*% t(X.new)
  l <- ncol(Y.new)
  h <- qr(X.new)$rank
  m <- nrow(Y.new)
  S <- sum((Y.new - P %*% Y.new)^2) / (m * l - h * l)
  mu <- 2 * S * (l + h)
  decomp <- Re(eigen(t(Y.new) %*% P %*% Y.new)$values)
  khat <- length(which(decomp > mu | decomp == mu))
  return(khat)
}

#' Sparse cointegration function used in determine_rank function for Rank Selection Criterion
#' @param p number of lagged differences
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param r cointegration rank
#' @param alpha initial value for adjustment coefficients
#' @param Beta initial value for cointegrating vector
#' @param max.iter maximum number of iterations
#' @param conv convergence parameter
#' @param lambda.gamma tuning paramter short-run effects
#' @param lambda_beta tuning paramter cointegrating vector
#' @param rho.glasso tuning parameter inverse error covariance matrix
#' @param cutoff cutoff value time series cross-validation approach
#' @param tol tolerance parameter glmnet function
#' @return A list containing:
#' BETAhat: estimate of cointegrating vectors
#' ALPHAhat: estimate of adjustment coefficients
#' ZBETA: estimate of short-run effects
#' OMEGA: estimate of inverse covariance matrix
rankSelectionCriterion <- function(p, Y, X, Z, r, alpha = NULL, Beta, max.iter = 25, conv = 10^-3, lambda.gamma = 0.001, lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), rho.glasso = 0.5, cutoff = 0.8, intercept = F, tol = 1e-04) {
  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]

  # Starting value
  if (is.null(alpha)) {
    alpha <- matrix(rnorm(q * r), ncol = r, nrow = q)
    Beta <- matrix(rnorm(q * r), ncol = r, nrow = q)
    Pi <- alpha %*% t(Beta)
  } else {
    Pi <- alpha %*% t(Beta)
  }

  # Convergence parameters: initialization
  iter <- 1
  diff.obj <- 10 * conv
  Omega <- diag(1, q)
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter + 1)
  residuals <- calcResiduals(Y, X, Z, p, Beta, alpha, intercept)
  value.obj[1, ] <- (1 / n) * sum(diag(residuals %*% Omega %*% t(residuals))) - log(det(Omega))


  while ((iter < max.iter) & (diff.obj > conv)) {

    # Obtain Gamma, fixed value of tuning parameter
    gamma_fit <- nts.gamma.lassoreg(Y = Y, X = X, Z = Z,
                                          Pi = Pi, Omega = Omega,
                                          p = p, lambda.gamma = lambda.gamma, intercept=intercept, tol = tol, fixed=TRUE)
    # Obtain alpha
    alpha_fit <- nts.alpha.procrusted(Y = Y, X = X, Z = Z,
                                      zbeta = gamma_fit$gamma, Omega = Omega,
                                      r = r, P = gamma_fit$P, intercept=intercept, beta = Beta)
    # Obtain beta and omega
    beta_fit <- nts.beta(Y = Y, X = X, Z = Z,
                         zbeta = gamma_fit$gamma, alpha = alpha_fit$alpha, alphastar = alpha_fit$alpha_star,
                         lambda_grid = lambda_beta, rho_omega = rho.glasso,
                         rank = r, P = gamma_fit$P, cutoff = cutoff, intercept=intercept, tol = tol)


    # Check convergence
    beta.new <- beta_fit$beta
    residuals <- calcResiduals(Y, X, Z, gamma_fit$gamma, beta_fit$beta, alpha_fit$alpha, intercept)

    value.obj[1 + iter, ] <- (1 / n) * sum(diag((residuals) %*% beta_fit$omega %*% t(residuals))) - log(det(beta_fit$omega))
    diff.obj <- abs(value.obj[1 + iter, ] - value.obj[iter, ]) / abs(value.obj[iter, ])

    alpha <- alpha_fit$alpha
    Beta <- beta.new
    Pi <- alpha %*% t(beta.new) ###### Pi.init!!!
    Omega <- beta_fit$omega
    iter <- iter + 1
  }

  out <- list(beta = beta_fit$beta, alpha = alpha_fit$alpha, it = iter, zbeta = gamma_fit$gamma, omega = beta_fit$omega)
}
