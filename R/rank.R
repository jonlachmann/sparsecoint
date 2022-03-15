determine_rank <- function(data, r.init = NULL, p, max.iter.lasso = 3, conv.lasso = 10^-2, max.iter.r = 5, beta.init, alpha.init, rho.glasso = 0.1, lambda.gamma = matrix(seq(from = 0.1, to = 0.001, length = 10), nrow = 1), lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), glmnetthresh = 1e-04) {
  ### FUNCTION TO DETERMINE THE COINTEGRATION RANK USING THE RANK SELECTION CRITERION ###

  Y <- data$diff
  X <- data$diff_lag
  Z <- data$level

  # INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r.init: initial value of cointegration rank
  # p: number of lags to be included
  # max.iter.lasso: maximum number of iterations PML
  # conv.lasso: convergence parameter
  # max.iter.r: maximu number of iterations to compute cointegration rank
  # beta.init: initial value for beta
  # alpha.init: initial value for alpha
  # rho.glasso: tuning parameter for inverse covariance matrix
  # lambda.gamma: tuning parameter for GAMMA
  # lambda.beta: tuning parameter for BETA
  # glmnetthresh: tolerance parameter glmnet function

  # OUTPUT
  # rhat: estimated cointegration rank
  # it.r: number of iterations
  # rhat_iterations: estimate of cointegration rank in each iteration
  # mu: value of mu
  # decomp: eigenvalue decomposition

  # START CODE

  # Starting values
  if (is.null(beta.init) & is.null(alpha.init)) {
    q <- ncol(Y)
    beta.init <- matrix(1, ncol = 1, nrow = q)
    SVDinit <- svd(t(Z %*% beta.init) %*% (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)))
    alpha.init <- t(SVDinit$u %*% t(SVDinit$v))
    beta.init <- matrix(NA, ncol = 1, nrow = q)
    Yinit <- (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)) %*% alpha.init
    for (i in 1:1) {
      FIT <- lars(x = Z, y = Yinit[, i], type = "lasso", normalize = F, intercept = F)$beta[-1, ]
      beta.init[, i] <- FIT[nrow(FIT), ] # nrow(FIT)
    }
    beta.init <- apply(beta.init, 2, normalisationUnit)
  }

  # Preliminaries
  n <- nrow(Y)
  if (is.null(r.init)) {
    r.init <- ncol(Y)
  }
  diff.r <- 1
  it.r <- 1
  rhat_iterations <- matrix(ncol = max.iter.r, nrow = 1)


  while ((it.r < max.iter.r) & (diff.r > 0)) {
    # Initialization
    beta.init.fit <- matrix(0, ncol = r.init, nrow = ncol(Y))
    alpha.init.fit <- matrix(0, ncol = r.init, nrow = ncol(Y))
    if (r.init == 0) {
      beta.init.fit <- matrix(0, ncol = r.init, nrow = ncol(Y))
      alpha.init.fit <- matrix(0, ncol = r.init, nrow = ncol(Y))
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
      lasso_fit <- SparseCointegration_RSC(Y = Y, X = X, Z = Z, beta.init = beta.init.fit, alpha.init = alpha.init.fit, p = p, r = r.init, max.iter = max.iter.lasso, conv = conv.lasso, lambda.gamma = lambda.gamma, lambda_beta = lambda_beta, rho.glasso = rho.glasso, glmnetthresh = glmnetthresh)

      # Response  and predictors in penalized reduced rank regression
      Y.new <- Y - X %*% lasso_fit$ZBETA
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

      # Convergence checking
      rhat_iterations[, it.r] <- khat
      diff.r <- abs(r.init - khat)
      r.init <- khat
      it.r <- it.r + 1
    }
  }

  out <- list(rhat = khat, it.r = it.r, rhat_iterations = rhat_iterations, mu = mu, decomp = decomp)
}

SparseCointegration_RSC <- function(p, Y, X, Z, r, alpha.init = NULL, beta.init, max.iter = 25, conv = 10^-3, lambda.gamma = 0.001, lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), rho.glasso = 0.5, cutoff = 0.8, intercept = F, glmnetthresh = 1e-04) {
  ### Sparse cointegration function used in determine_rank function ###

  ## INPUT
  # p: number of lagged differences
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r: cointegration rank
  # alpha.init: initial value for adjustment coefficients
  # beta.init: initial value for cointegrating vector
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # lambda.gamma: tuning paramter short-run effects
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # glmnetthresh: tolerance parameter glmnet function

  ## OUTPUT
  # BETAhat: estimate of cointegrating vectors
  # ALPHAhat: estimate of adjustment coefficients
  # ZBETA: estimate of short-run effects
  # OMEGA: estimate of inverse covariance matrix


  ## START CODE

  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]

  # Starting value
  if (is.null(alpha.init)) {
    alpha.init <- matrix(rnorm(q * r), ncol = r, nrow = q)
    beta.init <- matrix(rnorm(q * r), ncol = r, nrow = q)
    Pi.init <- alpha.init %*% t(beta.init)
  } else {
    Pi.init <- alpha.init %*% t(beta.init)
  }

  # Convergence parameters: initialization
  it <- 1
  diff.obj <- 10 * conv
  Omega.init <- diag(1, q)
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter + 1)
  RESID <- Y - X %*% matrix(rep(diag(1, q), p - 1), ncol = q, byrow = T) - Z %*% beta.init %*% t(alpha.init)
  value.obj[1, ] <- (1 / n) * sum(diag(RESID %*% Omega.init %*% t(RESID))) - log(det(Omega.init))


  while ((it < max.iter) & (diff.obj > conv)) {

    # Obtain Gamma, fixed value of tuning parameter
    gamma_fit <- nts.gamma.fixed.lassoreg(Y = Y, X = X, Z = Z, Pi = Pi.init, p = p, Omega = Omega.init, lambda.gamma = lambda.gamma, glmnetthresh = glmnetthresh)
    # Obtain alpha
    alpha_fit <- nts.alpha.procrusted(Y = Y, X = X, Z = Z, zbeta = gamma_fit$zbeta, r = r, Omega = Omega.init, P = gamma_fit$P, beta = beta.init)
    # abtain beta and omega
    beta_fit <- nts.beta(Y = Y, X = X, Z = Z, zbeta = gamma_fit$zbeta, rank = r, P = gamma_fit$P, alpha = alpha_fit$ALPHA, alphastar = alpha_fit$ALPHAstar, lambda_grid = lambda_beta, rho.glasso = rho.glasso, cutoff = cutoff, glmnetthresh = glmnetthresh)


    # Check convergence
    beta.new <- beta_fit$beta
    if (intercept == T) {
      RESID <- Y - X %*% gamma_fit$zbeta[-1, ] - Z %*% beta_fit$beta %*% t(alpha_fit$ALPHA)
    } else {
      RESID <- Y - X %*% gamma_fit$zbeta - Z %*% beta_fit$beta %*% t(alpha_fit$ALPHA)
    }
    value.obj[1 + it, ] <- (1 / n) * sum(diag((RESID) %*% beta_fit$omega %*% t(RESID))) - log(det(beta_fit$omega))
    diff.obj <- abs(value.obj[1 + it, ] - value.obj[it, ]) / abs(value.obj[it, ])

    alpha.init <- alpha_fit$ALPHA
    beta.init <- beta.new
    Pi.init <- alpha.init %*% t(beta.new) ###### Pi.init!!!
    Omega.init <- beta_fit$omega
    it <- it + 1
  }

  out <- list(BETAhat = beta_fit$beta, ALPHAhat = alpha_fit$ALPHA, it = it, ZBETA = gamma_fit$zbeta, OMEGA = beta_fit$omega)
}