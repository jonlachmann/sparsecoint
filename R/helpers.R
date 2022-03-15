CVscore.Reg <- function(U, X.data, Y.data, Y.sd) {
  # AUXILIARY FUNCTION
  # Standardized Mean Squared  Error
  return(mean((Y.data - X.data %*% U)^2 / Y.sd))
}

normalisationUnit <- function(U) {
  # AUXILIARY FUNCTION
  # output: normalized vector U
  length.U <- as.numeric(sqrt(t(U) %*% U))
  if (length.U == 0) {
    length.U <- 1
  }
  Uunit <- U / length.U
}

availableLambda <- function (U) {
  if (all(is.na(U) == T)) {
    return(F)
  } else {
    return(T)
  }
}

lagNames <- function (varnames, p) {
  lagnames <- character()
  for (i in seq_len(p)) {
    lagnames <- c(lagnames, paste0(varnames, ".", i))
  }
  return(lagnames)
}

#' Restrict the lambda sequence: exclude all zero-solution
restrictLambda <- function (Ymatrixsd, Xmatrix, lambda, glmnetthresh) {
  lassofit <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = lambda, family = "gaussian", thresh = glmnetthresh)
  restricted <- matrix(lambda[1, which(lassofit$df != 0)], nrow = 1)
  if (length(which(lassofit$df != 0)) <= 1) {
    lassofit <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
    restricted <- matrix(lassofit$lambda, nrow = 1)
  }
  return(restricted)
}

#' Run cross validation to determine the optimal lambda value
crossValidate <- function (cutoff, Ymatrix, Xmatrix, lambdas, glmnetthresh) {
  # Time series cross-validation to determine value of the tuning parameter lambda_grid
  cutoff.n <- round(cutoff * nrow(Ymatrix))
  cvscore <- matrix(NA, ncol = ncol(lambdas), nrow = nrow(Ymatrix) - cutoff.n)

  for (i in cutoff.n:nrow(Ymatrix) - 1) { # Loop to calculate cross-validation score
    # Training Data
    Ytrain <- Ymatrix[1:i, , drop=F]
    Ytrain.sd <- sd(Ytrain)
    Ytrain.scaled <- Ytrain / Ytrain.sd
    Xtrain <- Xmatrix[1:i, ]

    # Test Data
    Ytest <- Ymatrix[i + 1, , drop=F]
    Xtest <- Xmatrix[i + 1, ]

    # Estimation
    lasso_fit <- glmnet(y = Ytrain.scaled, x = Xtrain, lambda = lambdas, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
    coefs_std <- matrix(lasso_fit$beta[, which(lasso_fit$df != 0)], nrow = ncol(Xtrain)) # BETA in standardized scale
    coefs_org <- coefs_std * Ytrain.sd

    coefs_cv <- apply(coefs_org, 2, CVscore.Reg, Y.data = Ytest, X.data = Xtest, Y.sd = Ytrain.sd) # CVscore
    cvscore[i - cutoff.n + 1, which(lasso_fit$df != 0)] <- coefs_cv
  }

  CVscore.AVAILABLE <- as.matrix(cvscore[, apply(cvscore, 2, availableLambda)])
  lambdas <- lambdas[apply(cvscore, 2, availableLambda)]
  lambda.opt <- lambdas[which.min(apply(CVscore.AVAILABLE, 2, mean, na.rm=TRUE))]

  return(lambda.opt)
}