CVscore.Reg <- function(U, X.data, Y.data, Y.sd) {
  # AUXILIARY FUNCTION
  # Standardized Mean Squared  Error
  return(mean((Y.data - X.data %*% U)^2 / Y.sd))
}

NORMALIZATION_UNIT <- function(U) {
  # AUXILIARY FUNCTION
  # output: normalized vector U
  length.U <- as.numeric(sqrt(t(U) %*% U))
  if (length.U == 0) {
    length.U <- 1
  }
  Uunit <- U / length.U
}

AVAILABLE_LAMBDA <- function(U) {
  # AUXILIARY FUNCTION
  if (all(is.na(U) == T)) {
    return(F)
  } else {
    return(T)
  }
}

AVERAGE <- function(U) {
  # AUXILIARY FUNCTION
  mean(U[which(is.na(U) == F)])
}