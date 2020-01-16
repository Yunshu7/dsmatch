estforboot_covQTT <- function(data, indices, p){
  data1 <- data[indices, ]

  dim2 <- dim(data1)[2]
  Y <- data1[, 1]
  A <- data1[, 2]
  Kiw.x <- data1[, 3]
  X <- data1[, 4:dim2]
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)

  lm.y <- Y
  lm.x <- cbind(X, t(apply(X, 1, combn, 2, prod)), X ^ 2, X ^ 3)

  lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
  # mu1.seive <- cbind(1, lm.x) %*% lm.out1$coefficients
  mu1.x <- cbind(1, lm.x[,which(!is.na(lm.out1$coefficients)) - 1]) %*% lm.out1$coefficients[!is.na(lm.out1$coefficients)]
  sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
  lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
  # mu0.seive <- cbind(1, lm.x) %*% lm.out0$coefficients
  mu0.x <- cbind(1, lm.x[,which(!is.na(lm.out0$coefficients)) - 1]) %*% lm.out0$coefficients[!is.na(lm.out0$coefficients)]
  sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

  mu.x <- A * mu1.x + (1 - A) * mu0.x

  f0 <- function(par) {
    mean(pnorm((par - mu0.x[loc.1]) / sqrt(sigsqhat0)) - p) +
      sum(Kiw.x[loc.0] * ((Y[loc.0] < par) - pnorm((par - mu0.x[loc.0]) / sqrt(sigsqhat0)))) / length(loc.1)
  }
  bootq0 <- bisect(f0, lo = min(Y), hi = max(Y), ytol = 1e-12, itmax = 100)

  f1 <- function(par) {
    mean(Y[loc.1] < par) - p
  }
  bootq1 <- bisect(f1, lo = min(Y[loc.1]), hi = max(Y[loc.1]), ytol = 1e-12, itmax = 100)

  return(bootq1 - bootq0)
}



