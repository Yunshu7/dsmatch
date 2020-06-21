# source("bisect.R")

estforboot_psQTE <- function(data, indices, p, model.ps, ps = NULL, lp.ps = NULL){
  data1 <- data[indices, ]
  
  dim2 <- dim(data1)[2]
  Y <- data1[, 1]
  A <- data1[, 2]
  Kiw.ps <- data1[, 3]
  Yiw.ps <- data1[, 4:5]
  X <- data1[, 6:dim2]
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)

  # if propensity score model is given, then fit the model using bootstrap data
  # especailly when linear predictors are given, we should boostrap these predictors too
  # if model is not given, directly boostrap the given score
  if (model.ps == "logit"){
    glm.out <- glm(A ~ X, family = binomial(link = "logit"))
    ps <- cbind(1, X) %*% glm.out$coefficients
  }else if (model.ps == "probit"){
    glm.out <- glm(A ~ X, family = binomial(link = "probit"))
    ps <- cbind(1, X) %*% glm.out$coefficients
  }else if (model.ps == "linpred"){
    if (is.vector(lp.ps)){
      lp.ps <- lp.ps[indices]
    }else{
      lp.ps <- lp.ps[indices, ]
    }
    glm.out <- glm(A ~ lp.ps, family = binomial)
    ps <- cbind(1, lp.ps) %*% glm.out$coefficients
  }else {
    ps <- ps[indices]
  }
  
  lm.y <- Y
  lm.x <- cbind(ps, ps ^ 2, ps ^ 3)
  lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
  mu1.ps <- cbind(1, lm.x) %*% lm.out1$coefficients
  sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
  lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
  mu0.ps <- cbind(1, lm.x) %*% lm.out0$coefficients
  sigsqhat0 <- mean((lm.out0$residuals) ^ 2)
  
  mu.ps <- A * mu1.ps + (1 - A) * mu0.ps
  
  # f0 <- function(par) {
  #   mean((1 - A) * (Kiw.ps + 1) * (Y < par)) - p
  # }
  # estq0 <- bisect(f0, lo = min(Yiw.ps[, 1]), hi = max(Yiw.ps[, 1]), ytol = 1e-12, itmax = 100)
  # 
  # f1 <- function(par) {
  #   mean((A) * (Kiw.ps + 1) * (Y < par)) - p
  # }
  # estq1 <- bisect(f1, lo = min(Yiw.ps[, 2]), hi = max(Yiw.ps[, 2]), ytol = 1e-12, itmax = 100)
  # 
  # 
  # boot.ps <- estq1 - estq0
  # 
  # return(boot.ps)
  f0 <- function(par) {
    mean(pnorm((par - mu0.ps) / sqrt(sigsqhat0)) - p) +
      mean((1 - A) * (Kiw.ps + 1) * ((Y < par) - pnorm((par - mu0.ps) / sqrt(sigsqhat0))))
  }
  bootq0 <- bisect(f0, lo = min(Yiw.ps[, 1]), hi = max(Yiw.ps[, 1]), ytol = 1e-12, itmax = 100)
  
  f1 <- function(par) {
    mean(pnorm((par - mu1.ps) / sqrt(sigsqhat1)) - p) +
      mean((A) * (Kiw.ps + 1) * ((Y < par) - pnorm((par - mu1.ps) / sqrt(sigsqhat1))))
  }
  bootq1 <- bisect(f1, lo = min(Yiw.ps[, 2]), hi = max(Yiw.ps[, 2]), ytol = 1e-12, itmax = 100)
  
  return(bootq1 - bootq0)
}



