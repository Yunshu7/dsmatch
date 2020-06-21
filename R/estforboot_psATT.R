estforboot_psATT <- function(data, indices, model.ps, ps = NULL, lp.ps = NULL){
  data1 <- data[indices, ]

  dim2 <- dim(data1)[2]
  Y <- data1[, 1]
  A <- data1[, 2]
  Kiw.ps <- data1[, 3]
  X <- data1[, 4:dim2]
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
  lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
  mu0.ps <- cbind(1, lm.x) %*% lm.out0$coefficients

  mu.ps <- A * mu1.ps + (1 - A) * mu0.ps

  boot.ps <- sum(A * (mu1.ps - mu0.ps) + (A - (1 - A) * Kiw.ps) * (Y - mu.ps)) / length(loc.1)

  return(boot.ps)
}



