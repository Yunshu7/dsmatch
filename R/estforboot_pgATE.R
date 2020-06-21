estforboot_pgATE <- function(data, indices, model.pg, pg = NULL, lp.pg = NULL){
  data1 <- data[indices, ]
  
  dim2 <- dim(data1)[2]
  Y <- data1[, 1]
  A <- data1[, 2]
  Kiw.pg <- data1[, 3]
  Yiw.pg <- data1[, 4:5]
  X <- data1[, 6:dim2]
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)
  
  # if prognostic score model is given, then fit the model using bootstrap data
  # especailly when linear predictors are given, we should boostrap these predictors too
  # if model is not given, directly boostrap the given score
  if (model.pg == "glm"){
    lm1.out <- lm(Y[loc.1] ~ X[loc.1, ])
    lm0.out <- lm(Y[loc.0] ~ X[loc.0, ])
    mu1 <- cbind(1, X) %*% lm1.out$coefficients
    mu0 <- cbind(1, X) %*% lm0.out$coefficients
  }else if (model.pg == "glm_logit"){
    glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "logit"))
    mu1 <- cbind(1, X) %*% glm.out1$coefficients
    glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "logit"))
    mu0 <- cbind(1, X) %*% glm.out0$coefficients
  }else if (model.pg == "glm_probit"){
    glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "probit"))
    mu1 <- cbind(1, X) %*% glm.out1$coefficients
    glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "probit"))
    mu0 <- cbind(1, X) %*% glm.out0$coefficients
  }else if (model.pg == "linpred"){
    if (is.vector(lp.pg)){
      lp.pg <- lp.pg[indices]
      lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1])
      lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0])
      mu1 <- cbind(1, lp.pg) %*% lm1.out$coefficients
      mu0 <- cbind(1, lp.pg) %*% lm0.out$coefficients
    }else{
      lp.pg <- lp.pg[indices, ]
      lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1, ])
      lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0, ])
      mu1 <- cbind(1, lp.pg) %*% lm1.out$coefficients
      mu0 <- cbind(1, lp.pg) %*% lm0.out$coefficients
    }
  }else{
    mu0 = pg[indices, 1]
    mu1 = pg[indices, 2]
  }
  
  lm.y <- Y
  lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1)
  lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
  mu1.pg <- cbind(1, lm.x) %*% lm.out1$coefficients
  lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
  mu0.pg <- cbind(1, lm.x) %*% lm.out0$coefficients
  
  mu.pg <- A * mu1.pg + (1 - A) * mu0.pg
  
  boot.pg <- mean(mu1.pg - mu0.pg + (2 * A - 1) * (Kiw.pg + 1) * (Y - mu.pg))
  
  return(boot.pg)
}



