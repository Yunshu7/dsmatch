estforboot_pgQTT <- function(data, indices, p, model.pg, pg = NULL, lp.pg = NULL){
  data1 <- data[indices, ]

  dim2 <- dim(data1)[2]
  Y <- data1[, 1]
  A <- data1[, 2]
  Kiw.pg <- data1[, 3]
  X <- data1[, 4:dim2]
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)

  if (model.pg == "glm" && is.null(pg)){
    # lm1.out <- lm(Y[loc.1] ~ X[loc.1, ])
    # mu1 <- cbind(1, X)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
    lm0.out <- lm(Y[loc.0] ~ X[loc.0, ])
    mu0 <- cbind(1, X)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
  }else if (model.pg == "glm_logit" && is.null(pg)){
    # glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "logit"))
    # mu1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
    glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "logit"))
    mu0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
  }else if (model.pg == "glm_probit" && is.null(pg)){
    # glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "probit"))
    # mu1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
    glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "probit"))
    mu0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
  }else if (model.pg == "linpred"){
    # if linear predictors are valid, then fit a linear model
    if(is.vector(lp.pg)){
      # lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1])
      lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0])
      # mu1 <- cbind(1, lp.pg)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
      mu0 <- cbind(1, lp.pg)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    }else{
      # lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1, ])
      lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0, ])
      # mu1 <- cbind(1, lp.pg)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
      mu0 <- cbind(1, lp.pg)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    }
  }else{
    # mu0 = pg[indices, 1]
    # mu1 = pg[indices, 2]
    mu0 = pg[indices]
  }

  lm.y <- Y
  lm.x <- cbind(mu0, mu0 ^ 2, mu0 ^3)
  lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
  mu1.pg <- cbind(1, lm.x) %*% lm.out1$coefficients
  sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
  lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
  mu0.pg <- cbind(1, lm.x) %*% lm.out0$coefficients
  sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

  mu.pg <- A * mu1.pg + (1 - A) * mu0.pg

  f0 <- function(par) {
    mean(pnorm((par - mu0.pg[loc.1]) / sqrt(sigsqhat0)) - p) +
      sum(Kiw.pg[loc.0] * ((Y[loc.0] < par) - pnorm((par - mu0.pg[loc.0]) / sqrt(sigsqhat0)))) / length(loc.1)
  }
  bootq0 <- bisect(f0, lo = min(Y), hi = max(Y), ytol = 1e-12, itmax = 100)

  f1 <- function(par) {
    mean(Y[loc.1] < par) - p
  }
  bootq1 <- bisect(f1, lo = min(Y[loc.1]), hi = max(Y[loc.1]), ytol = 1e-12, itmax = 100)

  return(bootq1 - bootq0)
}



