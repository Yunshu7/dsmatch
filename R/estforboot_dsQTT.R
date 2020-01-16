estforboot_dsQTT <- function(data, indices, p, model.ps, ps = NULL, lp.ps = NULL, model.pg, pg = NULL, lp.pg = NULL){
  data1 <- data[indices, ]

  dim2 <- dim(data1)[2]
  Y <- data1[, 1]
  A <- data1[, 2]
  Kiw.ds <- data1[, 3]
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
    lp.ps <- lp.ps[indices, ]
    glm.out <- glm(A ~ lp.ps, family = binomial)
    ps <- cbind(1, lp.ps) %*% glm.out$coefficients
  }else {
    ps <- ps[indices]
  }

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
    lp.pg <- lp.pg[indices, ]
    lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1, ])
    lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0, ])
    mu1 <- cbind(1, lp.pg) %*% lm1.out$coefficients
    mu0 <- cbind(1, lp.pg) %*% lm0.out$coefficients
  }else if(model.pg == "zir_logit"){
    # location where outcome is not zero
    # each for whole set, treatment group and control group
    loc.r <- which(Y != 0)
    loc.1r <- intersect(loc.r, loc.1)
    loc.0r <- intersect(loc.r, loc.0)

    # change non-zero outcomes into 1
    # then run a glm model with a binomial family
    # use linear predictors as part of the prognostic score
    # calculate non-zero probability to estimate \mu
    glm.y1b <- Y
    glm.y1b[loc.1r] <- 1
    glm.y1b <- glm.y1b[loc.1]
    glm.x1b <- X[loc.1, ]
    glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "logit"))
    p1 <- cbind(1, X) %*% glm.out1$coefficients
    pb1 <- 1 / (1 + exp(-p1))

    glm.y0b <- Y
    glm.y0b[loc.0r] <- 1
    glm.y0b <- glm.y0b[loc.0]
    glm.x0b <- X[loc.0, ]
    glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "logit"))
    p0 <- cbind(1, X) %*% glm.out0$coefficients
    pb0 <- 1 / (1 + exp(-p0))

    # run regression model for non-zero outcomes
    lm1.out <- lm(Y[loc.1r] ~ X[loc.1r, ])
    lm0.out <- lm(Y[loc.0r] ~ X[loc.0r, ])
    mu1 <- cbind(1, X) %*% lm1.out$coefficients
    mu0 <- cbind(1, X) %*% lm0.out$coefficients
  }else if(model.pg == "zir_probit"){
    # location where outcome is not zero
    # each for whole set, treatment group and control group
    loc.r <- which(Y != 0)
    loc.1r <- intersect(loc.r, loc.1)
    loc.0r <- intersect(loc.r, loc.0)

    # change non-zero outcomes into 1
    # then run a glm model with a binomial family
    # use linear predictors as part of the prognostic score
    # calculate non-zero probability to estimate \mu
    glm.y1b <- Y
    glm.y1b[loc.1r] <- 1
    glm.y1b <- glm.y1b[loc.1]
    glm.x1b <- X[loc.1, ]
    glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "probit"))
    p1 <- cbind(1, X) %*% glm.out1$coefficients
    pb1 <- pnorm(p1)

    glm.y0b <- Y
    glm.y0b[loc.0r] <- 1
    glm.y0b <- glm.y0b[loc.0]
    glm.x0b <- X[loc.0, ]
    glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "probit"))
    p0 <- cbind(1, X) %*% glm.out0$coefficients
    pb0 <- pnorm(p0)

    # run regression model for non-zero outcomes
    lm1.out <- lm(Y[loc.1r] ~ X[loc.1r, ])
    lm0.out <- lm(Y[loc.0r] ~ X[loc.0r, ])
    mu1 <- cbind(1, X) %*% lm1.out$coefficients
    mu0 <- cbind(1, X) %*% lm0.out$coefficients
  }else{
    mu0 = pg[indices, 1]
    mu1 = pg[indices, 2]
  }

  if (!grepl("zir", model.pg)) {
    lm.y <- Y
    lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1,
      ps, ps ^ 2, ps * mu0, ps * mu1, ps * mu0 * mu1)
    lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
    mu1.ds <- cbind(1, lm.x) %*% lm.out1$coefficients
    sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
    lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
    mu0.ds <- cbind(1, lm.x) %*% lm.out0$coefficients
    sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

    f0 <- function(par) {
      mean(pnorm((par - mu0.ds[loc.1]) / sqrt(sigsqhat0)) - p) +
        sum(Kiw.ds[loc.0] * ((Y[loc.0] < par) - pnorm((par - mu0.ds[loc.0]) / sqrt(sigsqhat0)))) / length(loc.1)
    }
    bootq0 <- bisect(f0, lo = min(Y), hi = max(Y), ytol = 1e-12, itmax = 100)

    f1 <- function(par) {
      mean(Y[loc.1] < par) - p
    }
    bootq1 <- bisect(f1, lo = min(Y[loc.1]), hi = max(Y[loc.1]), ytol = 1e-12, itmax = 100)

  }else{
    lm.y <- Y
    lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1,
      p0, p1, p0 ^ 2, p1 ^ 2, p0 * p1,
      ps, ps ^ 2, ps * mu0, ps * mu1,
      p0 * mu0, p0 * mu1, p1 * mu0, p1 * mu1, p0 * ps, p1 * ps)

    lm.out1 <- lm(lm.y[loc.1r] ~ lm.x[loc.1r, ])
    mu1.seive <- cbind(1, lm.x) %*% lm.out1$coefficients
    sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
    lm.out0 <- lm(lm.y[loc.0r] ~ lm.x[loc.0r, ])
    mu0.seive <- cbind(1, lm.x) %*% lm.out0$coefficients
    sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

    if(model.pg == "zir_logit"){
      glm.x1b <- lm.x[loc.1, ]
      glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "logit"))
      p1.seive <- cbind(1, lm.x) %*% glm.out1$coefficients
      pb1.seive <- 1 / (1 + exp(-p1.seive))

      glm.x0b <- lm.x[loc.0, ]
      glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "logit"))
      p0.seive <- cbind(1, lm.x) %*% glm.out0$coefficients
      pb0.seive <- 1 / (1 + exp(-p0.seive))
    }else if(model.pg == "zir_probit"){
      glm.x1b <- lm.x[loc.1, ]
      glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "probit"))
      p1.seive <- cbind(1, lm.x) %*% glm.out1$coefficients
      pb1.seive <- pnorm(p1.seive)

      glm.x0b <- lm.x[loc.0, ]
      glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "probit"))
      p0.seive <- cbind(1, lm.x) %*% glm.out0$coefficients
      pb0.seive <- pnorm(p0.seive)
    }

    mu1.ds <- mu1.seive
    mu0.ds <- mu0.seive
    pb1.ds <- pb1.seive
    pb0.ds <- pb0.seive

    f0 <- function(par) {
      mean(pnorm((par - mu0.ds[loc.1]) / sqrt(sigsqhat0)) * pb0.ds[loc.1] + (par >= 0) * (1 - pb0.ds[loc.1]) - p) +
        sum(Kiw.ds * ((Y[loc.0] < par) - (pnorm((par - mu0.ds[loc.0]) / sqrt(sigsqhat0)) * pb0.ds[loc.0] + (par >= 0) * (1 - pb0.ds[loc.0])))) / length(loc.1)
    }
    bootq0 <- bisect(f0, lo = min(Y), hi = max(Y), ytol = 1e-12, itmax = 100)

    f1 <- function(par) {
      mean(Y[loc.1] < par) - p
    }
    bootq1 <- bisect(f1, lo = min(Y[loc.1]), hi = max(Y[loc.1]), ytol = 1e-12, itmax = 100)
  }

  return(bootq1 - bootq0)
}



