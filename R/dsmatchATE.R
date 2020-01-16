#' Double Score Matching Estimator for Average Treatment Effect.
#'
#' \code{dsmatchATE} applys matching algorithm to estimate average
#' treatment effect based on propensity score and prognostic score.
#' Classical matching algortihms, including propensity score matching,
#' prognositc score matching and matching directly on covarites are
#' also contained for comparison. Covariate balance results are also
#' provided.
#'
#' For both propensity socre and prognostic score, user should either
#' select a model or provide the score directly. If linear predictors
#' are used to fit a logistic model for propensity score or a linear
#' model for prognostic score, they should be determined by
#' \code{lp.ps} or \code{lp.pg} argument. If \code{model.ps} (and
#' \code{lp.ps} if linear predictors are used) is given, then
#' \code{ps} does not need to be specified, and vice versa. However,
#' if propensity socre is given by \code{ps} while \code{model.ps}
#' is chosen at the same time, the model will be ignored and matching
#' will be based on the score given by the user directly. A warning
#' will be thrown if this situation happens. Similar results for
#' prognostic score.
#'
#' A special model for prognostic score is the zero inflated regression
#'  model, which fits a logistic model for the probability to be zero,
#'  and a regression model for the non-zero values.
#'
#'
#' @param Y Outcome as numeric vector.
#' @param X Covarites as numeric vector or matrix.
#' @param A Treatment assignment as numeric vector with \code{1} stands for
#' treatment group and \code{0} stands for control group.
#' @param method Matching method to use, including \code{"dsm"} as double
#' score matching, \code{"ps"} as propensity score matching, \code{"pg"}
#' as prognostic score matching and \code{"cov"} as matching on covariates
#' directly.
#' @param model.ps Fitted model for propensity score, including
#' \code{"logit"} as logistic model, \code{"probit"} as probit model,
#' \code{"linpred"} as logistic model with linear predictors specified by
#' \code{"lp.ps"}. Don't need to be specified if \code{ps} is given.
#' @param ps Propensity score as numeric vector given by user. Don't need
#' to be specified if \code{model.ps} is given.
#' @param lp.ps Linear predictors for propensity score as numeric vector or
#' matrix. Don't need to be specified if \code{model.ps} is not
#' \code{"linpred"}.
#' @param model.pg Fitted model for prognostic score, including
#' \code{"glm"} as linear model for continuous outcome, \code{"glm_logit"}
#' as logistic model for binary outcome, \code{"glm_probit"} as probit model
#' for binary outcome, \code{"linpred"} as linear model for continuous
#' outcome with linear predictors specified by \code{"lp.pg"},
#' \code{"zir_logit"} as zero inflated model using logistic model to fit
#' non-zero probability, \code{"zir_probit"} as zero inflated model using
#' probit model to fit non-zero probability. Don't need to be specified if
#' \code{pg} is given.
#' @param pg Prognostic score as numeric matrix given by user. The first
#' column is potential outcome for control group and the second column is
#' potential outcome for treatment group. Don't need to be specified if
#' \code{model.pg} is given.
#' @param lp.pg Linear predictors for prognostic score as numeric vector or
#' matrix. Don't need to be specified if \code{model.pg} is not
#' \code{"linpred"}.
#' @param cov.balance A logical scalar for whether covariance balance
#' results should be shown.
#' @param varest A logical scalar for whether variance of estimator should
#' be estimated.
#' @param boots A numeric scalar for number of bootstrap relicates in
#' variance estimation. Don't need to be specified if \code{varest} is
#' \code{F}.
#' @param mc A logical scalar for whether multiple cores are used in
#' variance estimation. Don't need to be specified if \code{varest} is
#' \code{F}.
#' @param ncpus A numeric scalar for number of cores used in
#' variance estimation. Don't need to be specified if \code{varest} is
#' \code{F}.
#'
#' @return Results are put in a list:
#'   \item{est.ds}{Point estimate of ATE if matching is based on
#'   double score}
#'   \item{est.ps}{Point estimate of ATE if matching is based on
#'   propensity score}
#'   \item{est.pg}{Point estimate of ATE if matching is based on
#'   prognostic score}
#'   \item{est.x}{Point estimate of ATE if matching is based on
#'   covarites directly}
#'   \item{boot.var}{Variance of estimator estimated by bootstrap.
#'   Meaningless if \code{varest} if \code{F}.}
#'   \item{bootq1}{0.025 quantile of estimator estimated by bootstrap.
#'   Meaningless if \code{varest} if \code{F}.}
#'   \item{bootq2}{0.975 quantile of estimator estimated by bootstrap.
#'   Meaningless if \code{varest} if \code{F}.}
#'
#' @examples
#' # import lalonde data from package "Matching"
#' library(Matching)
#' data("lalonde")
#' Y = lalonde[,"re78"]
#' X = lalonde[,c("age","educ","black","hisp","married","nodegr","re74","re75","u74","u75")]
#' X = as.matrix(X)
#' A = lalonde[,"treat"]
#'
#' # linear predictors using in the algorithm
#' # take logarithm for income and standardize covariates
#' Z = X
#' Z[,"re74"] = log(Z[,"re74"] + 1)
#' Z[,"re75"] = log(Z[,"re75"] + 1)
#' Z[,"age"] = (Z[,"age"] - mean(Z[,"age"])) / sd(Z[,"age"])
#' Z[,"education"] = (Z[,"education"] - mean(Z[,"education"])) / sd(Z[,"education"])
#' Z[,"re74"] = (Z[,"re74"] - mean(Z[,"re74"])) / sd(Z[,"re74"])
#' Z[,"re75"] = (Z[,"re75"] - mean(Z[,"re75"])) / sd(Z[,"re75"])
#' Z = cbind(Z, Z[,"age"]^2, Z[,"educ"]^2, Z[,"re74"]^2, Z[,"re75"]^2)
#'
#' # estimate ATT using four matching methods
#' set.seed(1)
#' dsmatchATE(Y, X, A, method = "dsm", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#' dsmatchATE(Y, X, A, method = "ps", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#' dsmatchATE(Y, X, A, method = "pg", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#' dsmatchATE(Y, X, A, method = "cov", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#'
#' # estimate QTE using double score matching
#' p = 0.3
#' set.seed(1)
#' res <- dsmatchQTE(Y, X, A, p, method = "dsm", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T)
#' res
#' # Wald interval for QTE
#' res$est.ds + qnorm(0.025) * sqrt(res$bootvar)
#' res$est.ds - qnorm(0.025) * sqrt(res$bootvar)
#'
#' @export
dsmatchATE = function(Y, X, A, method = "dsm",
  model.ps = "other", ps = NULL, lp.ps = NULL,
  model.pg = "other", pg = NULL, lp.pg = NULL,
  cov.balance = F, varest = F, boots = 100, mc = F, ncpus = 4){

  # sort A for multiple treatments situation (although not implemented yet...)
  if (is.unsorted(A)) {
    temp <- sort(A, index.return = TRUE)
    A <- A[temp$ix]
    X <- X[temp$ix, ]
    Y <- Y[temp$ix]
    if(!is.null(ps)){
      ps <- ps[temp$ix]
    }
    if(!is.null(pg)){
      pg <- pg[temp$ix, ]
    }
    if(!is.null(lp.ps)){
      if(is.vector(lp.ps)){
        lp.ps <- lp.ps[temp$ix]
      }else{
        lp.ps <- lp.ps[temp$ix, ]
      }
    }
    if(!is.null(lp.pg)){
      if(is.vector(lp.pg)){
        lp.pg <- lp.pg[temp$ix]
      }else{
        lp.pg <- lp.pg[temp$ix, ]
      }
    }
  }

  n <- length(Y)
  loc.1 <- which(A == 1)
  loc.0 <- which(A == 0)

  if (method == "dsm"){
    # if prognostic score model is not zero inflated regression model
    # i.e. prognostic score has two dimensions
    if (!grepl("zir", model.pg)) {
      # if propensity score model is given and score is missing, then fit the model
      # if propensity score model and score are both given, then we will use the given score and throw out a warning
      # if propensity score model is not given but score is invalid, then stop
      # otherwise, directly use score from the arguments
      if (model.ps == "logit" && is.null(ps)){
        glm.out <- glm(A ~ X, family = binomial(link = "logit"))
        ps <- cbind(1, X)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
      }else if (model.ps == "probit" && is.null(ps)){
        glm.out <- glm(A ~ X, family = binomial(link = "probit"))
        ps <- cbind(1, X)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
      }else if (model.ps == "linpred" && is.null(ps)){
        # if linear predictors are valid, then fit a linear model
        if(is.lp(lp.ps, n)){
          glm.out <- glm(A ~ lp.ps, family = binomial)
          ps <- cbind(1, lp.ps)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
        }else{
          stop("invalid linear predictors for propensity score")
        }
      }else if (!(model.ps == "other") && !(is.null(ps))){
        warning("propensity score is given while a fitting model is chosen")
      }else if (is.ps(ps, n) == F){
        stop("propensity score should be given or the score is invalid")
      }

      # if prognostic score model is given and score is missing, then fit the model
      # if prognostic score model and score are both given, then we will use the given score and throw out a warning
      # if prognostic score model is not given but score is invalid, then stop
      # otherwise, directly use score from the arguments
      if (model.pg == "glm" && is.null(pg)){
        lm1.out <- lm(Y[loc.1] ~ X[loc.1, ])
        lm0.out <- lm(Y[loc.0] ~ X[loc.0, ])
        mu1 <- cbind(1, X)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
        mu0 <- cbind(1, X)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
      }else if (model.pg == "glm_logit" && is.null(pg)){
        glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "logit"))
        mu1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "logit"))
        mu0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
      }else if (model.pg == "glm_probit" && is.null(pg)){
        glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "probit"))
        mu1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "probit"))
        mu0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
      }else if (model.pg == "linpred" && is.null(pg)){
        # if linear predictors are valid, then fit a linear model
        if(is.lp(lp.pg, n)){
          if(is.vector(lp.pg)){
            lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1])
            lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0])
            mu1 <- cbind(1, lp.pg)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
            mu0 <- cbind(1, lp.pg)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
          }else{
            lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1, ])
            lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0, ])
            mu1 <- cbind(1, lp.pg)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
            mu0 <- cbind(1, lp.pg)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
          }
        }else{
          stop("invalid linear predictors for propensity score")
        }
      }else if (!(model.pg == "other") && !(is.null(pg))){
        warning("prognostic score is given while a fitting model is chosen")
      }else if (is.pg(pg, n) == F){
        stop("prognostic score should be given or the score is invalid")
      }else{
        mu0 = pg[,1]
        mu1 = pg[,2]
      }

      # combine propensity score and prognostic score together
      # to construct double score, then standardize them
      doublescore <- cbind(mu0, mu1, ps)
      dmean <- apply(doublescore, 2, mean)
      dse <- apply(doublescore, 2, sd)
      doublescore <- cbind((mu0 - dmean[1]) / dse[1],
        (mu1 - dmean[2]) / dse[2],
        (ps - dmean[3]) / dse[3])

      # apply method of seive to obtain estimates of \mu
      lm.y <- Y
      lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1,
        ps, ps ^ 2, ps * mu0, ps * mu1, ps * mu0 * mu1)
      lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
      mu1.seive <- cbind(1, lm.x)[,which(!is.na(lm.out1$coefficients))] %*% lm.out1$coefficients[which(!is.na(lm.out1$coefficients))]
      sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
      lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
      mu0.seive <- cbind(1, lm.x)[,which(!is.na(lm.out0$coefficients))] %*% lm.out0$coefficients[which(!is.na(lm.out0$coefficients))]
      sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

      mu1.ds <- mu1.seive
      mu0.ds <- mu0.seive

      trtnumber <- length(unique(A)) # number of treatment levels
      trtlevels <- unique(A) # all treatment levels
      pertrtlevelnumber <- table(A) # number of observations by treatment level

      Kiw.ds <- matrix(NA,n,1)  #Kiw is vector of number of times unit i used as a match
      Yiw.ds <- matrix(NA,n,trtnumber) #Yiw is the full imputed data set

      # use double score to match
      for (kk in 1:trtnumber) {
        thistrt <- trtlevels[kk]
        if (kk == 1) {
          fromto <- 1:pertrtlevelnumber[1]
        }
        if (kk > 1) {
          fromto <- (1:pertrtlevelnumber[kk]) + sum(pertrtlevelnumber[1:(kk - 1)])
        }
        A1 <- A != thistrt
        out1 <- Matching::Match(Y = Y, Tr = A1, X = doublescore,
          distance.tolerance = 0, ties = FALSE, Weight = 2)
        # bias correction for mu
        if (thistrt == 1) {
          bias.mu1 <- sum(mu1.ds[out1$index.control] - mu1.ds[out1$index.treated]) / n
        }
        if (thistrt == 0) {
          bias.mu0 <- sum(mu0.ds[out1$index.control] - mu0.ds[out1$index.treated]) / n
        }
        mdata1 <- out1$mdata
        Kiw.ds[fromto, 1] <- table(factor(out1$index.control, levels = fromto))
        Yiw.ds[which(A == thistrt), kk] <- Y[which(A == thistrt)]
        Yiw.ds[which(A != thistrt), kk] <- mdata1$Y[which(mdata1$Tr == 0)]
      }
      est.ds <- mean(Yiw.ds[, 2]) - mean(Yiw.ds[, 1]) - (bias.mu1 - bias.mu0)

      if (cov.balance){
        # duplicate data generated by matching
        dup.trt <- matrix(0, 1, ncol(X))
        dup.ctl <- matrix(0, 1, ncol(X))

        for (d in loc.1) {
          for (k in 1:(Kiw.ds[d] + 1)) {
            dup.trt <- rbind(dup.trt, X[d, ])
          }
        }

        for (d in loc.0) {
          for (k in 1:(Kiw.ds[d] + 1)) {
            dup.ctl <- rbind(dup.ctl, X[d, ])
          }
        }

        # remove the first row
        dup.trt <- dup.trt[-1,]
        dup.ctl <- dup.ctl[-1,]

        for (c in 1:ncol(X)) {
          cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
          res = matrix(0, 3, 2)
          colnames(res) <- c("Before Matching", "After Matching")
          rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
          res[1, 1] <- mean(X[loc.1, c])
          res[2, 1] <- mean(X[loc.0, c])
          res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
          res[1, 2] <- mean(dup.trt[, c])
          res[2, 2] <- mean(dup.ctl[, c])
          res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
          print(res)
        }
      }

      #variance estimation
      data <- cbind(Y, A, Kiw.ds, Yiw.ds, X)
      if (varest == 1) {
        if (mc == 1) {
          results <- boot::boot(data = data, statistic = estforboot_dsATE, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg,
            parallel = "multicore", ncpus = ncpus)
        }
        if (mc == 0) {
          results <- boot::boot(data = data, statistic = estforboot_dsATE, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg)
        }

        bootvar = var(results$t, na.rm = T)
        bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
        bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
      }else {
        bootvar <- bootq1 <- bootq2 <- rep(1, 2)
      }

      return(list(est.ds = est.ds, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2))

    }else{
      # if propensity score model is given and score is missing, then fit the model
      # if propensity score model and score are both given, then we will use the given score and throw out a warning
      # if propensity score model is not given but score is invalid, then stop
      # otherwise, directly use score from the arguments
      if (model.ps == "logit" && is.null(ps)){
        glm.out <- glm(A ~ X, family = binomial(link = "logit"))
        ps <- cbind(1, X)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
      }else if (model.ps == "probit" && is.null(ps)){
        glm.out <- glm(A ~ X, family = binomial(link = "probit"))
        ps <- cbind(1, X)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
      }else if (model.ps == "linpred" && is.null(ps)){
        # if linear predictors are valid, then fit a linear model
        if(is.lp(lp.ps, n)){
          glm.out <- glm(A ~ lp.ps, family = binomial)
          ps <- cbind(1, lp.ps)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
        }else{
          stop("invalid linear predictors for propensity score")
        }
      }else if (!(model.ps == "other") && !(is.null(ps))){
        warning("propensity score is given while a fitting model is chosen")
      }else if (is.ps(ps, n) == F){
        stop("propensity score should be given or the score is invalid")
      }

      # if prognostic score model is zero inflated regression model
      # i.e. prognostic score has four dimensions

      # location where outcome is not zero
      # each for whole set, treatment group and control group
      loc.r <- which(Y != 0)
      loc.1r <- intersect(loc.r, loc.1)
      loc.0r <- intersect(loc.r, loc.0)

      # change non-zero outcomes into 1
      # then run a glm model with a binomial family
      # use linear predictors as part of the prognostic score
      # calculate non-zero probability to estimate \mu
      if(model.pg == "zir_logit"){
        glm.y1b <- Y
        glm.y1b[loc.1r] <- 1
        glm.y1b <- glm.y1b[loc.1]
        glm.x1b <- X[loc.1, ]
        glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "logit"))
        p1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        pb1 <- 1 / (1 + exp(-p1))

        glm.y0b <- Y
        glm.y0b[loc.0r] <- 1
        glm.y0b <- glm.y0b[loc.0]
        glm.x0b <- X[loc.0, ]
        glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "logit"))
        p0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
        pb0 <- 1 / (1 + exp(-p0))
      }else if(model.pg == "zir_probit"){
        glm.y1b <- Y
        glm.y1b[loc.1r] <- 1
        glm.y1b <- glm.y1b[loc.1]
        glm.x1b <- X[loc.1, ]
        glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "probit"))
        p1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        pb1 <- pnorm(p1)

        glm.y0b <- Y
        glm.y0b[loc.0r] <- 1
        glm.y0b <- glm.y0b[loc.0]
        glm.x0b <- X[loc.0, ]
        glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "probit"))
        p0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
        pb0 <- pnorm(p0)
      }

      # run regression model for non-zero outcomes
      lm1.out <- lm(Y[loc.1r] ~ X[loc.1r, ])
      lm0.out <- lm(Y[loc.0r] ~ X[loc.0r, ])
      mu1 <- cbind(1, X)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
      mu0 <- cbind(1, X)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]

      # combine propensity score and prognostic score together
      # to construct double score, then standardize them
      doublescore <- cbind(p0, p1, mu0, mu1, ps)
      dmean <- apply(doublescore, 2, mean)
      dse <- apply(doublescore, 2, sd)
      doublescore <- cbind((p0 - dmean[1]) / dse[1],
        (p1 - dmean[2]) / dse[2],
        (mu0 - dmean[3]) / dse[3],
        (mu1 - dmean[4]) / dse[4],
        (ps - dmean[5]) / dse[5])

      # apply method of seive to obtain estimates of \mu
      # note that \mu has two parts: nonzero probability and nonzero value
      # we take their multiplication so that the estimator is unbiased of \mu
      lm.y <- Y
      lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1,
        p0, p1, p0 ^ 2, p1 ^ 2, p0 * p1,
        ps, ps ^ 2, ps * mu0, ps * mu1,
        p0 * mu0, p0 * mu1, p1 * mu0, p1 * mu1, p0 * ps, p1 * ps)

      lm.out1 <- lm(lm.y[loc.1r] ~ lm.x[loc.1r, ])
      mu1.seive <- cbind(1, lm.x)[,which(!is.na(lm.out1$coefficients))] %*% lm.out1$coefficients[which(!is.na(lm.out1$coefficients))]
      sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
      lm.out0 <- lm(lm.y[loc.0r] ~ lm.x[loc.0r, ])
      mu0.seive <- cbind(1, lm.x)[,which(!is.na(lm.out0$coefficients))] %*% lm.out0$coefficients[which(!is.na(lm.out0$coefficients))]
      sigsqhat0 <- mean((lm.out0$residuals) ^ 2)
      if(model.pg == "zir_logit"){
        glm.x1b <- lm.x[loc.1, ]
        glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "logit"))
        p1.seive <- cbind(1, lm.x)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        pb1.seive <- 1 / (1 + exp(-p1.seive))

        glm.x0b <- lm.x[loc.0, ]
        glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "logit"))
        p0.seive <- cbind(1, lm.x)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
        pb0.seive <- 1 / (1 + exp(-p0.seive))
      }else if(model.pg == "zir_probit"){
        glm.x1b <- lm.x[loc.1, ]
        glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "probit"))
        p1.seive <- cbind(1, lm.x)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        pb1.seive <- pnorm(p1.seive)

        glm.x0b <- lm.x[loc.0, ]
        glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "probit"))
        p0.seive <- cbind(1, lm.x)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
        pb0.seive <- pnorm(p0.seive)
      }

      mu1.ds <- mu1.seive * pb1.seive
      mu0.ds <- mu0.seive * pb0.seive

      trtnumber <- length(unique(A)) # number of treatment levels
      trtlevels <- unique(A) # all treatment levels
      pertrtlevelnumber <- table(A) # number of observations by treatment level

      Kiw.ds <- matrix(NA,n,1)  #Kiw is vector of number of times unit i used as a match
      Yiw.ds <- matrix(NA,n,trtnumber) #Yiw is the full imputed data set

      # use double score to match
      for (kk in 1:trtnumber) {
        thistrt <- trtlevels[kk]
        if (kk == 1) {
          fromto <- 1:pertrtlevelnumber[1]
        }
        if (kk > 1) {
          fromto <- (1:pertrtlevelnumber[kk]) + sum(pertrtlevelnumber[1:(kk - 1)])
        }
        A1 <- A != thistrt
        out1 <- Matching::Match(Y = Y, Tr = A1, X = doublescore,
          distance.tolerance = 0, ties = FALSE, Weight = 2)
        # bias correction for mu
        if (thistrt == 1) {
          bias.mu1 <- sum(mu1.ds[out1$index.control] - mu1.ds[out1$index.treated]) / n
        }
        if (thistrt == 0) {
          bias.mu0 <- sum(mu0.ds[out1$index.control] - mu0.ds[out1$index.treated]) / n
        }
        mdata1 <- out1$mdata
        Kiw.ds[fromto, 1] <- table(factor(out1$index.control, levels = fromto))
        Yiw.ds[which(A == thistrt), kk] <- Y[which(A == thistrt)]
        Yiw.ds[which(A != thistrt), kk] <- mdata1$Y[which(mdata1$Tr == 0)]
      }
      est.ds <- mean(Yiw.ds[, 2]) - mean(Yiw.ds[, 1]) - (bias.mu1 - bias.mu0)

      if (cov.balance){
        # duplicate data generated by matching
        dup.trt <- matrix(0, 1, ncol(X))
        dup.ctl <- matrix(0, 1, ncol(X))

        for (d in loc.1) {
          for (k in 1:(Kiw.ds[d] + 1)) {
            dup.trt <- rbind(dup.trt, X[d, ])
          }
        }

        for (d in loc.0) {
          for (k in 1:(Kiw.ds[d] + 1)) {
            dup.ctl <- rbind(dup.ctl, X[d, ])
          }
        }

        # remove the first row
        dup.trt <- dup.trt[-1,]
        dup.ctl <- dup.ctl[-1,]

        for (c in 1:ncol(X)) {
          cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
          res = matrix(0, 3, 2)
          colnames(res) <- c("Before Matching", "After Matching")
          rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
          res[1, 1] <- mean(X[loc.1, c])
          res[2, 1] <- mean(X[loc.0, c])
          res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
          res[1, 2] <- mean(dup.trt[, c])
          res[2, 2] <- mean(dup.ctl[, c])
          res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
          print(res)
        }
      }

      #variance estimation
      data <- cbind(Y, A, Kiw.ds, Yiw.ds, X)
      if (varest == 1) {
        if (mc == 1) {
          results <- boot::boot(data = data, statistic = estforboot_dsATE, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg,
            parallel = "multicore", ncpus = ncpus)
        }
        if (mc == 0) {
          results <- boot::boot(data = data, statistic = estforboot_dsATE, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg)
        }

        bootvar = var(results$t, na.rm = T)
        bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
        bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
      }else {
        bootvar <- bootq1 <- bootq2 <- rep(1, 2)
      }

      return(list(est.ds = est.ds, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2))
    }


  }else if(method == "ps"){
    # if propensity score model is given and score is missing, then fit the model
    # if propensity score model and score are both given, then we will use the given score and throw out a warning
    # if propensity score model is not given but score is invalid, then stop
    # otherwise, directly use score from the arguments
    if (model.ps == "logit" && is.null(ps)){
      glm.out <- glm(A ~ X, family = binomial(link = "logit"))
      ps <- cbind(1, X)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
    }else if (model.ps == "probit" && is.null(ps)){
      glm.out <- glm(A ~ X, family = binomial(link = "probit"))
      ps <- cbind(1, X)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
    }else if (model.ps == "linpred" && is.null(ps)){
      # if linear predictors are valid, then fit a linear model
      if(is.lp(lp.ps, n)){
        glm.out <- glm(A ~ lp.ps, family = binomial)
        ps <- cbind(1, lp.ps)[,which(!is.na(glm.out$coefficients))] %*% glm.out$coefficients[which(!is.na(glm.out$coefficients))]
      }else{
        stop("invalid linear predictors for propensity score")
      }
    }else if (!(model.ps == "other") && !(is.null(ps))){
      warning("propensity score is given while a fitting model is chosen")
    }else if (is.ps(ps, n) == F){
      stop("propensity score should be given or the score is invalid")
    }

    trtnumber <- length(unique(A)) # number of treatment levels
    trtlevels <- unique(A) # all treatment levels
    pertrtlevelnumber <- table(A) # number of observations by treatment level

    Kiw.ps <- matrix(NA,n,1)  #Kiw is vector of number of times unit i used as a match
    Yiw.ps <- matrix(NA,n,trtnumber) #Yiw is the full imputed data set
    for (kk in 1:trtnumber) {
      thistrt <- trtlevels[kk]
      if (kk == 1) {
        fromto <- 1:pertrtlevelnumber[1]
      }
      if (kk > 1) {
        fromto <- (1:pertrtlevelnumber[kk]) + sum(pertrtlevelnumber[1:(kk - 1)])
      }
      A1 <- A != thistrt
      out1 <- Matching::Match(Y = Y, Tr = A1, X = ps,
        distance.tolerance = 0, ties = FALSE, Weight = 2)
      mdata1 <- out1$mdata
      Kiw.ps[fromto, 1] <- table(factor(out1$index.control, levels = fromto))
      Yiw.ps[which(A == thistrt), kk] <- Y[which(A == thistrt)]
      Yiw.ps[which(A != thistrt), kk] <- mdata1$Y[which(mdata1$Tr == 0)]
    }
    est.ps <- mean(Yiw.ps[, 2]) - mean(Yiw.ps[, 1])

    if (cov.balance){
      # duplicate data generated by matching
      dup.trt <- matrix(0, 1, ncol(X))
      dup.ctl <- matrix(0, 1, ncol(X))

      for (d in loc.1) {
        for (k in 1:(Kiw.ps[d] + 1)) {
          dup.trt <- rbind(dup.trt, X[d, ])
        }
      }

      for (d in loc.0) {
        for (k in 1:(Kiw.ps[d] + 1)) {
          dup.ctl <- rbind(dup.ctl, X[d, ])
        }
      }

      # remove the first row
      dup.trt <- dup.trt[-1,]
      dup.ctl <- dup.ctl[-1,]

      for (c in 1:ncol(X)) {
        cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
        res = matrix(0, 3, 2)
        colnames(res) <- c("Before Matching", "After Matching")
        rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
        res[1, 1] <- mean(X[loc.1, c])
        res[2, 1] <- mean(X[loc.0, c])
        res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
        res[1, 2] <- mean(dup.trt[, c])
        res[2, 2] <- mean(dup.ctl[, c])
        res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
        print(res)
      }
    }

    #variance estimation
    data <- cbind(Y, A, Kiw.ps, Yiw.ps, X)
    if (varest == 1) {
      if (mc == 1) {
        results <- boot::boot(data = data, statistic = estforboot_psATE, R = boots,
          model.ps = model.ps, ps = ps, lp.ps = lp.ps,
          parallel = "multicore", ncpus = ncpus)
      }
      if (mc == 0) {
        results <- boot::boot(data = data, statistic = estforboot_psATE, R = boots,
          model.ps = model.ps, ps = ps, lp.ps = lp.ps)
      }

      bootvar = var(results$t, na.rm = T)
      bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
      bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
    }else {
      bootvar <- bootq1 <- bootq2 <- rep(1, 2)
    }

    return(list(est.ps = est.ps, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2))

  }else if(method == "pg"){
    # if prognostic score model is given and score is missing, then fit the model
    # if prognostic score model and score are both given, then we will use the given score and throw out a warning
    # if prognostic score model is not given but score is invalid, then stop
    # otherwise, directly use score from the arguments
    if (model.pg == "glm" && is.null(pg)){
      lm1.out <- lm(Y[loc.1] ~ X[loc.1, ])
      lm0.out <- lm(Y[loc.0] ~ X[loc.0, ])
      mu1 <- cbind(1, X)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
      mu0 <- cbind(1, X)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
    }else if (model.pg == "glm_logit" && is.null(pg)){
      glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "logit"))
      mu1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
      glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "logit"))
      mu0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
    }else if (model.pg == "glm_probit" && is.null(pg)){
      glm.out1 <- glm(Y[loc.1] ~ X[loc.1, ], family = binomial(link = "probit"))
      mu1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
      glm.out0 <- glm(Y[loc.0] ~ X[loc.0, ], family = binomial(link = "probit"))
      mu0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
    }else if (model.pg == "linpred" && is.null(pg)){
      # if linear predictors are valid, then fit a linear model
      if(is.lp(lp.pg, n)){
        if(is.vector(lp.pg)){
          lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1])
          lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0])
          mu1 <- cbind(1, lp.pg)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
          mu0 <- cbind(1, lp.pg)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
        }else{
          lm1.out <- lm(Y[loc.1] ~ lp.pg[loc.1, ])
          lm0.out <- lm(Y[loc.0] ~ lp.pg[loc.0, ])
          mu1 <- cbind(1, lp.pg)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
          mu0 <- cbind(1, lp.pg)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]
        }
      }else{
        stop("invalid linear predictors for propensity score")
      }
    }else if (!(model.pg == "other") && !(is.null(pg))){
      warning("prognostic score is given while a fitting model is chosen")
    }else if (is.pg(pg, n) == F){
      stop("prognostic score should be given or the score is invalid")
    }else{
      mu0 = pg[,1]
      mu1 = pg[,2]
    }

    lm.y <- Y
    lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1)
    lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
    mu1.seive <- cbind(1, lm.x)[,which(!is.na(lm.out1$coefficients))] %*% lm.out1$coefficients[which(!is.na(lm.out1$coefficients))]
    sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
    lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
    mu0.seive <- cbind(1, lm.x)[,which(!is.na(lm.out0$coefficients))] %*% lm.out0$coefficients[which(!is.na(lm.out0$coefficients))]
    sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

    mu1.pg <- mu1.seive
    mu0.pg <- mu0.seive

    coeff1.pg <- lm.out1$coefficients
    coeff0.pg <- lm.out0$coefficients

    trtnumber <- length(unique(A)) # number of treatment levels
    trtlevels <- unique(A) # all treatment levels
    pertrtlevelnumber <- table(A) # number of observations by treatment level

    Kiw.pg <- matrix(NA, n, 1)  #Kiw is vector of number of times unit i used as a match
    Yiw.pg <- matrix(NA, n, trtnumber) #Yiw is the full imputed data set
    for (kk in 1:trtnumber) {
      thistrt <- trtlevels[kk]
      if (kk == 1) {
        fromto <- 1:pertrtlevelnumber[1]
      }
      if (kk > 1) {
        fromto <- (1:pertrtlevelnumber[kk]) + sum(pertrtlevelnumber[1:(kk - 1)])
      }
      A1 <- A != thistrt
      out1 <- Matching::Match(Y = Y, Tr = A1, X = cbind(mu0, mu1),
        distance.tolerance = 0, ties = FALSE, Weight = 2)
      mdata1 <- out1$mdata
      if (thistrt == 1) {
        bias.mu1 <- sum(mu1.pg[out1$index.control] - mu1.pg[out1$index.treated]) / n
      }
      if (thistrt == 0) {
        bias.mu0 <- sum(mu0.pg[out1$index.control] - mu0.pg[out1$index.treated]) / n
      }

      Kiw.pg[fromto, 1] <- table(factor(out1$index.control, levels = fromto))
      Yiw.pg[which(A == thistrt), kk] <- Y[which(A == thistrt)]
      Yiw.pg[which(A != thistrt), kk] <- mdata1$Y[which(mdata1$Tr == 0)]
    }

    est.pg <- mean(Yiw.pg[, 2]) - mean(Yiw.pg[, 1]) - (bias.mu1 - bias.mu0)

    if (cov.balance){
      # duplicate data generated by matching
      dup.trt <- matrix(0, 1, ncol(X))
      dup.ctl <- matrix(0, 1, ncol(X))

      for (d in loc.1) {
        for (k in 1:(Kiw.pg[d] + 1)) {
          dup.trt <- rbind(dup.trt, X[d, ])
        }
      }

      for (d in loc.0) {
        for (k in 1:(Kiw.pg[d] + 1)) {
          dup.ctl <- rbind(dup.ctl, X[d, ])
        }
      }

      # remove the first row
      dup.trt <- dup.trt[-1,]
      dup.ctl <- dup.ctl[-1,]

      for (c in 1:ncol(X)) {
        cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
        res = matrix(0, 3, 2)
        colnames(res) <- c("Before Matching", "After Matching")
        rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
        res[1, 1] <- mean(X[loc.1, c])
        res[2, 1] <- mean(X[loc.0, c])
        res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
        res[1, 2] <- mean(dup.trt[, c])
        res[2, 2] <- mean(dup.ctl[, c])
        res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
        print(res)
      }
    }

    #variance estimation
    data <- cbind(Y, A, Kiw.pg, Yiw.pg, X)
    if (varest == 1) {
      if (mc == 1) {
        results <- boot::boot(data = data, statistic = estforboot_pgATE, R = boots,
          model.pg = model.pg, pg = pg, lp.pg = lp.pg,
          parallel = "multicore", ncpus = ncpus)
      }
      if (mc == 0) {
        results <- boot::boot(data = data, statistic = estforboot_pgATE, R = boots,
          model.pg = model.pg, pg = pg, lp.pg = lp.pg)
      }

      bootvar = var(results$t, na.rm = T)
      bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
      bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
    }else {
      bootvar <- bootq1 <- bootq2 <- rep(1, 2)
    }

    return(list(est.pg = est.pg, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2))

  }else if(method == "cov"){
    lm.y <- Y
    lm.x <- cbind(X, t(apply(X, 1, combn, 2, prod)), X ^ 2, X ^ 3)
    # lm.x <- cbind(X)

    lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
    mu1.seive <- cbind(1, lm.x)[,which(!is.na(lm.out1$coefficients))] %*% lm.out1$coefficients[which(!is.na(lm.out1$coefficients))]
    sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
    lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
    mu0.seive <- cbind(1, lm.x)[,which(!is.na(lm.out0$coefficients))] %*% lm.out0$coefficients[which(!is.na(lm.out0$coefficients))]
    sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

    mu1.x <- mu1.seive
    mu0.x <- mu0.seive

    trtnumber <- length(unique(A)) # number of treatment levels
    trtlevels <- unique(A) # all treatment levels
    pertrtlevelnumber <- table(A) # number of observations by treatment level

    Kiw.x <- matrix(NA, n, 1)  # Kiw is vector of number of times unit i used as a match
    Yiw.x <- matrix(NA, n, trtnumber) # Yiw is the full imputed data set

    for (kk in 1:trtnumber) {
      thistrt <- trtlevels[kk]
      if (kk == 1) {
        fromto <- 1:pertrtlevelnumber[1]
      }
      if (kk > 1) {
        fromto <- (1:pertrtlevelnumber[kk]) + sum(pertrtlevelnumber[1:(kk - 1)])
      }
      A1 <- A != thistrt
      out1 <- Matching::Match(Y = Y, Tr = A1, X = X,
        distance.tolerance = 0, ties = FALSE, Weight = 2)
      if (thistrt == 1) {
        bias.mu1 <- sum(mu1.x[out1$index.control] - mu1.x[out1$index.treated]) /n
      }
      if (thistrt == 0) {
        bias.mu0 <- sum(mu0.x[out1$index.control] - mu0.x[out1$index.treated]) /n
      }
      mdata1 <- out1$mdata
      Kiw.x[fromto, 1] <- table(factor(out1$index.control, levels = fromto))
      Yiw.x[which(A == thistrt), kk] <- Y[which(A == thistrt)]
      Yiw.x[which(A != thistrt), kk] <- mdata1$Y[which(mdata1$Tr == 0)]
    }

    est.x <- mean(Yiw.x[, 2]) - mean(Yiw.x[, 1]) - (bias.mu1 - bias.mu0)

    if (cov.balance){
      # duplicate data generated by matching
      dup.trt <- matrix(0, 1, ncol(X))
      dup.ctl <- matrix(0, 1, ncol(X))

      for (d in loc.1) {
        for (k in 1:(Kiw.x[d] + 1)) {
          dup.trt <- rbind(dup.trt, X[d, ])
        }
      }

      for (d in loc.0) {
        for (k in 1:(Kiw.x[d] + 1)) {
          dup.ctl <- rbind(dup.ctl, X[d, ])
        }
      }

      # remove the first row
      dup.trt <- dup.trt[-1,]
      dup.ctl <- dup.ctl[-1,]

      for (c in 1:ncol(X)) {
        cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
        res = matrix(0, 3, 2)
        colnames(res) <- c("Before Matching", "After Matching")
        rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
        res[1, 1] <- mean(X[loc.1, c])
        res[2, 1] <- mean(X[loc.0, c])
        res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
        res[1, 2] <- mean(dup.trt[, c])
        res[2, 2] <- mean(dup.ctl[, c])
        res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
        print(res)
      }
    }

    #variance estimation
    data <- cbind(Y, A, Kiw.x, Yiw.x, X)
    if (varest == 1) {
      if (mc == 1) {
        results <- boot::boot(data = data, statistic = estforboot_covATE, R = boots,
          parallel = "multicore", ncpus = ncpus)
      }
      if (mc == 0) {
        results <- boot::boot(data = data, statistic = estforboot_covATE, R = boots)
      }

      bootvar = var(results$t, na.rm = T)
      bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
      bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
    }else {
      bootvar <- bootq1 <- bootq2 <- rep(1, 2)
    }

    return(list(est.x = est.x, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2))

  }else{
    stop("invalid method")
  }

}











