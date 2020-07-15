#' Double Score Matching Estimator for Average Treatment Effect for the Treated.
#'
#' \code{dsmatchATT} applys matching algorithm to estimate average
#' treatment effect for the treated based on propensity score and prognostic score.
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
#' model, which fits a logistic model for the probability to be zero,
#' and a regression model for the non-zero values.
#'
#' Some parameters in \code{Match} in \code{Matching} has already been
#' set as parameters in the function, such as \code{caliper} and
#' \code{replace}. Addtional parameters for \code{Match} function can
#' also be assigned by \code{...} except that \code{tie} and
#' \code{Weight} has already been specified in the function.
#'
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
#' @param caliper A scalar or vector denoting the caliper(s) which should
#' be used when matching. A caliper is the distance which is acceptable
#' for any match. Observations which are outside of the caliper are dropped.
#' If a scalar caliper is provided, this caliper is used for all covariates
#' in X. If a vector of calipers is provided, a caliper value should be
#' provided for each covariate in X. See function \code{Match} in
#' \code{Matching} for more details.
#' @param replace A logical flag for whether matching should be done with
#' replacement. Note that if FALSE, the order of matches generally matters.
#' See function \code{Match} in \code{Matching} for more details.
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
#' @param ... Additional parameters for \code{Match} in \code{Matching}.
#'
#' @return Results are put in a list:
#'   \item{est.ds}{Point estimate of ATT if matching is based on
#'   double score}
#'   \item{est.ps}{Point estimate of ATT if matching is based on
#'   propensity score}
#'   \item{est.pg}{Point estimate of ATT if matching is based on
#'   prognostic score}
#'   \item{est.x}{Point estimate of ATT if matching is based on
#'   covarites directly}
#'   \item{boot.var}{Variance of estimator estimated by bootstrap.
#'   Meaningless if \code{varest} if \code{F}.}
#'   \item{bootq1}{0.025 quantile of estimator estimated by bootstrap.
#'   Meaningless if \code{varest} if \code{F}.}
#'   \item{bootq2}{0.975 quantile of estimator estimated by bootstrap.
#'   Meaningless if \code{varest} if \code{F}.}
#'   \item{cov.bal}{standard difference in mean for all covariates}
#'   \item{matching.detail}{returned object from function \code{Match}
#'   in package \code{Matching}}
#'   \item{matching.rate}{matching rate due to caliper or replacement}
#'
#' @examples
#' # import lalonde data from package "lalonde"
#' library(lalonde)
#' nsw <- lalonde::nsw
#' cps3 <- lalonde::cps_controls3
#'
#' # combine datasets
#' nsw <- nsw[,-1]
#' cps3 <- cps3[,c(2,3,4,5,6,7,8,10,11)]
#' lalonde <- rbind(nsw, cps3)
#'
#' # preprocessing of data
#' Y = lalonde[,"re78"]
#' Y = as.matrix(Y)
#' Y = as.vector(Y)
#' X = lalonde[,c("age","education","black","hispanic","married","nodegree","re75")]
#' X = as.matrix(X)
#' A = lalonde[,"treat"]
#' A = as.matrix(A)
#' A = as.vector(A)
#'
#' # linear predictors using in the algorithm
#' # take logarithm for income and standardize covariates
#' Z = X
#' Z[,"re75"] = log(Z[,"re75"] + 1)
#' Z[,"age"] = (Z[,"age"] - mean(Z[,"age"])) / sd(Z[,"age"])
#' Z[,"education"] = (Z[,"education"] - mean(Z[,"education"])) / sd(Z[,"education"])
#' Z[,"re75"] = (Z[,"re75"] - mean(Z[,"re75"])) / sd(Z[,"re75"])
#' Z = cbind(X, Z[,"age"]^2, Z[,"education"]^2, Z[,"re75"]^2)
#'
#' # estimate ATT using four matching methods
#' set.seed(1)
#' dsmatchATT(Y, X, A, method = "dsm", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#' dsmatchATT(Y, X, A, method = "ps", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#' dsmatchATT(Y, X, A, method = "pg", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#' dsmatchATT(Y, X, A, method = "cov", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T, cov.balance = T)
#'
#' # estimate QTT using double score matching
#' p = 0.3
#' set.seed(1)
#' res <- dsmatchQTT(Y, X, A, p, method = "dsm", model.ps = "linpred", lp.ps = Z, model.pg = "linpred", lp.pg = Z, varest = T)
#' res
#' # Wald interval for QTT
#' res$est.ds + qnorm(0.025) * sqrt(res$bootvar)
#' res$est.ds - qnorm(0.025) * sqrt(res$bootvar)
#'
#' @export
dsmatchATT = function(Y, X, A, method = "dsm",
  model.ps = "other", ps = NULL, lp.ps = NULL,
  model.pg = "other", pg = NULL, lp.pg = NULL,
  caliper = NULL, replace = T, cov.balance = F,
  varest = F, boots = 100, mc = F, ncpus = 4, ...){

  if(!is.matrix(X)){
    X = as.matrix(X)
  }

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
  n1 <- length(loc.1)

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
        # mu0 = pg[,1]
        # mu1 = pg[,2]
        mu0 = pg
      }

      # combine propensity score and prognostic score together
      # to construct double score, then standardize them
      # doublescore <- cbind(mu0, mu1, ps)
      # dmean <- apply(doublescore, 2, mean)
      # dse <- apply(doublescore, 2, sd)
      # doublescore <- cbind((mu0 - dmean[1]) / dse[1],
      #   (mu1 - dmean[2]) / dse[2],
      #   (ps - dmean[3]) / dse[3])
      doublescore <- cbind(mu0, ps)
      dmean <- apply(doublescore, 2, mean)
      dse <- apply(doublescore, 2, sd)
      doublescore <- cbind((mu0 - dmean[1]) / dse[1],
        (ps - dmean[2]) / dse[2])

      # apply method of seive to obtain estimates of \mu
      lm.y <- Y
      # lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1,
      #   ps, ps ^ 2, ps * mu0, ps * mu1, ps * mu0 * mu1)
      lm.x <- cbind(mu0, mu0 ^ 2, ps, ps ^ 2, ps * mu0)
      lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
      mu1.seive <- cbind(1, lm.x)[,which(!is.na(lm.out1$coefficients))] %*% lm.out1$coefficients[which(!is.na(lm.out1$coefficients))]
      sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
      lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
      mu0.seive <- cbind(1, lm.x)[,which(!is.na(lm.out0$coefficients))] %*% lm.out0$coefficients[which(!is.na(lm.out0$coefficients))]
      sigsqhat0 <- mean((lm.out0$residuals) ^ 2)

      mu1.ds <- mu1.seive
      mu0.ds <- mu0.seive

      Kiw.ds <- rep(0,n)  #Kiw is vector of number of times unit i used as a match

      if(is.null(caliper)){
        thistrt <- 0
        A1 <- A != thistrt
        out1 <- Matching::Match(Y = Y, Tr = A1, X = doublescore,
                                distance.tolerance = 0, ties = FALSE, Weight = 2,
                                caliper = NULL, replace = replace, ...)
        matching.rate = 1 - out1$ndrops / n1
        # bias correction for mu
        bias.mu0 <- sum(mu0.ds[out1$index.control] - mu0.ds[out1$index.treated]) / n1
        mdata1 <- out1$mdata
        Kiw.ds[loc.0] <- table(factor(out1$index.control, levels = loc.0))
        est.ds <- mean(Y[loc.1]) - mean(mdata1$Y[which(mdata1$Tr == 0)]) + bias.mu0

        # covariate balance checking
        # duplicate data generated by matching
        dup.trt <- X[loc.1, ]
        dup.ctl <- matrix(0, 1, ncol(X))

        for (d in loc.0) {
          if(Kiw.ds[d] > 0){
            for (k in 1:Kiw.ds[d]) {
              dup.ctl <- rbind(dup.ctl, X[d, ])
            }
          }
        }

        # remove the first row
        dup.ctl <- dup.ctl[-1,]

        # record standard difference in mean for all covariates
        sdm.cov = rep(0, ncol(X))

        for (c in 1:ncol(X)) {
          res = matrix(0, 3, 2)
          colnames(res) <- c("Before Matching", "After Matching")
          rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
          res[1, 1] <- mean(X[loc.1, c])
          res[2, 1] <- mean(X[loc.0, c])
          res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
          res[1, 2] <- mean(dup.trt[, c])
          res[2, 2] <- mean(dup.ctl[, c])
          res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
          if (cov.balance){
            cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
            print(res)
          }
          sdm.cov[c] = res[3, 2]
        }
      }else{
        # thistrt <- 0
        # A1 <- A != thistrt
        # out1 <- Matching::Match(Y = Y, Tr = A1, X = doublescore,
        #                         distance.tolerance = 0, ties = FALSE, Weight = 2,
        #                         caliper = NULL, replace = replace, ...)
        #
        # # screening based on caliper
        # # at least one score need to fall in the caliper
        # index.control.keep = c()
        # index.treated.keep = c()
        # for (i in 1:length(out1$index.treated)) {
        #   i0 = out1$index.control[i]
        #   i1 = out1$index.treated[i]
        #   ismuin = abs(doublescore[i0,1] - doublescore[i1, 1]) <= caliper
        #   ispsin = abs(doublescore[i0,2] - doublescore[i1, 2]) <= caliper
        #   if(ismuin | ispsin){
        #     index.control.keep = c(index.control.keep, i0)
        #     index.treated.keep = c(index.treated.keep, i1)
        #   }
        # }
        # matching.rate = length(index.control.keep) / n1
        #
        # # bias correction for mu
        # bias.mu0 <- sum(mu0.ds[index.control.keep] - mu0.ds[index.treated.keep]) / length(index.control.keep)
        # Kiw.ds[loc.0] <- table(factor(index.control.keep, levels = loc.0))
        # est.ds <- mean(Y[index.treated.keep]) - mean(Y[index.control.keep]) + bias.mu0

        out <- dsmatchATT_caliper(Y = Y, A = A, X = doublescore, caliper = caliper, replace = replace)
        matching.rate = 1 - out$ndrops / n1
        bias.mu0 <- sum(mu0.ds[out$index.control] - mu0.ds[out$index.treated]) / length(out$index.control)
        Kiw.ds[loc.0] <- table(factor(out$index.control, levels = loc.0))
        est.ds <- mean(Y[out$index.treated]) - mean(Y[out$index.control]) + bias.mu0

        # covariate balance checking
        # duplicate data generated by matching
        # dup.trt <- X[index.treated.keep, ]
        dup.trt <- X[out$index.treated, ]
        dup.ctl <- matrix(0, 1, ncol(X))

        for (d in loc.0) {
          if(Kiw.ds[d] > 0){
            for (k in 1:Kiw.ds[d]) {
              dup.ctl <- rbind(dup.ctl, X[d, ])
            }
          }
        }

        # remove the first row
        dup.ctl <- dup.ctl[-1,]

        # record standard difference in mean for all covariates
        sdm.cov = rep(0, ncol(X))

        for (c in 1:ncol(X)) {
          res = matrix(0, 3, 2)
          colnames(res) <- c("Before Matching", "After Matching")
          rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
          res[1, 1] <- mean(X[loc.1, c])
          res[2, 1] <- mean(X[loc.0, c])
          res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
          res[1, 2] <- mean(dup.trt[, c])
          res[2, 2] <- mean(dup.ctl[, c])
          res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
          if (cov.balance){
            cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
            print(res)
          }
          sdm.cov[c] = res[3, 2]
        }
      }

      #variance estimation
      data <- cbind(Y, A, Kiw.ds, X)
      if (varest == 1) {
        if (mc == 1) {
          results <- boot::boot(data = data, statistic = estforboot_dsATT, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg,
            parallel = "multicore", ncpus = ncpus)
        }
        if (mc == 0) {
          results <- boot::boot(data = data, statistic = estforboot_dsATT, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg)
        }

        bootvar = var(results$t, na.rm = T)
        bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
        bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
      }else {
        bootvar <- bootq1 <- bootq2 <- NA
      }

      return(list(est.ds = est.ds, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2, cov.bal = sdm.cov, matching.detail = out1, matching.rate = matching.rate))

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
      # i.e. prognostic score has two dimensions

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
        # glm.y1b <- Y
        # glm.y1b[loc.1r] <- 1
        # glm.y1b <- glm.y1b[loc.1]
        # glm.x1b <- X[loc.1, ]
        # glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "logit"))
        # p1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        # pb1 <- 1 / (1 + exp(-p1))

        glm.y0b <- Y
        glm.y0b[loc.0r] <- 1
        glm.y0b <- glm.y0b[loc.0]
        glm.x0b <- X[loc.0, ]
        glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "logit"))
        p0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
        pb0 <- 1 / (1 + exp(-p0))
      }else if(model.pg == "zir_probit"){
        # glm.y1b <- Y
        # glm.y1b[loc.1r] <- 1
        # glm.y1b <- glm.y1b[loc.1]
        # glm.x1b <- X[loc.1, ]
        # glm.out1 <- glm(glm.y1b ~ glm.x1b, family = binomial(link = "probit"))
        # p1 <- cbind(1, X)[,which(!is.na(glm.out1$coefficients))] %*% glm.out1$coefficients[which(!is.na(glm.out1$coefficients))]
        # pb1 <- pnorm(p1)

        glm.y0b <- Y
        glm.y0b[loc.0r] <- 1
        glm.y0b <- glm.y0b[loc.0]
        glm.x0b <- X[loc.0, ]
        glm.out0 <- glm(glm.y0b ~ glm.x0b, family = binomial(link = "probit"))
        p0 <- cbind(1, X)[,which(!is.na(glm.out0$coefficients))] %*% glm.out0$coefficients[which(!is.na(glm.out0$coefficients))]
        pb0 <- pnorm(p0)
      }

      # run regression model for non-zero outcomes
      # lm1.out <- lm(Y[loc.1r] ~ X[loc.1r, ])
      lm0.out <- lm(Y[loc.0r] ~ X[loc.0r, ])
      # mu1 <- cbind(1, X)[,which(!is.na(lm1.out$coefficients))] %*% lm1.out$coefficients[which(!is.na(lm1.out$coefficients))]
      mu0 <- cbind(1, X)[,which(!is.na(lm0.out$coefficients))] %*% lm0.out$coefficients[which(!is.na(lm0.out$coefficients))]

      # combine propensity score and prognostic score together
      # to construct double score, then standardize them
      # doublescore <- cbind(p0, p1, mu0, mu1, ps)
      # dmean <- apply(doublescore, 2, mean)
      # dse <- apply(doublescore, 2, sd)
      # doublescore <- cbind((p0 - dmean[1]) / dse[1],
      #   (p1 - dmean[2]) / dse[2],
      #   (mu0 - dmean[3]) / dse[3],
      #   (mu1 - dmean[4]) / dse[4],
      #   (ps - dmean[5]) / dse[5])
      doublescore <- cbind(p0, mu0, ps)
      dmean <- apply(doublescore, 2, mean)
      dse <- apply(doublescore, 2, sd)
      doublescore <- cbind((p0 - dmean[1]) / dse[1],
                           (mu0 - dmean[2]) / dse[2],
                           (ps - dmean[3]) / dse[3])

      # apply method of seive to obtain estimates of \mu
      # note that \mu has two parts: nonzero probability and nonzero value
      # we take their multiplication so that the estimator is unbiased of \mu
      lm.y <- Y
      # lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1,
      #   p0, p1, p0 ^ 2, p1 ^ 2, p0 * p1,
      #   ps, ps ^ 2, ps * mu0, ps * mu1,
      #   p0 * mu0, p0 * mu1, p1 * mu0, p1 * mu1, p0 * ps, p1 * ps)
      lm.x <- cbind(mu0, mu0 ^ 2, p0, p0 ^ 2,
                    ps, ps ^ 2, ps * mu0,
                    p0 * mu0, p0 * ps)

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

      Kiw.ds <- rep(0,n)  #Kiw is vector of number of times unit i used as a match
      thistrt <- 0
      A1 <- A != thistrt
      out1 <- Matching::Match(Y = Y, Tr = A1, X = doublescore,
        distance.tolerance = 0, ties = FALSE, Weight = 2,
        caliper = caliper, replace = replace, ...)
      # bias correction for mu
      bias.mu0 <- sum(mu0.ds[out1$index.control] - mu0.ds[out1$index.treated]) / n1
      mdata1 <- out1$mdata
      Kiw.ds[loc.0] <- table(factor(out1$index.control, levels = loc.0))
      est.ds <- mean(Y[loc.1]) - mean(mdata1$Y[which(mdata1$Tr == 0)]) + bias.mu0

      # covariate balance checking
      # duplicate data generated by matching
      dup.trt <- X[loc.1, ]
      dup.ctl <- matrix(0, 1, ncol(X))

      for (d in loc.0) {
        if(Kiw.ds[d] > 0){
          for (k in 1:Kiw.ds[d]) {
            dup.ctl <- rbind(dup.ctl, X[d, ])
          }
        }
      }

      # remove the first row
      dup.ctl <- dup.ctl[-1,]

      # record standard difference in mean for all covariates
      sdm.cov = rep(0, ncol(X))

      for (c in 1:ncol(X)) {
        res = matrix(0, 3, 2)
        colnames(res) <- c("Before Matching", "After Matching")
        rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
        res[1, 1] <- mean(X[loc.1, c])
        res[2, 1] <- mean(X[loc.0, c])
        res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
        res[1, 2] <- mean(dup.trt[, c])
        res[2, 2] <- mean(dup.ctl[, c])
        res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
        if (cov.balance){
          cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
          print(res)
        }
        sdm.cov[c] = res[3, 2]
      }

      #variance estimation
      data <- cbind(Y, A, Kiw.ds, X)
      if (varest == 1) {
        if (mc == 1) {
          results <- boot::boot(data = data, statistic = estforboot_dsATT, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg,
            parallel = "multicore", ncpus = ncpus)
        }
        if (mc == 0) {
          results <- boot::boot(data = data, statistic = estforboot_dsATT, R = boots,
            model.ps = model.ps, ps = ps, lp.ps = lp.ps,
            model.pg = model.pg, pg = pg, lp.pg = lp.pg)
        }

        bootvar = var(results$t, na.rm = T)
        bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
        bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
      }else {
        bootvar <- bootq1 <- bootq2 <- NA
      }

      return(list(est.ds = est.ds, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2, cov.bal = sdm.cov, matching.detail = out1))
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

    Kiw.ps <- rep(0,n)  #Kiw is vector of number of times unit i used as a match
    thistrt <- 0
    A1 <- A != thistrt
    out1 <- Matching::Match(Y = Y, Tr = A1, X = ps,
      distance.tolerance = 0, ties = FALSE, Weight = 2,
      caliper = caliper, replace = replace, ...)
    # # bias correction for mu
    # bias.mu0 <- sum(mu0.ds[out1$index.control] - mu0.ds[out1$index.treated]) / n1
    mdata1 <- out1$mdata
    Kiw.ps[loc.0] <- table(factor(out1$index.control, levels = loc.0))
    est.ps <- mean(Y[loc.1]) - mean(mdata1$Y[which(mdata1$Tr == 0)])
    matching.rate = 1 - out1$ndrops / n1

    # covariate balance checking
    # duplicate data generated by matching
    dup.trt <- X[loc.1, ]
    dup.ctl <- matrix(0, 1, ncol(X))

    for (d in loc.0) {
      if(Kiw.ps[d] > 0){
        for (k in 1:Kiw.ps[d]) {
          dup.ctl <- rbind(dup.ctl, X[d, ])
        }
      }
    }

    # remove the first row
    dup.ctl <- dup.ctl[-1,]

    # record standard difference in mean for all covariates
    sdm.cov = rep(0, ncol(X))

    for (c in 1:ncol(X)) {
      res = matrix(0, 3, 2)
      colnames(res) <- c("Before Matching", "After Matching")
      rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
      res[1, 1] <- mean(X[loc.1, c])
      res[2, 1] <- mean(X[loc.0, c])
      res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
      res[1, 2] <- mean(dup.trt[, c])
      res[2, 2] <- mean(dup.ctl[, c])
      res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
      if (cov.balance){
        cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
        print(res)
      }
      sdm.cov[c] = res[3, 2]
    }

    #variance estimation
    data <- cbind(Y, A, Kiw.ps, X)
    if (varest == 1) {
      if (mc == 1) {
        results <- boot::boot(data = data, statistic = estforboot_psATT, R = boots,
          model.ps = model.ps, ps = ps, lp.ps = lp.ps,
          parallel = "multicore", ncpus = ncpus)
      }
      if (mc == 0) {
        results <- boot::boot(data = data, statistic = estforboot_psATT, R = boots,
          model.ps = model.ps, ps = ps, lp.ps = lp.ps)
      }

      bootvar = var(results$t, na.rm = T)
      bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
      bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
    }else {
      bootvar <- bootq1 <- bootq2 <- NA
    }

    return(list(est.ps = est.ps, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2, cov.bal = sdm.cov, matching.detail = out1, matching.rate = matching.rate))

  }else if(method == "pg"){
    # if prognostic score model is given and score is missing, then fit the model
    # if prognostic score model and score are both given, then we will use the given score and throw out a warning
    # if prognostic score model is not given but score is invalid, then stop
    # otherwise, directly use score from the arguments
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
    }else if (model.pg == "linpred" && is.null(pg)){
      # if linear predictors are valid, then fit a linear model
      if(is.lp(lp.pg, n)){
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
        stop("invalid linear predictors for propensity score")
      }
    }else if (!(model.pg == "other") && !(is.null(pg))){
      warning("prognostic score is given while a fitting model is chosen")
    }else if (is.pg(pg, n) == F){
      stop("prognostic score should be given or the score is invalid")
    }else{
      # mu0 = pg[,1]
      # mu1 = pg[,2]
      mu0 = pg
    }

    # lm.y <- Y
    # lm.x <- cbind(mu0, mu1, mu0 ^ 2, mu1 ^ 2, mu0 * mu1)
    # lm.out1 <- lm(lm.y[which(A == 1)] ~ lm.x[which(A == 1), ])
    # mu1.seive <- cbind(1, lm.x)[,which(!is.na(lm.out1$coefficients))] %*% lm.out1$coefficients[which(!is.na(lm.out1$coefficients))]
    # sigsqhat1 <- mean((lm.out1$residuals) ^ 2)
    # lm.out0 <- lm(lm.y[which(A == 0)] ~ lm.x[which(A == 0), ])
    # mu0.seive <- cbind(1, lm.x)[,which(!is.na(lm.out0$coefficients))] %*% lm.out0$coefficients[which(!is.na(lm.out0$coefficients))]
    # sigsqhat0 <- mean((lm.out0$residuals) ^ 2)
    #
    # mu1.pg <- mu1.seive
    # mu0.pg <- mu0.seive
    #
    # coeff1.pg <- lm.out1$coefficients
    # coeff0.pg <- lm.out0$coefficients

    Kiw.pg <- rep(0, n)
    thistrt <- 0
    A1 <- A != thistrt
    out1 <- Matching::Match(Y = Y, Tr = A1, X = mu0,
      distance.tolerance = 0, ties = FALSE, Weight = 2,
      caliper = caliper, replace = replace, ...)
    # # bias correction for mu
    # bias.mu0 <- sum(mu0.pg[out1$index.control] - mu0.pg[out1$index.treated]) / n1
    mdata1 <- out1$mdata
    Kiw.pg[loc.0] <- table(factor(out1$index.control, levels = loc.0))
    est.pg <- mean(Y[loc.1]) - mean(mdata1$Y[which(mdata1$Tr == 0)])
    matching.rate = 1 - out1$ndrops / n1

    # covariate balance checking
    # duplicate data generated by matching
    dup.trt <- X[loc.1, ]
    dup.ctl <- matrix(0, 1, ncol(X))

    for (d in loc.0) {
      if(Kiw.pg[d] > 0){
        for (k in 1:Kiw.pg[d]) {
          dup.ctl <- rbind(dup.ctl, X[d, ])
        }
      }
    }

    # remove the first row
    dup.ctl <- dup.ctl[-1,]

    # record standard difference in mean for all covariates
    sdm.cov = rep(0, ncol(X))

    for (c in 1:ncol(X)) {
      res = matrix(0, 3, 2)
      colnames(res) <- c("Before Matching", "After Matching")
      rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
      res[1, 1] <- mean(X[loc.1, c])
      res[2, 1] <- mean(X[loc.0, c])
      res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
      res[1, 2] <- mean(dup.trt[, c])
      res[2, 2] <- mean(dup.ctl[, c])
      res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
      if (cov.balance){
        cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
        print(res)
      }
      sdm.cov[c] = res[3, 2]
    }

    #variance estimation
    data <- cbind(Y, A, Kiw.pg, X)
    if (varest == 1) {
      if (mc == 1) {
        results <- boot::boot(data = data, statistic = estforboot_pgATT, R = boots,
          model.pg = model.pg, pg = pg, lp.pg = lp.pg,
          parallel = "multicore", ncpus = ncpus)
      }
      if (mc == 0) {
        results <- boot::boot(data = data, statistic = estforboot_pgATT, R = boots,
          model.pg = model.pg, pg = pg, lp.pg = lp.pg)
      }

      bootvar = var(results$t, na.rm = T)
      bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
      bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
    }else {
      bootvar <- bootq1 <- bootq2 <- NA
    }

    return(list(est.pg = est.pg, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2, cov.bal = sdm.cov, matching.detail = out1, matching.rate = matching.rate))

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

    Kiw.x <- rep(0, n)
    thistrt <- 0
    A1 <- A != thistrt
    out1 <- Matching::Match(Y = Y, Tr = A1, X = X,
      distance.tolerance = 0, ties = FALSE, Weight = 2,
      caliper = caliper, replace = replace, ...)
    # bias correction for mu
    bias.mu0 <- sum(mu0.x[out1$index.control] - mu0.x[out1$index.treated]) / n1
    mdata1 <- out1$mdata
    Kiw.x[loc.0] <- table(factor(out1$index.control, levels = loc.0))
    est.x <- mean(Y[loc.1]) - mean(mdata1$Y[which(mdata1$Tr == 0)]) + bias.mu0
    matching.rate = 1 - out1$ndrops / n1

    # covariate balance checking
    # duplicate data generated by matching
    dup.trt <- X[loc.1, ]
    dup.ctl <- matrix(0, 1, ncol(X))

    for (d in loc.0) {
      if(Kiw.x[d] > 0){
        for (k in 1:Kiw.x[d]) {
          dup.ctl <- rbind(dup.ctl, X[d, ])
        }
      }
    }

    # remove the first row
    dup.ctl <- dup.ctl[-1,]

    # record standard difference in mean for all covariates
    sdm.cov = rep(0, ncol(X))

    for (c in 1:ncol(X)) {
      res = matrix(0, 3, 2)
      colnames(res) <- c("Before Matching", "After Matching")
      rownames(res) <- c("Treatment group mean....", "Control group mean......", "Standard diff. in mean..")
      res[1, 1] <- mean(X[loc.1, c])
      res[2, 1] <- mean(X[loc.0, c])
      res[3, 1] <- (res[1, 1] - res[2, 1]) / sd(X[, c])
      res[1, 2] <- mean(dup.trt[, c])
      res[2, 2] <- mean(dup.ctl[, c])
      res[3, 2] <- (res[1, 2] - res[2, 2]) / sd(X[, c])
      if (cov.balance){
        cat("\n***** (V", c, ") ", colnames(X)[c], " *****\n", sep = "")
        print(res)
      }
      sdm.cov[c] = res[3, 2]
    }

    #variance estimation
    data <- cbind(Y, A, Kiw.x, X)
    if (varest == 1) {
      if (mc == 1) {
        results <- boot::boot(data = data, statistic = estforboot_covATT, R = boots,
          parallel = "multicore", ncpus = ncpus)
      }
      if (mc == 0) {
        results <- boot::boot(data = data, statistic = estforboot_covATT, R = boots)
      }

      bootvar = var(results$t, na.rm = T)
      bootq1 = quantile(results$t, probs = 0.025, type = 5, na.rm = TRUE)
      bootq2 = quantile(results$t, probs = 0.975, type = 5, na.rm = TRUE)
    }else {
      bootvar <- bootq1 <- bootq2 <- NA
    }

    return(list(est.x = est.x, bootvar = bootvar, bootq1 = bootq1, bootq2 = bootq2, cov.bal = sdm.cov, matching.detail = out1, matching.rate = matching.rate))

  }else{
    stop("invalid method")
  }

}











