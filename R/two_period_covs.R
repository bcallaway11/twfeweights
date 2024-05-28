#-------------------------------------------------------------------------------
#
# File contains functions to compute implicit weights for the setting with
# two time periods.  The two cases covered are:
#
# * implicit TWFE regression weigthts
# * implicit AIPW weights
#
#-------------------------------------------------------------------------------

#' two_period_covs_obj
#' 
#' an object that holds the results of the `two_period_covs` function
#' 
#' @format A list with the following components:
#' \item{est}{the two period overall estimate, e.g., the alpha coefficient from the two period TWFE regression}
#' \item{weights}{the implicit weights}
#' \item{cov_balance_df}{a data frame with the balance statistics for the covariates}
#' 
#' @export
two_period_covs_obj <- function(est, weights, cov_balance_df) {
  out <- list(
    est=est,
    weights=weights,
    cov_balance_df=cov_balance_df
  )
  class(out) <- "two_period_covs_obj"
  out
}


#' @title Two Period TWFE with Covariates as Re-weighting
#'
#' @description A function that computes implicit weights from a TWFE 
#' regression with two periods and covariates.
#' 
#' @param yname the name of the outcome variable
#' @param tname the name of the time variable
#' @param idname the name of the id variable
#' @param gname the name of the group variable.  Units in the untreated group 
#'  should have this coded to be 0.  For treated units, this can be any 
#'  positive number though all treated units should have the same value of 
#'  group.
#' @param xformula a formula for the covariates to include in the model
#' @param data the data frame
#' @param weightsname optional name of the sampling weights variable
#' @param ... additional arguments, currently unused but allows same arguments 
#'  as would be used in `did::att_gt`
#'  
#' @export
two_period_reg_weights <- function(yname,
                                    tname,
                                    idname,
                                    gname,
                                    xformula=~1,
                                    time_invariant_x_formula=NULL,
                                    additional_covariates_formula=NULL,
                                    additional_dcovariates_formula=NULL,
                                    data,
                                    weightsname=NULL,
                                    ...) {
  
  if (length(unique(data[,tname])) != 2) {
    stop("Only two periods are allowed.")
  }
  
  if (length(unique(data[,gname])) != 2) {
    stop("Only two groups are allowed.")
  }
  
  # check for balanced panel
  if (nrow(BMisc::makeBalancedPanel(data, idname, tname)) != nrow(data)) {
    stop("function only works for balanced panel data")
  }
  
  # sort data by id and then by time period
  data <- data[order(data[,idname], data[,tname]),]

  .n <- length(unique(data[,idname]))
  .tlist <- sort(unique(data[,tname]))
  maxT <- max(.tlist)
  minT <- min(.tlist)
    
  # check if sampling weights are time constant for particular units
  if (!is.null(weightsname)) {
    dweights <- data[ data[,tname]==maxT, weightsname ] - data[ data[,tname]==minT, weightsname ]
    if ( !all(dweights==0) ) {
      stop("sampling weights are not time constant for all units")
    }
  }
  
  # ----------------------------------------------------------------------------
  # run the regression and calculate the weights
  # ----------------------------------------------------------------------------
  # add the variables that change over time
  dvars <- lapply(c(yname, BMisc::rhs.vars(xformula)), function(varname) {
    data[ data[,tname]==maxT, varname] - data[ data[,tname]==minT, varname]  
  })
  .df_wide <- do.call(cbind.data.frame, dvars)
  colnames(.df_wide) <- paste0(".d",c(yname, BMisc::rhs.vars(xformula)))
  # this will be the formula for the covariates in the regression later
  reg_xformula <- BMisc::toformula("",paste0("", colnames(.df_wide))[-1])
  
  # add the time invariant variables
  if( !is.null(time_invariant_x_formula) ) {
    # try to account for factors here
    .df_wide_temp <- model.matrix(time_invariant_x_formula, data[ data[,tname]==minT, ])[,-1]
    # recode formula using new names of variables
    reg_xformula <- BMisc::addCovToFormla(colnames(.df_wide_temp), reg_xformula)
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
  }
  
  # add the treatment indicator
  .df_wide$.D <- ifelse(data[ data[,tname]==maxT, gname ]==0, 0, 1)
  
  .D <- .df_wide$.D
  .dy <- .df_wide[, paste0(".d", yname)]
  if (!is.null(weightsname)) {
    .sampling_weights <- data[ data[,tname]==maxT, weightsname ] / mean(data[ data[,tname]==maxT, weightsname ])
  } else {
    .sampling_weights <- rep(1, .n)
  }
  .df_wide$.sampling_weights <- .sampling_weights
  p <- weighted.mean(.D, w=.sampling_weights)
  dreg <- lm(BMisc::toformula(".D", BMisc::rhs.vars(reg_xformula)), data=.df_wide, weights=.sampling_weights)
  dresid <- .D - predict(dreg)
  twfe_reg <- lm(BMisc::toformula(paste0(".d",yname), c(".D",BMisc::rhs.vars(reg_xformula))), data=.df_wide, weights=.sampling_weights)
  twfe_alp <- coef(twfe_reg)[2]
  pm <- .D - (1-.D)
  .regression_weights <- dresid/weighted.mean(dresid^2, w=.sampling_weights) # from FWL
  
  alp_using_weights <- weighted.mean(.regression_weights*.dy, w=.sampling_weights) # should be equal to alpha from TWFE
  if (abs(twfe_alp - alp_using_weights) > 1e-6) {
    warning("Looks like something has gone wrong...")
  }
  
  # calculate balance for all supplied covariates
  # data.table::setDT(.df_wide)
  # .df_wide <- data.table::dcast(data, get(idname) ~ get(tname), value.var=BMisc::rhs.vars(xformula))
  .df_wide_temp <- data[ data[,tname]==minT, BMisc::rhs.vars(xformula), drop=FALSE]
  colnames(.df_wide_temp) <- paste0(BMisc::rhs.vars(xformula), "_", minT)
  # if requested, additional covariates from first period
  if ( !is.null(additional_covariates_formula) ) {
    # try to account for factors here
    .df_wide_temp2 <- model.matrix(additional_covariates_formula, data[ data[,tname]==minT, ])[,-1,drop=FALSE]
    #.df_wide_temp2 <- data[ data[,tname]==minT, BMisc::rhs.vars(additional_covariates_formula), drop=FALSE]
    colnames(.df_wide_temp2) <- paste0(colnames(.df_wide_temp2), "_", minT)
    .df_wide_temp <- cbind.data.frame(.df_wide_temp, .df_wide_temp2)
  }
  # if requested, the change in additional covariates over time
  if ( !is.null(additional_dcovariates_formula) ) {
    .df_wide_temp3 <- data[ data[,tname]==maxT, BMisc::rhs.vars(additional_dcovariates_formula), drop=FALSE] - 
      data[ data[,tname]==minT, BMisc::rhs.vars(additional_dcovariates_formula), drop=FALSE]
    colnames(.df_wide_temp3) <- paste0("d",BMisc::rhs.vars(additional_dcovariates_formula))
    .df_wide_temp <- cbind.data.frame(.df_wide_temp, .df_wide_temp3)
  }
  # merge the additional covariates with the main data and drop the 
  # treatment indicator and the sampling weights
  .df_wide <- cbind.data.frame(.df_wide[,1:(ncol(.df_wide)-2)], .df_wide_temp)
  
  this_covs <- colnames(.df_wide)
  cov_balance <- lapply(this_covs, function(covariate) {
    unweighted_treat <- weighted.mean(.df_wide[,covariate][.D==1], w=.sampling_weights[.D==1])
    unweighted_untreat <- weighted.mean(.df_wide[,covariate][.D==0], w=.sampling_weights[.D==0])
    unweighted_diff <- unweighted_treat - unweighted_untreat
    
    weighted_treat <- weighted.mean(pm[.D==1]*.regression_weights[.D==1]*.df_wide[,covariate][.D==1], w=.sampling_weights[.D==1])*p
    weighted_untreat <- weighted.mean(pm[.D==0]*.regression_weights[.D==0]*.df_wide[,covariate][.D==0], w=.sampling_weights[.D==0])*(1-p)
    weighted_diff <- weighted_treat - weighted_untreat
    
    sd <- sqrt(weighted.mean((.df_wide[,covariate] - weighted.mean(.df_wide[,covariate], w=.sampling_weights))^2, w=.sampling_weights))
    
    list(covariate=covariate, 
         unweighted_treat=unweighted_treat,
         unweighted_untreat=unweighted_untreat,
         unweighted_diff=unweighted_diff, 
         weighted_treat=weighted_treat,
         weighted_untreat=weighted_untreat,
         weighted_diff=weighted_diff,
         sd=sd)
  })
  cov_balance_df <- do.call(rbind.data.frame, cov_balance)
  
  
  # re-normalize the weights to be in line with other functions
  .regression_weights[.D==1] <- .regression_weights[.D==1]*.sampling_weights[D==1]
  .regression_weights[.D==1] <- .regression_weights[.D==1]/mean(.regression_weights[.D==1])
  .regression_weights[.D==0] <- .regression_weights[.D==0]*.sampling_weights[D==0]
  .regression_weights[.D==0] <- .regression_weights[.D==0]/mean(.regression_weights[.D==0])
  out <- two_period_covs_obj(twfe_alp, .regression_weights, cov_balance_df)
  out
}

#' a function to plot the balance statistics from the `two_period_covs` function
#' 
#' @param two_period_covs_obj the object returned from the `two_period_covs` function
#' @param ... additional arguments to pass to `ggplot2::geom_point`
#' 
#' @export
ggtwfeweights.two_period_covs_obj <- function(two_period_covs_obj) {
  cov_balance_df <- two_period_covs_obj$cov_balance_df
  cov_balance_df$covariate <- factor(cov_balance_df$covariate, levels = rev(unique(cov_balance_df$covariate)))
  ggplot(cov_balance_df, aes(y=covariate, x=abs(weighted_diff/sd))) + 
    geom_point(color="steelblue", size=3, shape=16) +
    geom_point(aes(x=abs(unweighted_diff/sd)), color="red", size=3, shape=1) +
    theme_bw()
}


#' @title Two Period AIPW Weights
#' 
#' @description Computes implicit weights coming from AIPW estimation 
#'  of ATT with two periods.
#'  
#' @inheritParams two_period_reg_weights
#' @export
two_period_aipw_weights <- function(yname,
                                    tname,
                                    idname,
                                    gname,
                                    xformula=~1,
                                    d_covs_formula=NULL,
                                    pscore_formula=xformula,
                                    pscore_d_covs_formula=d_covs_formula,
                                    extra_balance_vars_formula=NULL,
                                    extra_balance_d_vars_formula=NULL,
                                    data,
                                    weightsname=NULL,
                                    ...) {
  

  if (length(unique(data[,tname])) != 2) {
    stop("Only two periods are allowed.")
  }
  
  if (length(unique(data[,gname])) != 2) {
    stop("Only two groups are allowed.")
  }
  
  # check for balanced panel
  if (nrow(BMisc::makeBalancedPanel(data, idname, tname)) != nrow(data)) {
    stop("function only works for balanced panel data")
  }
  
  # sort data by id and then by time period
  data <- data[order(data[,idname], data[,tname]),]
  
  n <- length(unique(data[,idname]))
  tlist <- sort(unique(data[,tname]))
  maxT <- max(tlist)
  minT <- min(tlist)
  
  # check if sampling weights are time constant for particular units
  if (!is.null(weightsname)) {
    dweights <- data[ data[,tname]==maxT, weightsname ] - data[ data[,tname]==minT, weightsname ]
    if ( !all(dweights==0) ) {
      stop("sampling weights are not time constant for all units")
    }
  }
  
  # ----------------------------------------------------------------------------
  # run the regression and calculate the weights
  # ----------------------------------------------------------------------------
  # add outcome
  .df_wide <- data[ data[,tname]==maxT, yname, drop=FALSE] - data[ data[,tname]==minT, yname, drop=FALSE]
  colnames(.df_wide) <- paste0(".d",yname)
  # add treatment
  .df_wide$.D <- ifelse(data[ data[,tname]==maxT, gname ]==0, 0, 1)
  # add the time invariant variables
  if (!is.null(xformula)) {
    # try to account for factors here
    .df_wide_temp <- model.matrix(xformula, data[ data[,tname]==minT, ])[,-1, drop=FALSE]
    colnames(.df_wide_temp) <- paste0(colnames(.df_wide_temp), "_", minT)
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    reg_xformula <- BMisc::toformula("", colnames(.df_wide_temp))
  }
  # add the variables that change over time
  if (!is.null(d_covs_formula)) {
    dvars <- lapply(c(BMisc::rhs.vars(d_covs_formula)), function(varname) {
      data[ data[,tname]==maxT, varname] - data[ data[,tname]==minT, varname]  
    })
    .df_wide_temp <- do.call(cbind.data.frame, dvars)
    colnames(.df_wide_temp) <- paste0(".d",c(BMisc::rhs.vars(d_covs_formula)))
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    # this will be the formula for the covariates in the regression later
    reg_xformula <- BMisc::addCovToFormla(colnames(.df_wide_temp), reg_xformula)
  }
  # add any additional/alternative variables for the propensity score
  if (pscore_formula != xformula) {
    .df_wide_temp <- model.matrix(pscore_formula, data[ data[,tname]==minT, ])[,-1, drop=FALSE]
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    pscore_xformula <- BMisc::toformula("", colnames(.df_wide_temp))
  } 
  if ( isTRUE(pscore_d_covs_formula != d_covs_formula) ) {
    dvars <- lapply(c(BMisc::rhs.vars(pscore_d_covs_formula)), function(varname) {
      data[ data[,tname]==maxT, varname] - data[ data[,tname]==minT, varname]  
    })
    .df_wide_temp <- do.call(cbind.data.frame, dvars)
    colnames(.df_wide_temp) <- paste0(".d",c(BMisc::rhs.vars(pscore_d_covs_formula)))
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    pscore_xformula <- BMisc::addCovToFormla(colnames(.df_wide_temp), pscore_xformula)
  }
  if (pscore_formula == xformula & !isFALSE(pscore_d_covs_formula == d_covs_formula) ) {
    pscore_xformula <- reg_xformula
  }
  
  D <- .df_wide$.D
  dy <- .df_wide[, paste0(".d", yname)]
  if (!is.null(weightsname)) {
    sampling_weights <- data[ data[,tname]==maxT, weightsname ] / mean(data[ data[,tname]==maxT, weightsname ])
  } else {
    sampling_weights <- rep(1, .n)
  }
  .df_wide$.sampling_weights <- sampling_weights
  p <- weighted.mean(D, w=sampling_weights)
  pm <- D - (1-D)
  full_pscore_formula <- BMisc::toformula("D", BMisc::rhs.vars(pscore_xformula))
  environment(full_pscore_formula) <- environment() # not sure why I need this
  pscore <- predict(glm(full_pscore_formula, 
                        data=.df_wide, 
                        family=binomial, 
                        weights=sampling_weights), 
                    type="response")
    
  ipw_weights0 <- (pscore/(1-pscore)*(1-p)/p)[D==0]
  ipw_weights1 <- rep(1, sum(D==1))
  # => normalize jointly instead of one by one
  # weighted.mean(dy[D==1], w=ipw_weights1*sampling_weights[D==1]) -
  #  weighted.mean(dy[D==0], w=ipw_weights0*sampling_weights[D==0])# this is ipw and confirmed same as drdid
  .df_wide$.pscore <- pscore
  .df_wide$.sampling_weights <- sampling_weights
  .df_wide$.odds_ratio <- pscore/(1-pscore)
  gamma0_tilde <- coef(lm(BMisc::toformula(".odds_ratio", BMisc::rhs.vars(reg_xformula)), 
                       data=.df_wide[D==0,],
                       weights=.sampling_weights))
  X <- model.matrix(BMisc::toformula("", BMisc::rhs.vars(reg_xformula)), data=.df_wide)
  X0 <- X[D==0,]
  X1 <- X[D==1,]
  Qxx <- t(sampling_weights[D==0]*X0)%*%X0/sum(D==0)
  gamma0 <- solve(Qxx)%*%colMeans(sampling_weights[D==1]*X1)/mean(sampling_weights[D==1])*p/(1-p)
  aipw_extra_weights0 <- (1-p)/p * X0 %*% (gamma0 - gamma0_tilde)
  aipw_weights0 <- ipw_weights0 + aipw_extra_weights0
  aipw_weights_att <- weighted.mean(dy[D==1], w=ipw_weights1*sampling_weights[D==1]) -
    weighted.mean(dy[D==0], w=aipw_weights0*sampling_weights[D==0])
  
  # check components
  # - avg change for treated group
  t1 <- weighted.mean(dy[D==1], w=ipw_weights1*sampling_weights[D==1]); t1 # yes!
  # - regression adjustment term
  t2 <- weighted.mean(dy[D==0], w=as.numeric(X0%*%gamma0)*sampling_weights[D==0]); t2# no!
  # - ipw term
  t3 <- weighted.mean(dy[D==0], w=ipw_weights0*sampling_weights[D==0]); t3 # yes!
  # - cross term, ipw and regression adjustment
  t4 <- weighted.mean(dy[D==0], w=as.numeric(X0%*%gamma0_tilde)*sampling_weights[D==0]); t4 # no!
  aipw_weights_att <- t1 - t2 - t3 + t4; t1-t2-t3+t4 # this should be equal to the final result
  ###
  
  aipw_weights_treat <- ipw_weights1*sampling_weights[D==1]
  aipw_weights_treat <- aipw_weights_treat/mean(aipw_weights_treat)
  aipw_weights_untreat1 <- as.numeric(X0%*%gamma0)*sampling_weights[D==0]
  aipw_weights_untreat1 <- aipw_weights_untreat1/mean(aipw_weights_untreat1)
  aipw_weights_untreat2 <- ipw_weights0*sampling_weights[D==0]
  aipw_weights_untreat2 <- aipw_weights_untreat2/mean(aipw_weights_untreat2)
  aipw_weights_untreat3 <- as.numeric(X0%*%gamma0_tilde)*sampling_weights[D==0]
  aipw_weights_untreat3 <- aipw_weights_untreat3/mean(aipw_weights_untreat3)
  aipw_weights <- rep(1, nrow(.df_wide))
  aipw_weights[D==1] <- aipw_weights_treat
  aipw_weights[D==0] <- (aipw_weights_untreat1 + aipw_weights_untreat2 - aipw_weights_untreat3)

  # can calculate ATT as follows
  mean(aipw_weights[D==1]*dy[D==1]) - mean(aipw_weights[D==0]*dy[D==0])
  # compute AIPW using ipw functions so that we can allow for reg and pscore
  # covariates to be different
  out_reg_formula <- BMisc::toformula(paste0(".d",yname), BMisc::rhs.vars(reg_xformula))
  out_reg <- lm(out_reg_formula, 
                data=.df_wide[D==0,], 
                weights=.sampling_weights)
  L0 <- predict(out_reg, newdata=.df_wide)
  
  att <- DRDID::std_ipw_did_panel(dy-L0, rep(0,nrow(.df_wide)), D, 
                                  covariates=model.matrix(BMisc::toformula("", BMisc::rhs.vars(reg_xformula)), data=.df_wide), 
                                  i.weights=sampling_weights)$ATT
  
  # DRDID::drdid_panel(dy, rep(0,nrow(.df_wide)), D, 
  #                   covariates=model.matrix(BMisc::toformula("", BMisc::rhs.vars(reg_xformula)), data=.df_wide), 
  #                   i.weights=sampling_weights)$ATT
  
  if (abs(aipw_weights_att - att) > 1e-6) {
    warning("Looks like something has gone wrong...")
  }
  
  # calculate balance for all supplied covariates
  # data.table::setDT(.df_wide)
  # .df_wide <- data.table::dcast(data, get(idname) ~ get(tname), value.var=BMisc::rhs.vars(xformula))
  #.df_wide_temp <- data[ data[,tname]==minT, BMisc::rhs.vars(xformula), drop=FALSE]
  #colnames(.df_wide_temp) <- paste0(BMisc::rhs.vars(xformula), "_", minT)
  
  # drop no longer needed columns
  .df_wide <- .df_wide[,-c(2,ncol(.df_wide):(ncol(.df_wide)-2))]
  
  # if requested, additional covariates from first period
  if ( !is.null(extra_balance_vars_formula) ) {
    # try to account for factors here
    .df_wide_temp2 <- model.matrix(extra_balance_vars_formula, data[ data[,tname]==minT, ])[,-1,drop=FALSE]
    #.df_wide_temp2 <- data[ data[,tname]==minT, BMisc::rhs.vars(additional_covariates_formula), drop=FALSE]
    colnames(.df_wide_temp2) <- paste0(colnames(.df_wide_temp2), "_", minT)
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp2)
  }
  # if requested, the change in additional covariates over time
  if ( !is.null(extra_balance_d_vars_formula) ) {
    .df_wide_temp3 <- data[ data[,tname]==maxT, BMisc::rhs.vars(extra_balance_d_vars_formula), drop=FALSE] - 
      data[ data[,tname]==minT, BMisc::rhs.vars(extra_balance_d_vars_formula), drop=FALSE]
    colnames(.df_wide_temp3) <- paste0("d",BMisc::rhs.vars(extra_balance_d_vars_formula))
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp3)
  }
  
  this_covs <- colnames(.df_wide)
  cov_balance <- lapply(this_covs, function(covariate) {
    unweighted_treat <- weighted.mean(.df_wide[,covariate][D==1], w=sampling_weights[D==1])
    unweighted_untreat <- weighted.mean(.df_wide[,covariate][D==0], w=sampling_weights[D==0])
    unweighted_diff <- unweighted_treat - unweighted_untreat
    
    weighted_treat <- weighted.mean(.df_wide[,covariate][D==1], w=ipw_weights1*sampling_weights[D==1]); weighted_treat # yes!
    # - regression adjustment term
    weighted_untreat1 <- weighted.mean(.df_wide[,covariate][D==0], w=as.numeric(X0%*%gamma0)*sampling_weights[D==0]); weighted_untreat1# no!
    # - ipw term
    weighted_untreat2 <- weighted.mean(.df_wide[,covariate][D==0], w=ipw_weights0*sampling_weights[D==0]); weighted_untreat2 # yes!
    # - cross term, ipw and regression adjustment
    weighted_untreat3 <- weighted.mean(.df_wide[,covariate][D==0], w=as.numeric(X0%*%gamma0_tilde)*sampling_weights[D==0]); weighted_untreat3 # no!
    weighted_untreat <- weighted_untreat1 + weighted_untreat2 - weighted_untreat3    
    
    weighted_diff <- weighted_treat - weighted_untreat
    
    sd <- sqrt(weighted.mean((.df_wide[,covariate] - weighted.mean(.df_wide[,covariate], w=sampling_weights))^2, w=sampling_weights))
    
    list(covariate=covariate, 
         unweighted_treat=unweighted_treat,
         unweighted_untreat=unweighted_untreat,
         unweighted_diff=unweighted_diff, 
         weighted_treat=weighted_treat,
         weighted_untreat=weighted_untreat,
         weighted_diff=weighted_diff,
         sd=sd)
  })
  cov_balance_df <- do.call(rbind.data.frame, cov_balance)
  
  out <- two_period_covs_obj(att, aipw_weights, cov_balance_df)
  out
}