# Code for computing TWFE weights with multiple
# periods and variation in treatment timing

#' @title gt_weights
#' @description a class for holding and working with group-time weights
#' @inheritParams implicit_twfe_weights_gt
#' @param treated a vector of units that make up the effective treated group for group g 
#'  and period tp
#' @param comparison a vector of units that make up the effective comparison 
#'  group for group g and period tp
#' @param base_period options are "first_period" or "gmin1".  If "first_period", 
#'  the outcome for the treated group is subtracted by the outcome in the first 
#'  period.  If "gmin1", the outcome for the treated group is subtracted by the
#'  outcome in the period right before treatment for group g.  The default 
#'  is "first_period".
#' @param weights_treated a vector of weights for the treated group g in time 
#'  period tp
#' @param weights_comparison a vector of weights for the comparison group for 
#'  group g in time period tp
#' @param weighted_outcome_treated the mean outcome for the treated group 
#'  after applying weights_treated to it
#' @param weighted_outcome_comparison the mean outcome for the comparison group 
#'  after applying weights_comparison to it
#' @param weighted_outcome_diff the difference between the mean outcomes for 
#'  the treated and comparison groups after applying weights to them
#' @param remainder when `base_period` is "gmin1", then the decompostion 
#'  doesn't exactly add up to alpha from the TWFE regression.  This term 
#'  contains the leftover part of the decomposition.  Typically, this is 
#'  very small.  
#' @param alpha_weight the contribution of this term to the estimate of alpha
#' @export
gt_weights <- function(g,
                       tp,
                       treated=NULL,
                       comparison=NULL,
                       weights_treated,
                       weights_comparison,
                       weighted_outcome_treated=NULL,
                       weighted_outcome_comparison=NULL,
                       weighted_outcome_diff=NULL,
                       base_period=NULL,
                       remainder=NULL,
                       cov_bal_df=NULL,
                       alpha_weight) {
  
  ret_gt_weights <- list(g=g,
                         tp=tp,
                         treated=treated,
                         comparison=comparison,
                         weights_treated=weights_treated,
                         weights_comparison=weights_comparison,
                         weighted_outcome_treated=weighted_outcome_treated,
                         weighted_outcome_comparison=weighted_outcome_comparison,
                         weighted_outcome_diff=weighted_outcome_diff,
                         base_period=base_period,
                         remainder=remainder,
                         cov_bal_df=cov_bal_df,
                         alpha_weight=alpha_weight)
  
  class(ret_gt_weights) <- "gt_weights"
  ret_gt_weights
}

#' @title implicit_twfe_weights_gt
#' @description Compute TWFE weights for a particular group
#'  and time period.  This function is developed for the 
#'  case where there is a binary treatment and staggered 
#'  treatment adoption.  The weights are for the coefficient 
#'  on the binary treatment variable and allows for 
#'  the regression to include additional time varying 
#'  covariates.  
#' 
#' This is mainly an internal function and data should already be sorted, 
#' balanced, and time periods and groups should be arranged like 1,2,3,etc.
#'
#' @param g the current group for which the weights will be computed
#' @param tp the time period for which the weights will be computed
#' @inheritParams did::att_gt
#' @param xformula formula for covariates that are included in the 
#'  TWFE regression.
#' @param xname by default, `implicit_twfe_weights_gt` applies the 
#'  implicit TWFE regression weights to the outcome.  In 
#'  cases where the researcher would like to apply the 
#'  implicit regression weights to some other variable, the 
#'  name of this variable can be passed in through this 
#'  argument. 
#'  
#'  This is typically used as a way to check how well 
#'  the TWFE regression balances levels of time-varying 
#'  covariates, time-invariant covariates, or other variables
#'  not included in the model.  This is often useful 
#'  as a diagnostic to assess how sensitive the TWFE regression 
#'  is to violations of linearity.
#'  
#'  @return a gt_weights object
#'  @export
implicit_twfe_weights_gt <- function(g,
                                     tp,
                                     yname,
                                     tname,
                                     idname,
                                     gname,
                                     base_period,
                                     xformula=NULL,
                                     data,
                                     weightsname) {
  
  TP <- data[,tname] # a vector of the same length as the full panel data, nT
  G <- data[,gname]  # a vector of the same length as the full panel data, nT
  thisG <- G[TP==tp] # a group vector of length n (# cross-sectional units)
  minT <- min(unique(TP))
  maxT <- max(unique(TP))
  nT <- length(unique(TP))
  if (!is.null(weightsname)) {
    sampling_weights <- data[,weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }
  
  if (g==0 | g==(maxT+1)) {
    stop("not computng this for untreated group")
  }
  
  Yit <- data[,yname]
  xformula <- BMisc::addCovToFormla("-1", xformula)
  Xit <- model.matrix(xformula, data=data)
  Dit <- 1*( (data[,tname] >= data[,gname]) & (data[,gname] != 0) )
  
  # double demean variables in the model
  ddotDit <- fixest::demean(Dit, data[,c(idname,tname)], weights=sampling_weights) %>% as.numeric
  ddotYit <- fixest::demean(Yit, data[,c(idname,tname)], weights=sampling_weights) %>% as.numeric
  ddotXit <- fixest::demean(Xit, data[,c(idname,tname)], weights=sampling_weights)
  
  # weights depend on linear projection of treatment
  # on covariates using all periods
  lp_treat <- lm(ddotDit ~ -1 + ddotXit, weights=sampling_weights)
  gam <- as.matrix(coef(lp_treat))
  
  idx <- G==g & TP==tp
  Yit <- get_Yit(data, tp, idname, yname, tname)[TP==tp]
  outcomeit <- Yit
  
  # this gets Y_it - Y_ig-1 for the g being considered here (if we want it)
  # note to self: TP==tp does nothing except get the length right
  Yigmin1 <- get_Yit(data, g-1, idname, yname, tname)[TP==tp] 
  
  # if you want to subtract off the Yi1, it will give you the same thing
  Yi1 <- get_Yi1(data, idname, yname, tname, gname)[TP==1]
  
  if (base_period == "first_period") {
    outcomeit <- Yit - Yi1
  } else if (base_period == "gmin1") {
    outcomeit <- Yit - Yigmin1
  } 
  
  gpart_w_inner <- (ddotDit - mv_mult(ddotXit, gam))[idx]
  gpart_w <- gpart_w_inner / weighted.mean(gpart_w_inner, w=sampling_weights[idx])
  gpart <- weighted.mean( gpart_w * outcomeit[thisG==g], w=sampling_weights[idx] )
  
  idxU <- G==0 & TP==tp
  upart_w_inner <- (ddotDit - mv_mult(ddotXit, gam))[idxU]
  upart_w <- upart_w_inner / weighted.mean(upart_w_inner, w=sampling_weights[idxU])
  upart <- weighted.mean( upart_w * outcomeit[thisG==0] , w=sampling_weights[idxU])
  
  reg_attgt <- gpart - upart
  
  alpha_weight <- combine_twfe_weights_gt(g,tp,gname,tname,xformula,idname,data,weightsname)
  
  remainder <- 0
  if (base_period == "gmin1") {
    # works but try to simplify
    # remainder <- mean( gpart_w * (Yigmin1-Yi1)[thisG==g] ) - mean( upart_w * (Yigmin1-Yi1)[thisG==0] )
    # also works but try to simplify more
    # remainder <- - mean( upart_w * (Yigmin1-Yi1)[thisG==0] )
    remainder <- - weighted.mean( upart_w * Yigmin1[thisG==0], w=sampling_weights[idxU])
  }
  
  out <- gt_weights(g=g,
                    tp=tp,
                    treated=outcomeit[thisG==g],
                    comparison=outcomeit[thisG==0],
                    base_period=base_period,
                    weights_treated=gpart_w,
                    weights_comparison=upart_w,
                    weighted_outcome_treated=gpart,
                    weighted_outcome_comparison=upart,
                    weighted_outcome_diff=reg_attgt,
                    remainder=remainder,
                    alpha_weight=alpha_weight)
  out
}

#' @title combine_twfe_weights_gt
#' @description a function that recovers the weights on each implicit 
#'  group-time average treatment effect from the TWFE regression
#' @inheritParams implicit_twfe_weights_gt
#' @export
combine_twfe_weights_gt <- function(g,
                                    tp,
                                    gname,
                                    tname,
                                    xformula,
                                    idname,
                                    data,
                                    weightsname=NULL) {
  
  G <- data[,gname]
  TP <- data[,tname]
  maxT <- max(TP)
  if (!is.null(weightsname)) {
    sampling_weights <- data[,weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }
  
  xformula <- BMisc::addCovToFormla("-1", xformula)
  Xit <- model.matrix(xformula, data=data)
  Dit <- 1*( (data[,tname] >= data[,gname]) & (data[,gname] != 0) )
  
  # double demean variables in the model
  ddotDit <- fixest::demean(Dit, data[,c(idname,tname)], weights=sampling_weights) %>% as.numeric()
  ddotXit <- fixest::demean(Xit, data[,c(idname,tname)], weights=sampling_weights)
  
  # weights depend on linear projection of treatment
  # on covariates using all periods
  lp_treat <- lm(ddotDit ~ -1 + ddotXit, weights=sampling_weights)
  gam <- as.matrix(coef(lp_treat))
  
  if (g==(maxT+1)) {
    g <- 0
  }
  
  idx <- G==g & TP==tp
  pg <- weighted.mean(G==g, w=sampling_weights)
  
  alp_den <- weighted.mean((ddotDit - mv_mult(ddotXit,gam))*(ddotDit), w=sampling_weights)
  weighted.mean( ( ddotDit - mv_mult(ddotXit, gam) )[idx] , w=sampling_weights[idx]) * pg / (alp_den*maxT)
  
  # 
  # 
  # G <- data[,gname]
  # TP <- data[,tname]
  # tlist <- sort(unique(TP))
  # local_g <- which(tlist == g) # the group in terms of t=1,2,3,...,T
  # minT <- min(unique(TP))
  # nT <- length(unique(TP))
  # pbarg <- mean( G[G!=0] == g)
  # ntreatedperiods <- nT - local_g + 1
  # w_gt <- pbarg/ntreatedperiods
  # if (tp < g) w_gt <- -w_gt
  # w_gt
}

#' @title implicit_twfe_weights 
#' @description a function to compute implicit TWFE weights when the TWFE 
#'  model includes covariates. 
#' 
#' @inheritParams implicit_twfe_weights_gt
#' @return a `decomposed_twfe` object
#' @export
implicit_twfe_weights <- function(yname,
                                  tname,
                                  idname,
                                  gname,
                                  base_period="first_period",
                                  xformula=NULL,
                                  extra_balance_vars_formula=NULL,
                                  data,
                                  weightsname=NULL) {
  
  # confirm balanced panel and sort data
  norig <- nrow(data) 
  data <- BMisc::makeBalancedPanel(data, idname=idname, tname=tname)
  data <- dplyr::arrange(data, idname, tname) %>% as.data.frame()
  if (nrow(data) != norig) {
    warning(paste0("Dropped ", norig-nrow(data), " rows of data to balance the data."))
  }
  
  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[,tname]))
  glist <- sort(unique(data[,gname]))
  data$.TP <- orig2t(data[,tname], tlist)
  data$.G <- orig2t(data[,gname], tlist)
  
  tlist <- sort(unique(data$.TP))
  glist <- sort(unique(data$.G))
  # glist <- glist[glist!=0] # drop never-treated group
  gt_mat <- as.data.frame(expand.grid(glist[-1],tlist))
  gt_mat <- gt_mat %>% arrange_all() %>% as.data.frame()
  colnames(gt_mat) <- c("group","time.period")
  
  twfe_gt <- lapply(1:nrow(gt_mat), function(i) {
    implicit_twfe_weights_gt(g=gt_mat[i,1],
                             tp=gt_mat[i,2],
                             yname=yname,
                             tname=".TP",
                             idname=idname,
                             gname=".G",
                             base_period=base_period,
                             xformula=xformula,
                             data=data,
                             weightsname=weightsname)
  })  
  
  out <- list()
  out$twfe_gt <- twfe_gt
  out$idname <- idname
  out$tname <- tname
  out$gname <- gname
  out$weightsname <- weightsname
  out$data <- data
  class(out) <- "decomposed_twfe"
  out
}

#'@title summary.decomposed_twfe
#'@description summarizes a two-way fixed effects regression 
#' decomposition
#'@param object a decomposed_twfe object
#'@param ... extra arguments, not used here
#'
#'@export
summary.decomposed_twfe <- function(object, ...) {
  object <- object$twfe_gt
  twfe_att_gt <- unlist(BMisc::getListElement(object, "weighted_outcome_diff"))
  group <- unlist(BMisc::getListElement(object, "g"))
  time.period <- unlist(BMisc::getListElement(object, "tp"))
  alpha_weight <- unlist(BMisc::getListElement(object, "alpha_weight"))
  post <- 1*(time.period >= group)
  print_df <- cbind.data.frame(group, time.period, post, twfe_att_gt, alpha_weight)
  alpha <- sum(alpha_weight*twfe_att_gt)
  pt_violations_bias <- sum(alpha_weight[post==0]*twfe_att_gt[post==0])
  remainders <- unlist(BMisc::getListElement(object, "remainder"))
  rem <- sum(remainders*alpha_weight)
  cat("alpha twfe:          ", round(alpha+rem,4), "\n")
  cat("alpha decomp:        ", round(alpha,4), "\n")
  cat("alpha post only:     ", round(sum(alpha_weight[post==1]*twfe_att_gt[post==1]),4), "\n")
  cat("pta violations bias: ", round(pt_violations_bias,4),"\n")
  cat("remainders:          " , round(rem,4), "\n")
  cat("\n")
  # this is basically just a check to see if the balance statistics have been computed
  # if not, we'll print the summary of the decomposition with respect 
  # to the actual outcome.  Otherwiser, we will summarize the covariate balance
  if (is.null(object[[1]]$cov_bal_df)) {
    cat("summary of twfe decomposition\n")
    print_df
  } else {
    cat("summary of twfe covariate balance:\n")
    cov_balance_summary_df <- mp_covariate_bal_summary_helper(object)
    cov_balance_summary_df
  }
}

#' @export
twfe_cov_bal <- function(decomposed_twfe, extra_balance_vars_formula) {
  idname <- decomposed_twfe$idname
  tname <- decomposed_twfe$tname
  gname <- decomposed_twfe$gname
  data <- decomposed_twfe$data
  weightsname <- decomposed_twfe$weightsname
  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[,tname]))
  glist <- sort(unique(data[,gname]))
  data$.TP <- orig2t(data[,tname], tlist)
  data$.G <- orig2t(data[,gname], tlist)
  tlist <- sort(unique(data$.TP))
  glist <- sort(unique(data$.G))
  
  # construct a data.frame that has the mean of each extra balance variable
  extra_vars_frame <- as.data.frame(model.matrix(BMisc::addCovToFormla("-1", extra_balance_vars_formula), data=data))
  mean_extra_var_list <- lapply(1:ncol(extra_vars_frame), function(j) {
    this_colname <- colnames(extra_vars_frame)[j]
    this_frame <- cbind.data.frame(extra_vars_frame[,j], data[,c(idname,tname,gname)])
    colnames(this_frame)[1] <- this_colname
    BMisc::get_Yibar(this_frame, idname=idname, yname=this_colname)
  })  
  mean_extra_vars <- do.call(cbind.data.frame, mean_extra_var_list)
  colnames(mean_extra_vars) <- paste0("mean_", colnames(extra_vars_frame))
  
  # compute balance statistics for each extra balance variable
  groups <- unlist(BMisc::getListElement(decomposed_twfe$twfe_gt, "g"))
  time_periods <- unlist(BMisc::getListElement(decomposed_twfe$twfe_gt, "tp"))
  gt_mat <- cbind(groups, time_periods)
  
  for (i in 1:nrow(gt_mat)) {
    this_balance_df <- twfe_cov_bal_gt(g = gt_mat[i, 1],
                                       tp = gt_mat[i, 2],
                                       decomposed_twfe = decomposed_twfe,
                                       X = mean_extra_vars)
    decomposed_twfe$twfe_gt[[i]]$cov_bal_df <- this_balance_df
  }
  
  # return the decomposed_twfe object with the balance statistics added
  decomposed_twfe
}

#' @export
twfe_cov_bal_gt <- function(g, tp, decomposed_twfe, X) {
  idname <- decomposed_twfe$idname
  tname <- decomposed_twfe$tname
  gname <- decomposed_twfe$gname
  data <- decomposed_twfe$data
  weightsname <- decomposed_twfe$weightsname
  groups <- unlist(BMisc::getListElement(decomposed_twfe$twfe_gt, "g"))
  time_periods <- unlist(BMisc::getListElement(decomposed_twfe$twfe_gt, "tp"))
  this_gt_idx <- which(groups==g & time_periods==tp)
  this_decomposed_twfe <- decomposed_twfe$twfe_gt[[this_gt_idx]]
  weights_treated <- this_decomposed_twfe$weights_treated
  weights_comparison <- this_decomposed_twfe$weights_comparison
  treated_idx <- data$.G==g & data$.TP==tp
  comparison_idx <- data$.G==0 & data$.TP==tp
  if(!is.null(weightsname)) {
    sampling_weights <- data[,weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }
  
  unweighted_covs_treated <- apply(X[treated_idx, , drop=FALSE], 2, function(x) {
    weighted.mean(x, w=sampling_weights[treated_idx])
  })
  unweighted_covs_comparison <- apply(X[comparison_idx, , drop=FALSE], 2, function(x) {
    weighted.mean(x, w=sampling_weights[comparison_idx])
  })
  unweighted_diff <- unweighted_covs_treated - unweighted_covs_comparison
  
  weighted_covs_treated <- apply(X[treated_idx, , drop=FALSE], 2, function(x) {
    weighted.mean(x*weights_treated, w=sampling_weights[treated_idx])
  })
  weighted_covs_comparison <- apply(X[comparison_idx, , drop=FALSE], 2, function(x) {
    weighted.mean(x*weights_comparison, w=sampling_weights[comparison_idx])
  })
  weighted_diff <- weighted_covs_treated - weighted_covs_comparison
  
  both_idx <- treated_idx | comparison_idx
  sd <- apply(X, 2, function(x) {
    sqrt(weighted.mean((x[both_idx] - weighted.mean(x[both_idx], w=sampling_weights[both_idx]))^2, w=sampling_weights[both_idx]))
  })  
  
  # return data frame that contains the weighted covariates for the treated and comparison groups
  cbind.data.frame(unweighted_covs_treated,
                   unweighted_covs_comparison,
                   unweighted_diff,
                   weighted_covs_treated, 
                   weighted_covs_comparison,
                   weighted_diff,
                   sd)
}

#' @export
mp_covariate_bal_summary_helper <- function(decomposed_twfe_obj) {
  
  # code to get this to work for both twfe and aipw
  # attribute names for twfe
  est_name <- "alpha_weight"
  # catch if this aipw instead
  if ("att_weight" %in% attributes(decomposed_twfe_obj[[1]])$names) {
    est_name <- "att_weight"
  }
  
  group <- unlist(BMisc::getListElement(decomposed_twfe_obj, "g"))
  time.period <- unlist(BMisc::getListElement(decomposed_twfe_obj, "tp"))
  alpha_weight <- unlist(BMisc::getListElement(decomposed_twfe_obj, est_name))
  post <- 1*(time.period >= group)
  ncovs <- nrow(decomposed_twfe_obj[[1]]$cov_bal_df)
  covnames <- rownames(decomposed_twfe_obj[[1]]$cov_bal_df)
  cov_bal_df_list <- BMisc::getListElement(decomposed_twfe_obj, "cov_bal_df") 
  cov_bal_list <- lapply(1:ncovs, function(i) {
    do.call(rbind.data.frame, lapply(cov_bal_df_list, function(this_df) this_df[i, , drop=FALSE]))
  })
  unweighted_diffs <- unlist(lapply(cov_bal_list, function(x) sum(post*alpha_weight*x$unweighted_diff)))
  weighted_diffs <- unlist(lapply(cov_bal_list, function(x) sum(post*alpha_weight*x$weighted_diff)))
  wtd_avg_sd <- unlist(lapply(cov_bal_list, function(x) sum(post*alpha_weight*x$sd)))
  # sanity check: for twfe should balance across all periods, b/c time-invariant - ok
  # lapply(cov_bal_list, function(x) sum(alpha_weight*x$weighted_diff))
  
  cov_balance_summary_df <- cbind.data.frame(unweighted=unweighted_diffs, weighted=weighted_diffs, avg_sd=wtd_avg_sd)
  rownames(cov_balance_summary_df) <- covnames
  cov_balance_summary_df
}

# TODO: possibly add the outcome to the plot, to make it the same as for the two period 
# case, though it is not 100% obvious what the easiest way to do this is
#' @export
ggtwfeweights.decomposed_twfe <- function(decomposed_twfe_obj, standardize=TRUE) {
  
  # code to make this run for both twfe and aipw
  # catch if we are in aipw case
  if ("aipw_gt" %in% attributes(decomposed_twfe_obj)$names) {
    decomposed_twfe_obj <- decomposed_twfe_obj$aipw_gt
  } else {
    # otherwise we are in twfe
    decomposed_twfe_obj <- decomposed_twfe_obj$twfe_gt
  }
  
  if (is.null(decomposed_twfe_obj[[1]]$cov_bal_df)) {
    stop("nothing to plot.  Run twfe_cov_bal first")
  }
  cov_balance_df <- mp_covariate_bal_summary_helper(decomposed_twfe_obj)
  cov_balance_df$covariate <- rownames(cov_balance_df) %>% as.factor()
  if (standardize) {
    cov_balance_df$weighted <- cov_balance_df$weighted/cov_balance_df$avg_sd
    cov_balance_df$unweighted <- cov_balance_df$unweighted/cov_balance_df$avg_sd
  }
  ggplot(cov_balance_df, aes(y=.data$covariate, x=abs(.data$weighted))) + 
    geom_point(color="steelblue", size=3, shape=16) +
    geom_point(aes(x=abs(.data$unweighted)), color="red", size=3, shape=1) +
    theme_bw()
}

#' @export
implicit_aipw_weights <- function(yname,
                                  tname,
                                  idname,
                                  gname,
                                  xformula=~1,
                                  d_covs_formula=NULL,
                                  pscore_formula=xformula,
                                  pscore_d_covs_formula=d_covs_formula,
                                  extra_balance_vars_formula=NULL,
                                  extra_balance_d_vars_formula=NULL,
                                  balance_d_vars_post=TRUE,
                                  data,
                                  weightsname=NULL,
                                  ...) {
  
  # confirm balanced panel and sort data
  norig <- nrow(data) 
  data <- BMisc::makeBalancedPanel(data, idname=idname, tname=tname)
  data <- dplyr::arrange(data, idname, tname) %>% as.data.frame()
  if (nrow(data) != norig) {
    warning(paste0("Dropped ", norig-nrow(data), " rows of data to balance the data."))
  }
  
  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[,tname]))
  glist <- sort(unique(data[,gname]))
  data$.TP <- orig2t(data[,tname], tlist)
  data$.G <- orig2t(data[,gname], tlist)
  
  tlist <- sort(unique(data$.TP))
  glist <- sort(unique(data$.G))
  maxT <- max(tlist)
  if (!is.null(weightsname)) {
    sampling_weights <- data[,weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }
  gt_mat <- as.data.frame(expand.grid(glist[-1],tlist))
  gt_mat <- gt_mat %>% arrange_all() %>% as.data.frame()
  colnames(gt_mat) <- c("group","time.period")
  # drop pre-treatment values of t
  gt_mat <- gt_mat[gt_mat$time.period >= gt_mat$group,]
  
  aipw_gt <- lapply(1:nrow(gt_mat), function(i) {
    this_data <- data[data$.G %in% c(0,gt_mat[i,1]) & data$.TP %in% c(gt_mat[i,1]-1, gt_mat[i,2]),]
    two_period_aipw_res <- two_period_aipw_weights(yname=yname,
                                                   tname=".TP",
                                                   idname=idname,
                                                   gname=".G",
                                                   xformula=xformula,
                                                   d_covs_formula=d_covs_formula,
                                                   pscore_formula=pscore_formula,
                                                   pscore_d_covs_formula=pscore_d_covs_formula,
                                                   extra_balance_vars_formula=extra_balance_vars_formula,
                                                   extra_balance_d_vars_formula=extra_balance_d_vars_formula,
                                                   balance_d_vars_post=balance_d_vars_post,
                                                   data=this_data,
                                                   weightsname=weightsname)
    # drop local covariate balance, as it is confusing in multiple period case
    two_period_aipw_res$cov_balance_df <- NULL
    two_period_aipw_res$g <- gt_mat[i,1]
    two_period_aipw_res$tp <- gt_mat[i,2]
    barpg <- weighted.mean(data[data$.G!=0,]$.G==gt_mat[i,1], w=sampling_weights[data$.G!=0])
    two_period_aipw_res$att_weight <- barpg / (maxT - gt_mat[i,1] + 1)
    two_period_aipw_res
  })  
  
  out <- list()
  out$aipw_gt <- aipw_gt
  out$idname <- idname
  out$tname <- tname
  out$gname <- gname
  out$weightsname <- weightsname
  out$data <- data
  class(out) <- "decomposed_aipw"
  out
}


aipw_cov_bal <- function(decomposed_aipw, extra_balance_vars_formula) {
  idname <- decomposed_aipw$idname
  tname <- decomposed_aipw$tname
  gname <- decomposed_aipw$gname
  data <- decomposed_aipw$data
  weightsname <- decomposed_aipw$weightsname
  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[,tname]))
  glist <- sort(unique(data[,gname]))
  data$.TP <- orig2t(data[,tname], tlist)
  data$.G <- orig2t(data[,gname], tlist)
  tlist <- sort(unique(data$.TP))
  glist <- sort(unique(data$.G))
  
  # construct a data.frame that has the mean of each extra balance variable
  extra_vars_frame <- as.data.frame(model.matrix(BMisc::addCovToFormla("-1", extra_balance_vars_formula), data=data))
  mean_extra_var_list <- lapply(1:ncol(extra_vars_frame), function(j) {
    this_colname <- colnames(extra_vars_frame)[j]
    this_frame <- cbind.data.frame(extra_vars_frame[,j], data[,c(idname,tname,gname)])
    colnames(this_frame)[1] <- this_colname
    BMisc::get_Yibar(this_frame, idname=idname, yname=this_colname)
  })  
  mean_extra_vars <- do.call(cbind.data.frame, mean_extra_var_list)
  colnames(mean_extra_vars) <- paste0("mean_", colnames(extra_vars_frame))
  
  # compute balance statistics for each extra balance variable
  groups <- unlist(BMisc::getListElement(decomposed_aipw$aipw_gt, "g"))
  time_periods <- unlist(BMisc::getListElement(decomposed_aipw$aipw_gt, "tp"))
  gt_mat <- cbind(groups, time_periods)
  
  for (i in 1:nrow(gt_mat)) {
    this_balance_df <- aipw_cov_bal_gt(g = gt_mat[i, 1],
                                       tp = gt_mat[i, 2],
                                       decomposed_aipw = decomposed_aipw,
                                       X = mean_extra_vars)
    decomposed_aipw$aipw_gt[[i]]$cov_bal_df <- this_balance_df
  }
  
  # return the decomposed_twfe object with the balance statistics added
  decomposed_aipw
}

#' @export
aipw_cov_bal_gt <- function(g, tp, decomposed_aipw, X) {
  idname <- decomposed_aipw$idname
  tname <- ".TP"
  gname <- ".G"
  data <- decomposed_aipw$data
  weightsname <- decomposed_aipw$weightsname
  groups <- unlist(BMisc::getListElement(decomposed_aipw$aipw_gt, "g"))
  time_periods <- unlist(BMisc::getListElement(decomposed_aipw$aipw_gt, "tp"))
  this_gt_idx <- which(groups==g & time_periods==tp)
  this_decomposed_aipw <- decomposed_aipw$aipw_gt[[this_gt_idx]]
  dy <- this_decomposed_aipw$dy
  D <- this_decomposed_aipw$D
  weights_treated <- this_decomposed_aipw$weights[D==1]
  weights_comparison <- this_decomposed_aipw$weights[D==0]
  # recovers estimate of ATT(g,t)
  # mean(weights_treated*dy[D==1]) - 
  #  mean(weights_untreated*dy[D==0])
  treated_idx <- data$.G==g & data$.TP==tp
  comparison_idx <- data$.G==0 & data$.TP==tp
  
  # sampling weights are already accounted for in this_decomposed_aipw$weights
  # so don't use them in weighted results below (but do use them for unweighted)
  if(!is.null(weightsname)) {
    sampling_weights <- data[,weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }
  
  unweighted_covs_treated <- apply(X[treated_idx, , drop=FALSE], 2, function(x) {
    weighted.mean(x, w=sampling_weights[treated_idx])
  })
  unweighted_covs_comparison <- apply(X[comparison_idx, , drop=FALSE], 2, function(x) {
    weighted.mean(x, w=sampling_weights[comparison_idx])
  })
  unweighted_diff <- unweighted_covs_treated - unweighted_covs_comparison
  
  weighted_covs_treated <- apply(X[treated_idx, , drop=FALSE], 2, function(x) {
    mean(x*weights_treated)
  })
  weighted_covs_comparison <- apply(X[comparison_idx, , drop=FALSE], 2, function(x) {
    mean(x*weights_comparison)
  })
  weighted_diff <- weighted_covs_treated - weighted_covs_comparison
  
  both_idx <- treated_idx | comparison_idx
  sd <- apply(X, 2, function(x) {
    sqrt(weighted.mean((x[both_idx] - weighted.mean(x[both_idx], w=sampling_weights[both_idx]))^2, w=sampling_weights[both_idx]))
  })  
  
  # return data frame that contains the weighted covariates for the treated and comparison groups
  cbind.data.frame(unweighted_covs_treated,
                   unweighted_covs_comparison,
                   unweighted_diff,
                   weighted_covs_treated, 
                   weighted_covs_comparison,
                   weighted_diff,
                   sd)
}


#'@export
summary.decomposed_aipw <- function(object, ...) {
  object <- object$aipw_gt
  aipw_att_gt <- unlist(BMisc::getListElement(object, "est"))
  group <- unlist(BMisc::getListElement(object, "g"))
  time.period <- unlist(BMisc::getListElement(object, "tp"))
  att_weight <- unlist(BMisc::getListElement(object, "att_weight"))
  post <- 1*(time.period >= group)
  print_df <- cbind.data.frame(group, time.period, post, aipw_att_gt, att_weight)
  atto <- sum(att_weight*aipw_att_gt)
  cat("overall att:          ", round(atto,4), "\n")
  cat("\n")
  # this is basically just a check to see if the balance statistics have been computed
  # if not, we'll print the summary of the decomposition with respect 
  # to the actual outcome.  Otherwiser, we will summarize the covariate balance
  if (is.null(object[[1]]$cov_bal_df)) {
    cat("summary of aipw decomposition\n")
    print_df
  } else {
    cat("summary of aipw covariate balance:\n")
    cov_balance_summary_df <- mp_covariate_bal_summary_helper(object)
    cov_balance_summary_df
  }
}

#' @export
ggtwfeweights.decomposed_aipw <- function(decomposed_aipw_obj) {
  # should work to just call the plot for twfe at this point
  ggtwfeweights.decomposed_twfe(decomposed_aipw_obj)
}