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
                            xname=NULL,
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
  if (!is.null(xname)) Xi <- data[,xname]

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
                             xname=NULL,
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
                             xname=xname,
                             data=data,
                             weightsname=weightsname)
  })  
  
  class(twfe_gt) <- "decomposed_twfe"
  twfe_gt
}

#'@title summary.decomposed_twfe
#'@description summarizes a two-way fixed effects regression 
#' decomposition
#'@param object a decomposed_twfe object
#'@param ... extra arguments, not used here
#'
#'@export
summary.decomposed_twfe <- function(object, ...) {
  twfe_att_gt <- unlist(BMisc::getListElement(object, "weighted_outcome_diff"))
  group <- unlist(BMisc::getListElement(object, "g"))
  time.period <- unlist(BMisc::getListElement(object, "tp"))
  alpha_weight <- unlist(BMisc::getListElement(object, "alpha_weight"))
  post <- 1*(time.period >= group)
  print_df <- cbind.data.frame(group, time.period, post, twfe_att_gt, alpha_weight)
  alpha <- sum(alpha_weight*twfe_att_gt)
  pt_violations_bias <- sum(alpha_weight[post==0]*twfe_att_gt[post==0])
  cat("alpha twfe:          ", round(alpha,4), "\n")
  cat("pta violations bias: ", round(pt_violations_bias,4),"\n")
  print_df
}



# ------------------------------------------------------------------------------
# old, delete at final version
# the code below if for the case where we use the same outside weights as for 
# ATT^o, but the drawback is that the inside weights do not have mean 1
# ------------------------------------------------------------------------------
# 
# twfe_weights_gt <- function(g,
#                             tp,
#                             yname,
#                             tname,
#                             idname,
#                             gname,
#                             xformula=NULL,
#                             xname=NULL,
#                             data) {
#   
#   # sort data
#   data <- dplyr::arrange(data, idname, tname)
#   
#   # setup 
#   G <- data[,gname]
#   TP <- data[,tname]
#   tlist <- sort(unique(TP))
#   preg <- tlist[tail(which(tlist < g),1)] # universal base period
#   Yit <- data[,yname]
#   xformula <- BMisc::addCovToFormla("-1", xformula)
#   Xit <- model.matrix(xformula, data=data)
#   Dit <- 1*( (data[,tname] >= data[,gname]) & (data[,gname] != 0) )
#   if (!is.null(xname)) Xi <- data[,xname]
#   minT <- min(unique(TP))
#   nT <- length(unique(TP))
#   
#   # double demean variables in the model
#   ddotDit <- fixest::demean(Dit, data[,c(idname,tname)]) %>% as.numeric
#   ddotYit <- fixest::demean(Yit, data[,c(idname,tname)]) %>% as.numeric
#   ddotXit <- fixest::demean(Xit, data[,c(idname,tname)])
#   
#   # weights depend on linear projection of treatment 
#   # on covariates using all periods
#   lp_treat <- lm(ddotDit ~ -1 + ddotXit)
#   gam <- as.matrix(coef(lp_treat))
#   
#   # pick the right periods and compute the weights
#   idx <- G==g & TP==tp
#   idx_preg <- G==g & TP==preg
#   pg <- mean(G==g)
#   pbarg <- mean( G[G!=0] == g)
#   pu <- mean(G==0)
#   local_g <- which(tlist == g)
#   time_weight <- (nT-local_g+1)/nT
#   pg_weight <- pg/pbarg
#   linear_pscore <- as.numeric(ddotXit%*%gam)
#   alp_den <- mean((ddotDit - linear_pscore)*(ddotDit))
#   weights_1 <- (ddotDit[idx] - linear_pscore[idx])*time_weight*pg_weight/alp_den
#   if (is.null(xname)) {
#     this_treated <- Yit[idx] - Yit[idx_preg]
#   } else {
#     this_treated <- Xi[idx]
#   }
#   term_1 <- mean( weights_1 * this_treated)
#   
#   idx2 <- G==0 & TP==tp
#   idx2_preg <- G==0 & TP==preg
#   weights_2 <- (ddotDit[idx2] - linear_pscore[idx2])*time_weight*pu/alp_den
#   if (is.null(xname)) {
#     this_comparison <- Yit[idx2] - Yit[idx2_preg]
#   } else {
#     this_comparison <- Xi[idx2]
#   }
#   term_2 <- mean(weights_2 * this_comparison)
#   
#   alpha_weight <- combine_twfe_weights_gt(g,tp,gname,tname,data)
#   
#   # make things negative if in pre-treatment period
#   if (tp < g) {
#     term_1 <- -term_1
#     term_2 <- -term_2
#   }
#   
#   ## browser()
#   ## # some debugging code
#   ## w1 <- ddotDit[idx] - linear_pscore[idx]
#   ## w2 <- ddotDit[idx2] - linear_pscore[idx2]
#   ## m1<-mean(w1)
#   ## m2<-mean(w2)
#   ## m1*pg/pbarg
#   ## m2*pu
#   ## idx3 <- (G==0 | tp < G) & TP==tp
#   ## idx3_preg <- (G==0 | tp < G) & TP==preg
#   ## w3 <- ddotDit[idx3] - linear_pscore[idx3]
#   ## m3 <- mean(w3)
#   ## pdt0 <- mean(G==0 | tp < G)
#   ## m1
#   ## m2*(pdt0/(1-pdt0))
#   
#   out <- gt_weights(g=g,
#                     tp=tp,
#                     treated=this_treated,
#                     comparison=this_comparison,
#                     weights_treated=weights_1,
#                     weights_comparison=-weights_2,
#                     weighted_outcome_treated=term_1,
#                     weighted_outcome_comparison=term_2,
#                     weighted_outcome_diff=(term_1 + term_2),
#                     alpha_weight=alpha_weight)
#   out
# }

# combine_twfe_weights_gt <- function(g,
#                                     tp,
#                                     gname,
#                                     tname,
#                                     data) {
#   G <- data[,gname]
#   TP <- data[,tname]
#   tlist <- sort(unique(TP))
#   local_g <- which(tlist == g) # the group in terms of t=1,2,3,...,T
#   minT <- min(unique(TP))
#   nT <- length(unique(TP))
#   pbarg <- mean( G[G!=0] == g)
#   ntreatedperiods <- nT - local_g + 1
#   w_gt <- pbarg/ntreatedperiods
#   if (tp < g) w_gt <- -w_gt
#   w_gt
# }