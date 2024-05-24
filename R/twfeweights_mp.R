# Code for computing TWFE weights with multiple
# periods and variation in treatment timing

#' @title gt_weights
#' @description a class for holding and working with group-time weights
#' @inheritParams twfe_weights_gt
#' @param treated a vector of units that make up the effective treated group for group g 
#'  and period tp
#' @param comparison a vector of units that make up the effective comparison 
#'  group for group g and period tp
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
                         alpha_weight=alpha_weight)
  class(ret_gt_weights) <- "gt_weights"
  ret_gt_weights
}

#' @title twfe_weights_gt
#' @description Compute TWFE weights for a particular group
#'  and time period.  This function is developed for the 
#'  case where there is a binary treatment and staggered 
#'  treatment adoption.  The weights are for the coefficient 
#'  on the binary treatment variable and allows for 
#'  the regression to include additional time varying 
#'  covariates.  
#'
#' @param g the current group for which the weights will be computed
#' @param tp the time period for which the weights will be computed
#' @inheritParams did::att_gt
#' @param xformla formula for covariates that are included in the 
#'  TWFE regression.
#' @param xname by default, `twfe_weights_gt` applies the 
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
twfe_weights_gt <- function(g,
                           tp,
                           yname,
                           tname,
                           idname,
                           gname,
                           xformla=NULL,
                           xname=NULL,
                           data) {
  
  # confirm balanced panel and sort data
  data <- BMisc::makeBalancedPanel(data, idname=idname, tname=tname)
  
  data <- dplyr::arrange(data, idname, tname)

  # setup 
  G <- data[,gname]
  TP <- data[,tname]
  tlist <- sort(unique(TP))
  preg <- tlist[tail(which(tlist < g),1)] # universal base period
  Yit <- data[,yname]
  xformla <- BMisc::addCovToFormla("-1", xformla)
  Xit <- model.matrix(xformla, data=data)
  Dit <- 1*( (data[,tname] >= data[,gname]) & (data[,gname] != 0) )
  if (!is.null(xname)) Xi <- data[,xname]
  minT <- min(unique(TP))
  nT <- length(unique(TP))
  
  # double demean variables in the model
  ddotDit <- fixest::demean(Dit, data[,c(idname,tname)]) %>% as.numeric
  ddotYit <- fixest::demean(Yit, data[,c(idname,tname)]) %>% as.numeric
  ddotXit <- fixest::demean(Xit, data[,c(idname,tname)])
  
  # weights depend on linear projection of treatment 
  # on covariates using all periods
  lp_treat <- lm(ddotDit ~ -1 + ddotXit)
  gam <- as.matrix(coef(lp_treat))
  
  # pick the right periods and compute the weights
  idx <- G==g & TP==tp
  idx_preg <- G==g & TP==preg
  pg <- mean(G==g)
  pbarg <- mean( G[G!=0] == g)
  pu <- mean(G==0)
  local_g <- which(tlist == g)
  time_weight <- (nT-local_g+1)/nT
  pg_weight <- pg/pbarg
  linear_pscore <- as.numeric(ddotXit%*%gam)
  alp_den <- mean((ddotDit - linear_pscore)*(ddotDit))
  weights_1 <- (ddotDit[idx] - linear_pscore[idx])*time_weight*pg_weight/alp_den
  if (is.null(xname)) {
    this_treated <- Yit[idx] - Yit[idx_preg]
  } else {
    this_treated <- Xi[idx]
  }
  term_1 <- mean( weights_1 * this_treated)
  
  idx2 <- G==0 & TP==tp
  idx2_preg <- G==0 & TP==preg
  weights_2 <- (ddotDit[idx2] - linear_pscore[idx2])*time_weight*pu/alp_den
  if (is.null(xname)) {
    this_comparison <- Yit[idx2] - Yit[idx2_preg]
  } else {
    this_comparison <- Xi[idx2]
  }
  term_2 <- mean(weights_2 * this_comparison)
  
  alpha_weight <- combine_twfe_weights_gt(g,tp,gname,tname,data)

  # make things negative if in pre-treatment period
  if (tp < g) {
    term_1 <- -term_1
    term_2 <- -term_2
  }

  ## browser()
  ## # some debugging code
  ## w1 <- ddotDit[idx] - linear_pscore[idx]
  ## w2 <- ddotDit[idx2] - linear_pscore[idx2]
  ## m1<-mean(w1)
  ## m2<-mean(w2)
  ## m1*pg/pbarg
  ## m2*pu
  ## idx3 <- (G==0 | tp < G) & TP==tp
  ## idx3_preg <- (G==0 | tp < G) & TP==preg
  ## w3 <- ddotDit[idx3] - linear_pscore[idx3]
  ## m3 <- mean(w3)
  ## pdt0 <- mean(G==0 | tp < G)
  ## m1
  ## m2*(pdt0/(1-pdt0))
  
  out <- gt_weights(g=g,
                    tp=tp,
                    treated=this_treated,
                    comparison=this_comparison,
                    weights_treated=weights_1,
                    weights_comparison=-weights_2,
                    weighted_outcome_treated=term_1,
                    weighted_outcome_comparison=term_2,
                    weighted_outcome_diff=(term_1 + term_2),
                    alpha_weight=alpha_weight)
  out
}

#' @title combine_twfe_weights_gt
#' @description a function that recovers the weights on each implicit 
#'  group-time average treatment effect from the TWFE regression
#' @inheritParams twfe_weights_gt
#' @export
combine_twfe_weights_gt <- function(g,
                                    tp,
                                    gname,
                                    tname,
                                    data) {
  G <- data[,gname]
  TP <- data[,tname]
  tlist <- sort(unique(TP))
  local_g <- which(tlist == g) # the group in terms of t=1,2,3,...,T
  minT <- min(unique(TP))
  nT <- length(unique(TP))
  pbarg <- mean( G[G!=0] == g)
  ntreatedperiods <- nT - local_g + 1
  w_gt <- pbarg/ntreatedperiods
  if (tp < g) w_gt <- -w_gt
  w_gt
}

#' @title all_twfe_weights 
#' @description a function to compute implicit TWFE weights 
#' 
#' @inheritParams twfe_weights_gt
#' @return a `decomposed_twfe` object
#' @export
all_twfe_weights <- function(yname,
                             tname,
                             idname,
                             gname,
                             xformla=NULL,
                             xname=NULL,
                             data) {
  tlist <- sort(unique(data[,tname]))
  glist <- sort(unique(data[,gname]))
  glist <- glist[glist!=0] # drop never-treated group
  gt_mat <- as.data.frame(expand.grid(glist,tlist))
  gt_mat <- gt_mat %>% arrange_all()
  colnames(gt_mat) <- c("group","time.period")
  
  twfe_gt <- lapply(1:nrow(gt_mat), function(i) {
    twfe_weights_gt(g=gt_mat[i,1],
                   tp=gt_mat[i,2],
                   yname=yname,
                   tname=tname,
                   idname=idname,
                   gname=gname,
                   xformla=xformla,
                   xname=xname,
                   data=df)
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

