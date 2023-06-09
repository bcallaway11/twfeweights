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
#' @export
gt_weights <- function(g,
                       tp,
                       treated=NULL,
                       comparison=NULL,
                       weights_treated,
                       weights_comparison,
                       weighted_outcome_treated=NULL,
                       weighted_outcome_comparison=NULL,
                       weighted_outcome_diff=NULL) {
  
  ret_gt_weights <- list(g=g,
                         tp=tp,
                         treated=treated,
                         comparison=comparison,
                         weights_treated=weights_treated,
                         weights_comparison=weights_comparison,
                         weighted_outcome_treated=weighted_outcome_treated,
                         weighted_outcome_comparison=weighted_outcome_comparison,
                         weighted_outcome_diff=weighted_outcome_diff)
  class(ret_gt_weights) <- "gt_weights"
  ret_gt_weights


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
#' @param xname by default, `twfe_weight_gt` applies the 
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
twfe_weights_gt <- function(g,
                           tp,
                           yname,
                           tname,
                           idname,
                           gname,
                           xformla=NULL,
                           xname=NULL,
                           data) {
  
  G <- data[,gname]
  TP <- data[,tname]
  Yit <- data[,yname]
  if (!is.null(xname)) Xi <- data[,xname]
  minT <- min(unique(TP))
  nT <- length(unique(TP))
  
  idx <- G==g & TP==tp
  idx_preg <- G==g & TP==(g-1)
  pg <- mean(G==g)
  pbarg <- mean( G[G!=0] == g)
  pu <- mean(G==0)
  time_weight <- (nT-(g-minT+1)+1)/nT
  pg_weight <- pg/pbarg
  if (is.null(xname)) {
    term_1 <- mean( (ddotDit[idx] - ddotXit[idx]*gam)*(Yit[idx] - Yit[idx_preg])*time_weight*pg_weight/alp_den )
  } else {
    term_1 <- mean( (ddotDit[idx] - ddotXit[idx]*gam)*Xi[idx]*time_weight*pg_weight/alp_den )
  }
  idx2 <- G==0 & TP==tp
  idx2_preg <- G==0 & TP==(g-1)
  if (is.null(xname)) {
    term_2 <- mean( (ddotDit[idx2] - ddotXit[idx2]*gam)*(Yit[idx2] - Yit[idx2_preg])*time_weight*pu/alp_den )
  } else {
    term_2 <- mean( (ddotDit[idx2] - ddotXit[idx2]*gam)*Xi[idx2]*time_weight*pu/alp_den )
  }
  term_1 + term_2
}

df$pop2006 <- get_Yi1(df, "id", "population", "year", "group")

years <- c(2006, 2007, 2008)
gt_mat <- expand.grid(years, years[-1])
gt_mat <- cbind.data.frame(gt_mat[,2], gt_mat[,1])
colnames(gt_mat) <- c("group","time.period")

twfe_gt <- sapply(1:nrow(gt_mat), function(i) {
  twfe_weight_gt(g=gt_mat[i,1],
                 tp=gt_mat[i,2],
                 yname="homicide_c",
                 tname="year",
                 idname="id",
                 gname="group",
                 xformla=NULL,
                 data=df)
})

combine_twfe_weights_gt <- function(g,tp,gname,data) {
  G <- data[,gname]
  pbarg <- mean( G[G!=0] == g)
  ntreatedperiods <- (nT-(g-minT+1)+1)
  pbarg/ntreatedperiods
}


twfe_gt_weight <- sapply(1:nrow(gt_mat), function(i) {
  combine_twfe_weights_gt(g=gt_mat[i,1],
                          tp=gt_mat[i,2],
                          gname="group",
                          data=df)
})
