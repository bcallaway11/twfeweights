#--------------------------------------------------------------------------
#
# All code in this file is for the case with no covariates.  
#
# The TWFE are essentially those from de Chaisemartin and d'Haultfoeuille (2020)
# though they emphasize decompositions in the sense that they do not rely 
# on imposing parallel trends and (i) the results add up to be exactly equal 
# to the TWFE regression coefficient (ii) violations of parallel trends 
# in pre-treatment contribute to the reported results.
#
#--------------------------------------------------------------------------

#' @title mp_weights_obj
#' 
#' @description A class to store weighted averages of group-time average 
#'  treatment effects
#'  
#' @param weights_df A data.frame containing the following elements:
#'  - group: group identifier
#'  - time.period: time period identifier
#'  - weight: weight for the group-time period
#'  - attgt: average treatment effect for the group-time period
#'  
#' @export
mp_weights_obj <- function(weights_df) {
  weights_df$post <- as.factor(1*(weights_df$time.period >= weights_df$group))
  out <- list()
  out$weights_df <- weights_df
  class(out) <- "mp_weights_obj"
  out
}


#' @title twfe_weights
#' @description A function to compute TWFE regression weights on ATT(g,t)
#'  in a staggered treatment adoption setting.
#' @param attgt_res Results from did::att_gt function
#'
#' @return a data.frame containing the TWFE weighs
#' @export
twfe_weights <- function(attgt_res) {

  if (attgt_res$DIDparams$control_group != "nevertreated") {
    stop("only \"nevertreated\" control group is supported.")
  }

  if (attgt_res$DIDparams$base_period != "universal") {
    stop("only \"universal\" base period is supported.")
  }

  if (attgt_res$DIDparams$xformla != ~1) {
    stop("including covariates is not supported.")
  }

  .gname <- attgt_res$DIDparams$gname
  .tname <- attgt_res$DIDparams$tname
  .data <- attgt_res$DIDparams$data
  .tlist <- attgt_res$DIDparams$tlist
  .glist <- c(0,attgt_res$DIDparams$glist)
  .attgt <- c(rep(0, length(.tlist)), attgt_res$att)
  .t <- c(.tlist, attgt_res$t)
  .group <- c(rep(0, length(.tlist)), attgt_res$group)
  
  # twfe weights
  Edt <- function(t) {
    mean( (t >= .data[,.gname]) & (.data[,.gname] !=0) )
  }
  mEdt <- mean(sapply(.tlist, Edt))
  pg2 <- sapply(.glist, function(g) mean(.data[,.gname]==g))
  maxT <- max(.tlist)
  wTWFE_num <- function(g,t,pre0=FALSE) {
    if ((t < g) & pre0) return(0) 
    if (g==0) {
      hgt <- -Edt(t) + mEdt
    } else {
      hgt <- 1*(t>=g) - (maxT-g+1)/length(.tlist) - Edt(t) + mEdt
    }
    hgt*pg2[.glist==g]
  }

  # calculate weights
  wTWFEgt_num <- sapply(1:length(.attgt), 
                        function(i) wTWFE_num(.group[i],
                                              .t[i],
                                              pre0=FALSE))
  wTWFEgt_den <- sapply(1:length(.attgt), 
                        function(i) wTWFE_num(.group[i],
                                              .t[i],
                                              pre0=FALSE))

  cond <- (.t >= .group & .group != 0)
  wTWFEgt <- wTWFEgt_num/sum(wTWFEgt_den[cond])

  out <- mp_weights_obj(cbind.data.frame(group=.group, time.period=.t, weight=wwTWFEgt, attgt=.attgt))
}

#' @title attO_weights
#' @description A function to compute weights on ATT(g,t) to deliver
#'  ATT^O as discussed in Callaway and Sant'Anna (2021)
#' @inheritParams twfe_weights
#' @param w Optional external unit-specific weights
#'
#' @return a data.frame containing attO weights
#' @export
attO_weights <- function(attgt_res, w=rep(1,nrow(attgt_res$DIDparams$data))) {

  .gname <- attgt_res$DIDparams$gname
  .tname <- attgt_res$DIDparams$tname
  .data <- attgt_res$DIDparams$data
  .tlist <- attgt_res$DIDparams$tlist
  .glist <- c(0,attgt_res$DIDparams$glist)
  .attgt <- c(rep(0, length(.tlist)), attgt_res$att)
  .t <- c(.tlist, attgt_res$t)
  .group <- c(rep(0, length(.tlist)), attgt_res$group)

  # overall attgt weights
  .ever_treated <- which(.data[,.gname] != 0)
  w <- w[.ever_treated]
  pg <- sapply(.glist, function(g) weighted.mean(.data[.ever_treated,][,.gname]==g, w=w))
  maxT <- max(.tlist)
  wO <- function(g,t) {
    1*(t >= g)*pg[.glist==g] / (maxT - g + 1)
  }
  # add weights to results
  wOgt <- sapply(1:length(.attgt),
                 function(i) wO(.group[i],
                                .t[i]))
  
  cbind.data.frame(group=.group, time.period=.t, weight=wOgt, attgt=.attgt)
}


#' @title att_simple_weights
#' @description A function to compute weights on ATT(g,t) to deliver
#'  ATT^{simple} as discussed in Callaway and Sant'Anna (2021)
#' @inheritParams twfe_weights
#' @param w Optional external unit-specific weights
#'
#' @return a data.frame containing att_simple weights
#' @export
att_simple_weights <- function(attgt_res, w=rep(1,nrow(attgt_res$DIDparams$data))) {
  
  .gname <- attgt_res$DIDparams$gname
  .tname <- attgt_res$DIDparams$tname
  .data <- attgt_res$DIDparams$data
  .tlist <- attgt_res$DIDparams$tlist
  .glist <- c(0,attgt_res$DIDparams$glist)
  .attgt <- c(rep(0, length(.tlist)), attgt_res$att)
  .t <- c(.tlist, attgt_res$t)
  .group <- c(rep(0, length(.tlist)), attgt_res$group)
  
  # overall attgt weights
  .ever_treated <- which(.data[,.gname] != 0)
  w <- w[.ever_treated]
  pg <- sapply(.glist, function(g) weighted.mean(.data[.ever_treated,][,.gname]==g, w=w))
  maxT <- max(.tlist)
  wO <- function(g,t) {
    1*(t >= g)*pg[.glist==g]
  }
  # add weights to results
  wsimplegt <- sapply(1:length(.attgt),
                 function(i) wO(.group[i],
                                .t[i]))
  # account for denominator just by normalizing weights to sum to 1
  wsimplegt <- wsimplegt / sum(wsimplegt) 
  
  cbind.data.frame(group=.group, time.period=.t, weight=wsimplegt, attgt=.attgt)
}

#' @title ggtwfeweights.mp_weights_obj
#' 
#' @description A function to plot weights on group-time average treatment effects
#' 
#' @param x An object of class mp_weights_obj
#' @param ... Additional arguments to be passed to ggplot
#' 
#' @export
ggtwfeweights.mp_weights_obj <- function(x, ...) {
  weights_df <- x$weights_df
  ggplot(data=weights_df,
         mapping=aes(x=weight, y=attgt, color=post)) +
    geom_hline(yintercept=0, linewidth=1.5) +
    geom_vline(xintercept=0, linewidth=1.5) + 
    geom_point(size=6) +
    theme_bw()
}