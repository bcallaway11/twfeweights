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

  cbind.data.frame(G=.group, TP=.t, wTWFEgt, attgt=.attgt)
  
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
  
  cbind.data.frame(G=.group, TP=.t, wOgt, attgt=.attgt)
}
