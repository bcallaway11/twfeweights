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