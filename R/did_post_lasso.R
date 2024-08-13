#' @title did_post_lasso
#' @description Function to estimate the ATT using post-Lasso and compute the asymptotic variance with normalized weights
#'
#' @param y1 Numeric vector of outcomes in the second time period
#' @param y0 Numeric vector of outcomes in the first time period
#' @param D Numeric vector of treatment indicators
#' @param covariates Matrix of covariates
#' @param i.weights Numeric vector of sampling weights, defaault is NULL, so that observations
#'  are equally weighted
#' @param inffunc Logical indicating whether to compute the influence function, default is TRUE
#' @param nfolds Number of folds for cross-validation, default is 5
#' @param alpha Elastic net mixing parameter, default is 1 for Lasso
#'
#' @return A list with the ATT estimate, standard error, selected variables for the outcome model,
#' selected variables for the propensity score model, and the influence function
#' @export
did_post_lasso <- function(
    y1,
    y0,
    D,
    covariates,
    pscore_covariates = covariates,
    i.weights = NULL,
    inffunc = TRUE,
    nfolds = 5) {
  # Calculate outcome
  outcome <- y1 - y0
  n <- length(outcome)
  X <- as.matrix(covariates)
  if (all(covariates[, 1] == 1)) {
    X <- X[, -1] # drop intercept, if it is there
  }
  folds <- sample(rep(1:nfolds, length.out = n))

  if (is.null(i.weights)) {
    i.weights <- rep(1, n) # Default weights to 1 if not provided
  }

  psi <- numeric(n)

  # Calculate outcome
  outcome <- y1 - y0
  n <- length(outcome)
  X <- as.matrix(covariates)
  if (all(covariates[, 1] == 1)) {
    X <- X[, -1] # drop intercept, if it is there
  }

  if (is.null(i.weights)) {
    i.weights <- rep(1, n) # Default weights to 1 if not provided
  }

  psi <- numeric(n)

  # Post-Lasso estimation for m_0(X)
  cv_fit_outcome <- cv.glmnet(X[D == 0, ],
    outcome[train_indices & D == 0],
    alpha = 1,
    weights = i.weights[D == 0]
  )
  lambda_outcome <- cv_fit_outcome$lambda.min
  lasso_outcome <- glmnet(X[D == 0, ],
    outcome[D == 0],
    alpha = 1,
    lambda = lambda_outcome,
    weights = i.weights[D == 0]
  )

  browser()

  selected_vars_outcome <- which(coef(lasso_outcome) != 0)[-1] - 1 # drop intercept

  if (ncol(pscore_covariates) > 1) { # if there are any covariates in the pscore model, then...

    for (k in 1:nfolds) {
      # Split data into folds
      train_indices <- folds != k
      test_indices <- folds == k

      # Post-Lasso estimation for p(X)
      cv_fit_propensity <- cv.glmnet(X[train_indices, ],
        D[train_indices],
        family = "binomial",
        alpha = 1,
        weights = i.weights[train_indices]
      )
      lambda_propensity <- cv_fit_propensity$lambda.min
      lasso_propensity <- glmnet(X[train_indices, ],
        D[train_indices],
        family = "binomial",
        alpha = 1,
        lambda = lambda_propensity,
        weights = i.weights[train_indices]
      )
      selected_vars_propensity <- which(coef(lasso_propensity) != 0)[-1] - 1 # drop intercept

      pscore_model <- glm(D ~ .,
        data = cbind.data.frame(
          D = D[train_indices],
          X[train_indices, selected_vars_propensity, drop = FALSE]
        ),
        family = "binomial", weights = i.weights[train_indices]
      )
      pscore <- predict(pscore_model, newdata = as.data.frame(X[test_indices, ]), type = "response")

      # Compute raw weights
      treated_raw_weights <- D[test_indices]
      treated_normalized_weights <- treated_raw_weights / mean(treated_raw_weights)
      untreated_raw_weights <- (1 - D[test_indices]) * pscore / (1 - pscore)
      untreated_normalized_weights <- untreated_raw_weights / mean(untreated_raw_weights)

      # Normalize the weights
      normalized_weights <- treated_normalized_weights - untreated_normalized_weights

      # Compute influence function for test fold
      for (i in seq_along(test_indices)) {
        index <- test_indices[i]
        psi[index] <- normalized_weights[i] * (outcome[index] - m0_hat[i])
      }
    }
  }

  # Compute ATT estimate
  att_hat <- sum(psi * i.weights) / sum(i.weights)

  # Compute the asymptotic variance
  var_att <- sum((psi^2) * i.weights) / (sum(i.weights)^2)
  se_att <- sqrt(var_att)


  # Return results
  list(
    ATT = att_hat,
    se = se_att,
    selected_vars_outcome = selected_vars_outcome,
    selected_vars_propensity = selected_vars_propensity,
    att.inf.func = psi
  )
}


#' @title did_post_lasso_ra
#' @description Function to estimate the ATT using post-Lasso and compute the asymptotic variance with normalized weights
#'  that only uses regression adjustment.  This function requires stronger assumptions on sparsity that `did_post_lasso`,
#'  but it does not require the estimation of the propensity score, which may be infeasible in applications with
#'  staggered treatment adoption and where the number of units per timing group is fairly small.
did_post_lasso_ra <- function(
    y1,
    y0,
    D,
    covariates,
    i.weights = NULL,
    inffunc = TRUE,
    alpha = 1) {
  browser()

  # effective_covariates <- cbind(1, select)
  # DRDID::reg_did_panel(y1, y0,  D, inffunc=TRUE)


  # Compute ATT estimate
  att_hat <- sum(psi * i.weights) / sum(i.weights)

  # Compute the asymptotic variance
  var_att <- sum((psi^2) * i.weights) / (sum(i.weights)^2)
  se_att <- sqrt(var_att)


  # Return results
  list(
    ATT = att_hat,
    se = se_att,
    selected_vars_outcome = selected_vars_outcome,
    att.inf.func = psi
  )
}
