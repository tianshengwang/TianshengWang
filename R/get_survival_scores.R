#' Compute doubly robust survival scores for iCSF algorithm.
#' Author: Erik Sverdrup
#' Compute doubly robust scores for 1) survival probability estimation \eqn{P[T_i > \text{horizon}]},
#' or 2) restricted mean survival time (RMST) \eqn{E[min(T_i, \text{horizon})]} under discrete-time right-censoring.
#'
#' Let \eqn{T_i} be the survival time for unit \eqn{i}. Due to censoring we only get to observe
#' \eqn{Y_i = min(T_i, C_i)} along with an event indicator \eqn{\Delta_i = 1(T_i < C_i)}, where
#' \eqn{C_i} is the censoring time. This function computes doubly robust scores \eqn{\hat \Gamma_i}
#' for the quantity 1) or 2).
#'
#' Let \eqn{H_i = min(Y_i, \text{horizon})}. Given estimates \eqn{\hat S_T(t; X_i)} and \eqn{\hat S_C(t; X_i)}
#' of the survival process \eqn{P[T_i > t | X_i]} and censoring process \eqn{P[C_i > t | X_i]},
#' the scores for 1) are given by:
#'
#' \deqn{\hat \Gamma_i =}
#' \deqn{\hat S_T(\text{horizon}; X_i) + }
#' \deqn{\sum_{t=1}^{H_i - 1} \frac{1}{\hat S_C(t; X_i)}
#'  \left(\frac{\hat S_T(\text{horizon}; X_i)}{\hat S_T(t; X_i)} - \frac{\hat S_T(\text{horizon}; X_i)}{\hat S_T(t - 1; X_i)}\right) +}
#' \deqn{\frac{\max(\Delta_i, 1(Y_i > \text{horizon}))}{\hat S_C(H_i; X_i)}
#'  \left(1(Y_i > \text{horizon}) -  \frac{\hat S_T(\text{horizon}; X_i)}{\hat S_T(H_i - 1; X_i)}\right).}
#'
#' Under regularity conditions, \eqn{\frac{1}{n} \sum_{i=1}^{n} \hat \Gamma_i} is an efficient estimate of \eqn{P[T_i > \text{horizon}]}.
#' The scores for 2) follow a similar construction.
#'
#' @param Y The event time.
#' @param D The event type \eqn{\Delta_i} (0: censored, 1: observed).
#' @param horizon A scalar time point which defines the estimand.
#' @param S.hat A matrix with estimates of the survival curves \eqn{P[T_i > t | X_i]}.
#' @param C.hat A matrix with estimates of the censoring curves \eqn{P[C_i > t | X_i]}.
#' @param Y.grid A vector with time grid \eqn{t = t_1, t_2, ..., t_M} S.hat and C.hat are estimated on.
#'
#' @return A list with scores along with an estimate of censoring probabilities at horizon `C.hat.H` that
#'  can be used to assess of positivity is an issue (i.e., if these probabilities are very low).
#'
#' @references Cui, Yifan, Michael R. Kosorok, Erik Sverdrup, Stefan Wager, and Ruoqing Zhu.
#'  "Estimating Heterogeneous Treatment Effects with Right-Censored Data via Causal Survival Forests".
#'  Journal of the Royal Statistical Society: Series B, 85(2), 2023.
#' @references Rubin, Daniel, and Mark van der Laan.
#'  "A doubly robust censoring unbiased transformation".
#'  The International Journal of Biostatistics, 3(1), 2007.
#' @references Wager, Stefan.
#'  "Causal Inference: A Statistical Learning Approach", 2024.
#'  (In chapter 14: doubly robust constructions for dynamic policies).
#'
#' @export# Compute scores for P[T > horizon] using the general recipe from the CSF paper.
get_survival_scores <- function(Y,
                                D,
                                horizon,
                                S.hat,
                                C.hat,
                                time.grid
                                ) {
  if (length(Y) != length(D) || length(Y) != NROW(S.hat)) {
    stop("Incompatible input sizes.")
  }
  if (!identical(dim(S.hat), dim(C.hat))) {
    stop("S.hat, C.hat should have same dimension.")
  }
  if (is.unsorted(time.grid, strictly = TRUE) || length(time.grid) != NCOL(S.hat)) {
    stop("time.grid should be a vector with length equal to ncol(S/C.hat).")
  }
  # TODO len(uniq(time.grid)), horizon > min time grid, div by zero
  # need min(Y) > min(time.grid) too

  # For every sample, get the time grid index corresponding to min(Yi, horizon)
  grid.index <- findInterval(pmin(Y, horizon), time.grid)
  # Get P[Ci > min(Yi, h) | Xi]
  C.Y.hat <- C.hat[cbind(seq_along(grid.index), grid.index)]

  # Construct the matrix Q(s, Xi) = P[Ti > horizon | Ti > s] = P[Ti > horizon | Xi] / P[Ti > s | Xi]
  horizon.grid.index <- findInterval(horizon, time.grid)
  
  
  #-copied from https://github.com/grf-labs/grf/blob/e2a2040690c3e461e793f98bce1c6f7163a8af7b/r-package/grf/R/causal_survival_forest.R#L514
  #to get Q.hat
  Y.diff <- diff(c(0, Y.grid))
  Q.hat <- matrix(NA, nrow(S.hat), ncol(S.hat))
  dot.products <- sweep(S.hat[, 1:(ncol(S.hat) - 1)], 2, Y.diff[2:ncol(S.hat)], "*")
  Q.hat[, 1] <- rowSums(dot.products)
  for (i in 2:(ncol(Q.hat) - 1)) {
    Q.hat[, i] <- Q.hat[, i - 1] - dot.products[, i - 1]
  }
  Q.hat <- Q.hat / S.hat
  Q.hat <- sweep(Q.hat, 2, Y.grid, "+") # Add back t
  Q.hat[, ncol(Q.hat)] <- max(Y.grid)
  #----------------------------------------------------------------------------------------------------------------------------------------
  
  Q.hat <- sweep(1 / S.hat, 1, S.hat[, horizon.grid.index], "*")
  # Get P[Ti > horizon | Ti > min(Yi, h)]
  Q.Y.hat <- Q.hat[cbind(seq_along(grid.index), grid.index)]

  log.surv.C <- -log(cbind(1, C.hat))
  dlambda.C.hat <- log.surv.C[, 2:(ncol(C.hat) + 1)] - log.surv.C[, 1:ncol(C.hat)]
  
  
  integrand <-  Q.hat * dlambda.C.hat / C.hat

  integral.correction <- rep(0, length(Y))
  for (sample in seq_along(Y)) {
    Yi.index <- grid.index[sample]
    integral.correction[sample] <- sum(integrand[sample, seq_len(Yi.index)])
  }

 
  
  
  Gamma.hat = ( D*pmin(Y, horizon) + (1 - D) * Q.Y.hat) / C.Y.hat - integral.correction
  
  list(Gamma.dr = Gamma.hat,
       C.hat.H = C.Y.hat)
}


# Try it out:
if (FALSE) {
  library(grf)
  n = 1000
  p = 5
  X = matrix(runif(n * p), n, p)
  failure.time = rexp(n, X[, 1] * 1 / 10)
  censor.time = rexp(n, X[, 2] * 1 / 30)
  Y = round(pmin(failure.time, censor.time), 1)
  D = as.integer(failure.time <= censor.time)
  
  Y.grid = sort(unique(Y)); time.grid = Y.grid
  sf.survival = survival_forest(X, Y, D,
                                failure.times = Y.grid,
                                num.trees = 250,
                                prediction.type = "Nelson-Aalen")
  S.hat = predict(sf.survival)$predictions
  sf.censor = survival_forest(X, Y, 1 - D,
                              failure.times = Y.grid,
                              num.trees = 250,
                              prediction.type = "Nelson-Aalen")
  C.hat = predict(sf.censor)$predictions
  
  horizon = 100
  dr = get_survival_scores(Y, D, horizon, S.hat, C.hat, Y.grid)
  estimate = mean(dr$Gamma.dr)
  print(estimate)
  
  # The "true" value:
  mean(rexp(1e6, X[, 1] * 1 / 10) > horizon)
}


