#' Construct doubly robust pseudo-outcomes for E[y(T_i)] and use those in place of Yi, i.e., Y_i^* = y(Ti)/PSi if Wi=1,  Y_i^* = -y(Ti)/(1-PSi) if Wi=0. 
#' Author: Erik Sverdrup
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
  D <- pmax(D, Y > horizon) #redefine "D" to be the "effective censoring indicator \Delta_i^h" like in the CSF paper
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


