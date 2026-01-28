#' Estimate normalized model`s weights based on a  Expectation–Maximization (EM) algorithm using a gamma distribution
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x  data frame containing the daily milking records for the individual
#' @param trait The name of the column containing daily milking records
#' @return A vector containing the model normalized weight based on the Bayeasian Model Average (BMA) appaoch proposed by Duan et al. (2006) with an adaptation for to use a gamma distribution
#' @keywords internal
BMAweight_gamma <- function(converged_models, x, trait) {

  # Step 1: Generate predictions for each model
  predictions <- lapply(converged_models, predict, newdata = x)
  predictions <- as.data.frame(predictions)

  k <- ncol(predictions)
  t <- nrow(predictions)

  observed <- x[, trait]

  weights <- rep(1 / k, k)
  variances <- colMeans((observed - predictions)^2)

  max_iter <- 100
  tol <- 1e-6
  log_likelihood <- -Inf

  for (iter in 1:max_iter) {
    z <- matrix(0, nrow = t, ncol = k)
    for (c in 1:k) {
      pred_mean <- predictions[, c]
      pred_var <- variances[c]

      shape <- (pred_mean^2) / pred_var
      scale <- pred_var / pred_mean

      shape[shape <= 0] <- 1e-6
      scale[scale <= 0] <- 1e-6

      g.dense <- dgamma(observed, shape = shape, scale = scale)

      z[, c] <- weights[c] * g.dense
    }

    z <- z / rowSums(z)

    z[is.na(z)]<-1e-6

    weights <- colMeans(z)

    alpha<-4

    weights<-weights^alpha

    weights <- weights / sum(weights)

    variances <- colSums(z * (observed - predictions)^2) / colSums(z)

    new_log_likelihood <- sum(log(rowSums(z)))

    if (abs(new_log_likelihood - log_likelihood) < tol) {
      break
    }
    log_likelihood <- new_log_likelihood
  }

  return(weights)
}
