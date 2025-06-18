#' Estimate normalized model weights based on the variance of the predictions
#'
#' This function estimates normalized model weights by evaluating the variance of each model's predictions.
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x A data frame containing the daily milking records for the individual
#' @return A vector containing the normalized weights for each model
#' @keywords internal
VarWeight<-function(converged_models, x){

  predictions <- lapply(converged_models, predict, newdata = x)

  variances <- sapply(predictions, var)
  weights_var <- 1 / variances
  normalized_weights <- weights_var / sum(weights_var)

  return(normalized_weights)

}
