#' Estimate normalized model`s weights based on the variance of the predictions
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x  data frame containing the daily milking records for the individual
#' @return A vector containing the model normalized weight
#' @keywords internal
VarWeight<-function(converged_models, x){

  predictions <- lapply(converged_models, predict, newdata = x)

  variances <- sapply(predictions, var)
  weights_var <- 1 / variances
  normalized_weights <- weights_var / sum(weights_var)

  return(normalized_weights)

}
