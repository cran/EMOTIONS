#' Estimate normalized model weights based on the cosine similarity of each model's predictions
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x A data frame containing the daily milking records for the individual
#' @return A vector containing the normalized model weights
#' @keywords internal
CosSquaredWeight<-function(converged_models, x){

  predictions <- lapply(converged_models, predict, newdata = x)

  pred_matrix <- do.call(cbind, predictions)

  avg_prediction <- rowMeans(pred_matrix)

  cosine_similarity <- apply(pred_matrix, 2, function(pred) {
    sum(pred * avg_prediction) / (sqrt(sum(pred^2)) * sqrt(sum(avg_prediction^2)))
  })

  alpha<-2

  cosine_squared_weights <- cosine_similarity^alpha

  normalized_weights <- cosine_squared_weights / sum(cosine_squared_weights)

  return(normalized_weights)

}
