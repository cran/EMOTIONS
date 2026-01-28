#' Estimate the Akaike information criterion (AIC), Bayeasian information criterion (BIC), root mean square percentage error (RMSPE) and mean squared error (MAE) for the fitted models
#'
#' @param converged_models A list containing the fitted models for the individual
#' @param x  data frame containing the daily milking records for the individual
#' @param trait The name of the column containing daily milking records
#' @importFrom quantreg AIC.rq
#' @return A data frame containing the AIC, BIC, RMSPE and MAE for each fitted model for the individual
#' @keywords internal
GetLacModelsMetrics<-function(converged_models, x, trait){

  metrics_list <- lapply(names(converged_models), function(model_name) {
    model <- converged_models[[model_name]]

    predictions <- predict(model, newdata = x)

    aic_value <- if (inherits(model, "rq")) {
      AIC.rq(model)
    } else if (inherits(model, c("lm", "nls"))) {
      AIC(model)
    } else {
      NA
    }

    bic_value <- if (inherits(model, "rq")) {
      AIC.rq(model, k=log(length(predictions)))
    } else if (inherits(model, c("lm", "nls"))) {
      AIC(model, k=log(length(predictions)))
    } else {
      NA
    }

    rmspe_value <- sqrt(mean((x[,trait] - predictions)^2))

    mae_value <- mean(abs(x[,trait] - predictions))

    data.frame(Model = model_name, AIC = aic_value, BIC=bic_value, RMSPE = rmspe_value, MAE = mae_value)
  })

  metrics_df <- do.call(rbind, metrics_list)

  return(metrics_df)

}
