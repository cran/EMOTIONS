#' A function to estimate resilience estimators (logarithm of variance, lag1 autocorrelation and skewness) based on daily milk production records
#' @param production_df The list containing the data frames with the daily production records (actual or predicted) obtained from the LacCurveFit function
#' @param dim_filter_range A vector containing the lower and upper limits to remove lactation records from the begin and end of the lactation, if needed. If it is not necessary to remove daily records, the first two values can be set as the minimum days in milk value and the last two as the maximum days in milk values
#' @param outlier_sd_threshold A threshold defining the maximum standard deviations to consider an individual resilience indicator value
#' @param weight The name of the column containing the selected ensemble prediction. The default is weight_AIC
#' @param trait The name of the column containing daily milking records
#' @param DIM The name of the column containing days in milk records
#' @param ID_col The name of the column containing the unique IDs of the individuals
#' @importFrom dplyr %>% bind_rows group_by summarize summarise filter select left_join
#' @importFrom parameters skewness
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect everything
#' @return A list containing the daily milk production values after filtering, the list of removed animals, and a data frame with the resilience indicators
#' @export
ResInd <- function(production_df,dim_filter_range = c(1, 7, 203, 210),
                                  outlier_sd_threshold = 4,
                                  weight="weight_AIC",trait,DIM,ID_col) {

  production_df<-do.call(rbind,production_df)

  production_df$CodGen <- production_df[,ID_col]

  production_df$DIM <- production_df[,DIM]

  production_df$dev <- production_df[,trait] - production_df[,weight]

  desviaciones_filtradas <- production_df %>%
    filter(!(DIM >= dim_filter_range[1] & DIM <= dim_filter_range[2]) &
             !(DIM >= dim_filter_range[3] & DIM <= dim_filter_range[4]))

  ri <- desviaciones_filtradas %>%
    group_by(CodGen) %>%
    summarize(
      log_varianza = log(var(dev, na.rm = TRUE)),
      autocorrelacion_lag1 = {
        dev_sin_na <- dev[!is.na(dev)]
        if (length(dev_sin_na) > 1) {
          acf_result <- acf(dev_sin_na, lag.max = 1, plot = FALSE)
          acf_result$acf[2]
        } else {
          NA
        }
      },
      skewness = skewness(dev, na.rm = TRUE)$Skewness
    )


  ri_summary <- ri %>%
    summarise(
      mean_log_varianza = mean(log_varianza, na.rm = TRUE),
      sd_log_varianza = sd(log_varianza, na.rm = TRUE),
      mean_autcorrelacion_lag1 = mean(autocorrelacion_lag1, na.rm = TRUE),
      sd_autocorrelacion_lag1 = sd(autocorrelacion_lag1, na.rm = TRUE),
      mean_skewness = mean(skewness, na.rm = TRUE),
      sd_skewness = sd(skewness, na.rm = TRUE)
    )

  ri_filtered <- ri %>%
    filter(
      log_varianza >= (ri_summary$mean_log_varianza - outlier_sd_threshold * ri_summary$sd_log_varianza) &
        log_varianza <= (ri_summary$mean_log_varianza + outlier_sd_threshold * ri_summary$sd_log_varianza) &
        autocorrelacion_lag1 >= (ri_summary$mean_autcorrelacion_lag1 - outlier_sd_threshold * ri_summary$sd_autocorrelacion_lag1) &
        autocorrelacion_lag1 <= (ri_summary$mean_autcorrelacion_lag1 + outlier_sd_threshold * ri_summary$sd_autocorrelacion_lag1) &
        skewness >= (ri_summary$mean_skewness - outlier_sd_threshold * ri_summary$sd_skewness) &
        skewness <= (ri_summary$mean_skewness + outlier_sd_threshold * ri_summary$sd_skewness)
    )

  removed_samples <- setdiff(ri$CodGen, ri_filtered$CodGen)


  media_prod_por_codgen <- production_df %>%
    group_by(CodGen) %>%
    summarise(mean_Prod = mean(.data[[trait]], na.rm = TRUE))

  ri_filtered <- ri_filtered %>%
    left_join(media_prod_por_codgen, by = "CodGen")


  dev_list <- production_df %>%
    select(CodGen, DIM, dev)

  rownames(dev_list) <- NULL


  ri_stats <- ri_filtered %>%
    summarise(
      log_varianza_mean = mean(log_varianza, na.rm = TRUE),
      log_varianza_sd = sd(log_varianza, na.rm = TRUE),
      log_varianza_min = min(log_varianza, na.rm = TRUE),
      log_varianza_max = max(log_varianza, na.rm = TRUE),
      autocorrelacion_lag1_mean = mean(autocorrelacion_lag1, na.rm = TRUE),
      autocorrelacion_lag1_sd = sd(autocorrelacion_lag1, na.rm = TRUE),
      autocorrelacion_lag1_min = min(autocorrelacion_lag1, na.rm = TRUE),
      autocorrelacion_lag1_max = max(autocorrelacion_lag1, na.rm = TRUE),
      skewness_mean = mean(skewness, na.rm = TRUE),
      skewness_sd = sd(skewness, na.rm = TRUE),
      skewness_min = min(skewness, na.rm = TRUE),
      skewness_max = max(skewness, na.rm = TRUE),
      mean_Prod_mean = mean(mean_Prod, na.rm = TRUE),
      mean_Prod_sd = sd(mean_Prod, na.rm = TRUE),
      mean_Prod_min = min(mean_Prod, na.rm = TRUE),
      mean_Prod_max = max(mean_Prod, na.rm = TRUE)
    ) %>%

    pivot_longer(
      cols = everything(),
      names_to = c("indicator", ".value"),
      names_pattern = "(.*)_(mean|sd|min|max)"
    )


  return(list(
    ri_filtered = ri_filtered,
    dev_list = dev_list,
    removed_samples = removed_samples,
    ri_stats = ri_stats
  ))
}
