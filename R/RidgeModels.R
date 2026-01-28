#' The function RidgeModels allows the visualization of the distribution of model's ranks across individuals using ridge density plots
#' @param LacCurveFit The object obtained from the LacCurveFit function
#' @param metric The name of the metric to be use to plot the model´s ranks
#' @return A ridge density plots for the models included in the ensemble
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_segment labs theme theme_minimal element_text
#' @importFrom dplyr arrange desc
#' @importFrom rlang sym
#' @importFrom ggridges  geom_density_ridges theme_ridges
#' @export
RidgeModels<-function(LacCurveFit,metric="AIC_rank"){

  data <- bind_rows(LacCurveFit$models_weight, .id = "ID")

  tab.count<-as.data.frame(table(data[,"Model"],data[,metric]))

  tab.count <- tab.count %>%
    arrange(Var2, desc(Freq)) %>% as.data.frame()

  levels.order<-tab.count[!duplicated(tab.count[,1]),1]

  data$Model<-factor(data$Model, levels = rev(levels.order))

  ggplot2::ggplot(data, aes(x = !!sym(metric), y = Model, fill = Model)) +
    geom_density_ridges(alpha = 0.7, scale = 1,bandwidth = 0.5) +
    labs(x = "Rank", y = "Model", title = "Distribution of Model Ranks Across Individuals") +
    theme_ridges() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), axis.title = element_text(size=15))

}
