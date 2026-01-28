#' Plot the actual daily milk daily production and the predicted values highlighting the detected milk loss events
#' @param data A data frame containing the observed and predicted daily milking records
#' @param ID The ID of the individual that will have the daily milking records plotted
#' @param res.milkloss The object with the output of milkloss_detect function
#' @param MY_col The name of the column containing the observed milk yield
#' @param MY_pred The name of the column containing the predicted milk yield
#' @param col The colors of the actual, predicted values, and milk loss events. In this order
#' @param id_col The name of the column containing the individual IDs
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2  scale_color_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme_bw
#' @return A plot with the actual daily milk daily production and the predicted values highlighting the detected milk loss events
#' @export

PlotMilkLoss <- function(data, ID, res.milkloss,
                         MY_col, MY_pred,
                         col = c("red", "blue", "darkgreen"),
                         id_col = "ID") {

  # Filter records for the chosen ID
  tmp <- data[which(data[, id_col] == ID), ]

  episodes_tbl <- res.milkloss$episodes[which(res.milkloss$episodes[, id_col] == ID), ] %>%
    tibble::as_tibble() %>%
    dplyr::mutate("{id_col}" := as.character(.data[[id_col]]))

  milk_tbl <- tmp %>%
    tibble::as_tibble() %>%
    dplyr::mutate("{id_col}" := as.character(.data[[id_col]]))

  shading_df <- milk_tbl %>%
    dplyr::inner_join(
      episodes_tbl,
      by = dplyr::join_by(
        !!rlang::sym(id_col),
        dplyr::between(DIM, start_DIM, end_DIM)
      ),
      relationship = "many-to-one"
    ) %>%
    dplyr::filter(.data[[MY_col]] <= .data[[MY_pred]]) %>%   # use dynamic cols here too
    dplyr::arrange(.data[[id_col]], episode_index, DIM)

  ggplot(tmp, aes(x = DIM)) +
    geom_ribbon(
      data = shading_df,
      aes(
        x = DIM,
        ymin = .data[[MY_col]],
        ymax = .data[[MY_pred]],
        group = episode_index
      ),
      fill = col[3], alpha = 0.35, inherit.aes = FALSE
    ) +
    geom_line(aes(y = .data[[MY_col]],  color = "MY_col"),  linewidth = 0.6) +
    geom_point(aes(y = .data[[MY_col]],  color = "MY_col"),  size = 1.6) +
    geom_point(aes(y = .data[[MY_pred]], color = "MY_pred"), size = 1.6, shape = 16) +
    scale_color_manual(
      NULL,
      values = c(
        MY_col  = col[1],
        MY_pred = col[2]
      )
    ) +
    labs(x = "DIM", y = "Milk yield") +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      axis.title  = element_text(size = 15)
    ) +
    scale_x_continuous(
      breaks = seq(0, max(tmp[, "DIM"], na.rm = TRUE), 5),
      labels = seq(0, max(tmp[, "DIM"], na.rm = TRUE), 5)
    )
}
