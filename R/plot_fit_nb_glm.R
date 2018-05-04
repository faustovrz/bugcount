#'  Plot negative binomial model fit
#' @param wf A whitefly count dataframe.
#' @param np.fitted a \code{fit.nb.glm()} list.
#' @param type \code{"CI"} Confidence interval dotplot
#'             \code{"density"} probability density joyplot
#' @return Nothing
#'
#' @examples
#' # plot_fit_nb_glm(wf_per_leaf,
#' #                 fit_nb_glm(nymphs ~ year + leaf, wf.per.leaf),
#' #                 type = "density")
#' @export

plot_fit_nb_glm <- function(wf, nb_fitted, type = "violin", xmax =15000, cld = TRUE){
  group <- nb_fitted$group
  n_grp <- length(levels(group))
  grp_col <- gg_color_hue(n_grp)[group]
  response <- nb_fitted$response
  predictor <- nb_fitted$predictor
  order <- nb_fitted$order
  wf[,predictor] <- factor(wf[, predictor], levels = order)
    if (type == "CI") {
    plot(nb_fitted$posthoc , col = grp_col, 
         lwd = 25, xlab = "Nymphs", 
         ylab = as.character(group),
         main = "Arithmetic Mean Confidence Interval \n from GLM")
      
    } else if (type == "joy") {
    
    # change factor for joyplot to appear in the same order as posthoc plot

    ggplot2::ggplot(wf, 
      ggplot2::aes_string(
        x = response,
        y = predictor,
        fill = predictor)) +
      ggplot2::xlim(0, xmax) +
      ggplot2::ggtitle("Density") + ggjoy::geom_joy() + 
      ggplot2::scale_fill_manual(values = grp_col) +
      ggplot2::theme(legend.position = "none")
    } else if (type == "violin") {
      if (cld == TRUE){
        group_text <- group
      } else{
        group_text <- ""
      }
    ggplot2::ggplot(data = wf) + 
      ggplot2::ylim(0,xmax) +
      ggplot2::geom_violin(
        ggplot2::aes_string(y = response, x = predictor),
        scale = "width") + 
      ggplot2::geom_text(data = nb_fitted$posthoc,
                         x = 1:length(order),
                         y = xmax,
                         label = group_text) +
      ggplot2::geom_point(data = nb_fitted$posthoc, 
                          ggplot2::aes_string(x = predictor, y = "response"), 
                          shape = 15, size = 1) +
      ggplot2::geom_errorbar(
        data = nb_fitted$posthoc,
        ggplot2::aes_string( x = predictor, ymin = "asymp.LCL", ymax = "asymp.UCL"),
        width =  0.2, 
        size  =  0.7, col = "black") +
      ggpubr::theme_pubr() +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::coord_flip()
  }# TODO: else{stop("")}
}


