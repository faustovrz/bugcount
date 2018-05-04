#'  Plot infestation X clone from GLM
#' @param cld \code{multcomp} common letter display \code{cld} object from
#'        posthoc comparison in \code{plot_fit.nb.glm()}
#' @return Nothing
#'
#' @examples
#' # plot_gg.infestationXclone(nb.glm.fit$posthoc)
#' @export

plot_gg_infestationXclone <- function(cld){
  title <- "Interaction: Infestation Regime x Clone"
  subtitle <- "(GLM Arithmetic Means + 95% Confidence Interval)"

  ggplot2::ggplot(cld,
    ggplot2::aes(
      x     = clone,
      y     = response,
      color = infestation,
      label = .group),
    log10 = "y") +
    ggplot2::geom_point(shape  = 15, size   = 1) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = clone,
        y = asymp.LCL,
        label = .group),
      nudge_y = -0.05) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin  =  asymp.LCL,
        ymax  =  asymp.UCL),
      width =  0.2,
      size  =  0.7) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "bottom" ) +
    ggplot2::scale_y_continuous(trans  = scales::log10_trans(),
                       limits = c(10,50000),
                       breaks = c(10,100,100,1000,10000)) +
    ggplot2::ggtitle(paste(title,subtitle, sep = "\n"))
}
