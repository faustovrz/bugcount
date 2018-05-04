#'Plot of efficacy probability density
#' just a wrapper for \code{reshape2::dcast} from \code{reshape2}
#' @param ef efficacy dataframe from \code{efficacy.aggregate()} 
#' @return List of wide formated dataframes of  efficacy by mean,
#'         geometric mean, median by clone for each experiment
#'
#' @examples
#' # ef_by_measure <- plot_ef_density(efficacy.df)

plot_ef_density <- function(ef,type, group){
  n_group <- length(unique(ef[,group]))
  ggplot2::ggplot(ef,aes_string(x = type, y = group, fill = group)) +
    ggplot2::ggtitle("Efficacy Distribution") + 
    ggplot2::geom_joy() +  
    ggplot2::xlim(-2, max(ef[,type])) +
    ggplot2::scale_fill_cyclical(values = gg_color_hue(n_group)) +
    ggplot2::xlab("Efficacy")  
}
