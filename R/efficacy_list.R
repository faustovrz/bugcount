#'  Vertical to horizontal data reformatting for efficacy calculations
#' just a wrapper for \code{reshape2::dcast} from \code{reshape2}
#' @param ef efficacy dataframe from \code{efficacy.aggregate()} 
#' @return List of wide formated dataframes of efficacy by mean,
#'         geometric mean median
#'
#' @examples
#' # ef.by.measure <- efficacy.list(efficacy.df)

efficacy_list <- function(ef){
list(mean = reshape2::dcast(ef,clone ~ experiment,value.var = "mean.eff"),
     geometric.mean = reshape2::dcast(ef,clone ~ experiment,value.var = "geometric.mean.eff"))
#   median = reshape2::dcast(wf.merge,clone ~ experiment,value.var = "median"))

}

