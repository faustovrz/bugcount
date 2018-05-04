#'  Aggregate whitefly data by single criterion
#' @param wf A whitefly count dataframe.
#' @param x criterion column for aggregation.
#' @param y variable column for aggregation.
#' @param by_stat column for ordering levels for graphs.
#' @return Dataframe of aggregates by specified {col} 
#'         values and levels ordered by \code{by_stat}
#'
#' @examples
#' # wf_aggregate(my.wf,x="clone", y="nymphs", by_stat ="geometric.mean" )
#' @export

wf_aggregate <- function(wf,x = "clone", y = "nymphs", by_stat = "mean"){
  form <- formula(paste(y, x, sep = " ~ "))
  wf_ag <- aggregate(form, 
                      data = wf,
                      FUN = function(x) {
                        c(n = length(x),
                          mean  = mean(x, na.rm = TRUE), 
                          gmean  = geometric_mean(x),
                          median = median(x, na.rm = TRUE))
                      }
              )
  # fix aggregate column names
  wf_ag <- fix_ag_colnames(wf_ag)

  # sort by_stat
  wf_ag <-  wf_ag[order(wf_ag[,paste(y,by_stat,sep = ".")]),]
                        
  # sort x levels by_stat
  # wf.sort <- function(wf,by_stat){}
  wf_ag[,x] <- factor(wf_ag[,x], levels = wf_ag[,x])
  wf_ag
}
