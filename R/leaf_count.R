#'  Counts per leaf
#' 
#' @param wf A dataframe of whitefly count counts per picture.
#' @return Dataframe containing whitefly counts per leaf
#' @examples
#' # per.leaf <- leaf.count(my.wf)
#' @export

leaf_count <- function(wf){
 aggregate(as.formula(paste("nymphs ~ experiment +",  leaf_criteria(wf))), 
                       data = wf, FUN = sum, na.rm = TRUE)
}

