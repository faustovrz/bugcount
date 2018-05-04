#'  Remove levels form factors of subsetted dataframes
#' 
#' @param df A dataframe
#' @return Dataframe with excess levels removed from factors
#' @examples
#' # wf2000 <- remove.levels(plant.count[plant.count$year = 2000, ]
#' @export

remove_levels <- function(df){
  df[] <- lapply(df, function(x) if(is.factor(x)) factor(x) else x)
  df
} 
