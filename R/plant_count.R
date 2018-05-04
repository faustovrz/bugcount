#'  Counts per plant
#' 
#' @param wf A dataframe of whitefly counts per picture
#' @return Dataframe containing whitefly counts per plant
#' @examples
#' # per_plant <- plant.count(my_wf)
#' @export

plant_count <- function(wf){
  wf_nleaf <- aggregate(as.formula(paste("leaf ~ ", 
                                         plant_criteria(wf))),
                        data = leaf_count(wf), FUN = length)

  wf_plant <- aggregate(as.formula(paste("nymphs ~ experiment +",
                                         plant_criteria(wf))),
                        data = leaf_count(wf), FUN = sum)
  wf_plant$nleaf <- wf_nleaf$leaf
  wf_plant[wf_plant$nleaf > 1,]
}

