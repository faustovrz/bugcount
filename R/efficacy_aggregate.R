#'  Efficacy calculations as ratio of count means, geometric means, or medians
#' @param wf A whitefly count dataframe.
#' @param control Control genotype string.
#'        Genotype Denominator in efficacy: e = u1/u0.
#' @return A dataframe with efficacy calculated per clone
#'
#' @examples
#' # efficacy.df <- efficacy.aggregate(wf.per.plant ~ control = "COL1468")

efficacy_aggregate <- function(wf, control = "COL1468"){
  # This formula describes experimenl design
  form <- formula(nymphs ~ exp_cross + experiment + clone + group)
  wf_ag <- aggregate(form, data = wf,
                     FUN = function(x) c(n = length(x),
                                       mean  = mean(x, na.rm = TRUE), 
                                       geometric.mean = geometric.mean(x),
                                       median = median(x,na.rm = TRUE)))
  # fix aggregate column names
  wf_ag <- fix_ag_colnames(wf_ag)

  wf_x <- wf_ag[wf_ag$group == "infestation_check" & wf_ag$clone == control,
                !(colnames(wf_ag) %in% c("clone", "group","exp.cross"))]
  colnames(wf_x) <- gsub("nymphs.","infestation_",
                         colnames(wf_x), fixed = TRUE)
  colnames(wf_x) <- gsub("nymphs.","infestation_",
                         colnames(wf_x), fixed = TRUE)
  wf_y <- wf_ag[wf_ag$group != "infestation_check",
                colnames(wf_ag) != "group"]
  colnames(wf_y) <- gsub("nymphs.","clone_",
                         colnames(wf_y), fixed = TRUE) 
  wf_merge <- merge(wf_x, wf_y,
              by = "experiment",
              all.y = TRUE)

  wf_y <- wf_ag[wf_ag$group != "infestation_check" & wf_ag$clone == control,
                !(colnames(wf_ag) %in% c("clone", "group","exp.cross"))]
  colnames(wf_y) <- gsub("nymphs.","control_",
                         colnames(wf_y),fixed=TRUE)
  
  wf_merge <- merge( wf_merge, wf_y,
                    by="experiment",
                    all.y = TRUE)
  wf_merge <- within(wf_merge,{
    mean.eff <- 1 - clone_mean / control_mean
    geometric_mean_eff <- 1 - clone_geometric_mean / control_geometric_mean
    median <- 1 - clone_median / control_median
  })
  wf_merge
}
