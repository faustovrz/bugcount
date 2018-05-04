# Compatibility of graphic device with windows OS.
if (.Platform$OS.type == "windows") {
  quartz <- function() windows()
}

# Basic Function Initialization ################################################
# Make vectors of descriptors from the expected column names of the whitefly
# data frame

exp_descriptor <- function(wf){
  if (any(!grepl("experiment",colnames(wf)))) {
    descriptor <- grep("exp_",colnames(wf), value = TRUE)
    c("experiment",exp_descriptor)
    } else {
      descriptor <- grep("exp_",colnames(wf), value = TRUE)
    }
  descriptor
}

exp_criteria <- function(wf){
  paste(exp.descriptor(wf), collapse = " + ")
  
} 

plant_descriptor <-  function(wf){
  c("group","propagation","substrate","clone", "rep","plant")
}

plant_criteria <- function(wf){
  paste(c(exp_descriptor(wf),plant_descriptor(wf)), collapse = " + ")
}


leaf_descriptor <- function(wf){
  c(plant_descriptor(wf),"leaf")
}

leaf_criteria <- function(wf){
  paste(c(exp_descriptor(wf), leaf_descriptor(wf)), collapse = " + ")
}

# fix aggregate column names
fix_ag_colnames <- function(ag){
  do.call(cbind.data.frame,ag)
}


#' Histogram panel for use un correlation matrix

panel_hist <- function(x, ...){ 
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

#'  Correlation panel for use un correlation matrix

panel_cor <- function(x, y, digits=2, prefix="")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  cex <- 0.4/strwidth(txt)
  
  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  
  text(0.5, 0.5, txt, cex = cex, col = "darkgrey")
  text(.65, .6, Signif, cex = cex , col = "darkgrey", pos = 4)
}


# Correlation panel for use un correlation matrix

panel_points <- function(x,y){
  points(x,y)
  abline(stats::lm(y ~ x), col = "red")
}



# Variance decomposition, Heritability, BLUP ###################################

clone_exp <- function(clone,wf){
  wf_exp <- unique(wf[,c("clone","experiment")])
  sort(as.vector(wf.exp$experiment[wf.exp$clone == clone]))
}

cross_exp <- function(cross,wf){
  wf_exp <- unique(wf[,c("exp_cross","experiment")])
  sort(as.vector(wf_exp$experiment[wf_exp$exp_cross == cross]))
}


select_complete_clones <- function(wf) {
  clones <- unique(as.character(wf$clone))
  experiments <- unique(as.character(wf$experiment))
  complete <- c()
  for (cross in unique(as.character(wf$exp_cross)) ) {
    complete_exp <- grep(cross, experiments, value = TRUE)
    
    if (length(complete.exp) > 0) {
      check_complete <- unlist(
        lapply(clones, 
               FUN = function(x){
                 all( complete_exp %in% clone_exp(x, wf))
               }))
      complete <- unique(c(complete, clones[check_complete]))
    }
  }
  sort(complete)
}
