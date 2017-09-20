library(car)
library(MASS)
library(ggjoy)
library(scales) 
library(reshape2)
library(psych)
library(lsmeans)
library(multcomp)

# Compatibility of graphic device with windows OS.
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

# Basic Function Initialization ################################################
# Make vectors of descriptors from the expected column names of the whitefly
# data frame

exp.descriptor <- function(wf){
  if (any(!grepl("experiment",colnames(wf)))){
    descriptor <- grep("exp.",colnames(wf), value=TRUE)
    c("experiment",exp.descriptor)
    } else {
      descriptor <- grep("exp.",colnames(wf), value=TRUE)
    }
  descriptor
}

exp.criteria <- function(wf){
  paste(exp.descriptor(wf), collapse= " + ")
  
} 

plant.descriptor <-  function(wf){
  c("group","propagation","substrate","clone", "rep","plant")
}

plant.criteria <- function(wf){
  paste(c(exp.descriptor(wf),plant.descriptor(wf)), collapse= " + ")
}


leaf.descriptor <- function(wf){
  c(plant.descriptor(wf),"leaf")
}

leaf.criteria <-function(wf){
  paste(c(exp.descriptor(wf),leaf.descriptor(wf)), collapse= " + ")
}

# Data Input     ##############################################################

#' Read whitefly count table
#' 
#' @param file A filename for the whitefly counts per picture.
#' @return Dataframe containing whitefly counts
#' @details 
#' Whitefly counts per picture. Each picture must have a unique combination of
#' experiment + replicate + plant + pic
#' experiment description consists of:
#' exp.year + exp.cross + exp.propagation + exp.substrate
#' @examples
#' my.wf <- read.wf("whitefly_counts.tab")

read.wf <- function(file,...){
  wf<- read.table(file="WF_consolidated.tab", sep="\t",header=TRUE)
  wf<-wf[!(is.na(wf$nymphs)),]
  # Add experiment as factor = exp.propagation + exp.substrate + exp.year
  wf$experiment <- apply( wf[, exp.descriptor(wf) ] ,
                          1 , paste , collapse = "_" )
  wf
  # TODO: order columns with espriment at the beginning. new function?
}

#' Aggregate whitefly counts per leaf
#' 
#' @param wf A dataframe of whitefly count counts per picture.
#' @return Dataframe containing whitefly counts per leaf
#' @details 
#' @examples
#' per.leaf <- leaf.count(my.wf)

leaf.count <- function(wf){
 aggregate(as.formula(paste ("nymphs ~ ",  leaf.criteria(wf)) ), 
                       data = wf, FUN=sum, na.rm=TRUE)
} 

#' Aggregate whitefly counts per plant
#' 
#' @param wf A dataframe of whitefly counts per picture
#' @return Dataframe containing whitefly counts per plant
#' @examples
#' per.plant <- plant.count(my.wf)

plant.count <- function(wf){
  
  wf.nleaf <- aggregate(as.formula(paste ("leaf ~ ",  plant.criteria(wf))),
                        data=leaf.count(wf), FUN=length)
  
  wf.plant <- aggregate(as.formula(paste ("nymphs ~ ",plant.criteria(wf))),
                        data=wf.leaf, FUN=sum)
  wf.plant$nleaf <- wf.nleaf$leaf
  wf.plant[wf.plant$nleaf>1,]
}


#' Geometric mean of a vector ignoring NAs.
#' 
#' @param x A vector.
#' @return Geometric mean of  \code{x} ignoring NAs.
#'
#' @examples
#' gm.mean(1:1000)

geometric.mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# fix aggregate column names
fix.ag.colnames <- function (ag){
  do.call(cbind.data.frame,ag)
}


#' Remove levels form factors of subsetted dataframes
#' 
#' @param df A dataframe
#' @return Dataframe with excess levels removed from factors
#' @examples
#' wf.2000 <- remove.levels(plant.count[plant.count$year = 2000, ]

remove.levels <-function(df){
  df[] <- lapply(df, function(x) if(is.factor(x)) factor(x) else x)
  df
} 

#' Aggregate whitefly data by single criterion
#' @param wf A whitefly count dataframe.
#' @param x criterion column for aggregation.
#' @param y variable column for aggregation.
#' @param by.stat column for ordering levels for graphs.
#' @return Dataframe of aggregates by specified {col} 
#'         values and levels ordered by \code{by.stat}
#'
#' @examples
#' wf.aggregate(my.wf,x="clone", y="nymphs", by.stat ="geometric.mean" )

wf.aggregate <- function (wf,x="clone", y="nymphs", by.stat ="mean"){
  form <- formula(paste(y,x, sep = " ~ "))
  wf.ag <- aggregate(form, data=wf,
                      FUN=function(x) c(n = length(x),
                                        mean  = mean(x, na.rm=TRUE), 
                                        gmean = geometric.mean(x),
                                        median = median(x,na.rm=TRUE)))
  # fix aggregate column names
  wf.ag <- fix.ag.colnames(wf.ag)

  # sort by.stat
  wf.ag <-  wf.ag[order(wf.ag[,paste(y,by.stat,sep = ".")]),]
                        
  # sort x levels by.stat
  # wf.sort <- function(wf,by.stat){}
  wf.ag[,x] <- factor(wf.ag[,x], levels=wf.ag[,x])
  wf.ag
}


# Distribution fit #############################################################

#' Equivalent to `MASS::fitdistr(x, densfun = "gamma")`, where x are first
#' rescaled to the appropriate scale for a gamma distribution.
#' Useful for fitting the gamma distribution to 
#' data which, when multiplied by a constant, follows this distribution
#' @param x A vector.
#' @return MASS::fitdistr with scaling multiplier
#'
#' @examples
#' fitgamma(1:1000)

fitgamma <- function(x) {
  if (!requireNamespace("MASS")) stop("Requires MASS package.")
  fit <- glm(formula = x ~ 1, family = Gamma)
  out <- MASS::fitdistr(x * coef(fit), "gamma")
  out$scaling_multiplier <- unname(coef(fit))
  out
}

#' Fit different distributions to count data
#' @param wf A whitefly count dataframe.
#' @return List of distribution parameters for diferent fits.
#'         and histogram data for plotting.
#'
#' @examples
#' fit.list <- count.fit(plant.count)

count.fit <- function(wf){
  poisson <- fitdistr(wf$nymphs, "Poisson")
  nbinom <- fitdistr(wf$nymphs, "Negative Binomial")
  #gammad <- fitgamma(wf.plant$nymphs+1)
  
  nymph.hist <-  hist(wf$nymphs, breaks=50, plot=FALSE)
  mids <- nymph.hist$mids
  
  list( hist =nymph.hist,
        nb.d = dnbinom(mids,
                       mu = nbinom$estimate[[2]],
                       size = nbinom$estimate[[1]]),
        
        lnorm.d = dlnorm(mids,
                         meanlog = mean(log(wf$nymphs+1)),
                         sdlog = sd(log((wf$nymphs+1)))),
        
        norm.d = dnorm(mids,
                        mean = mean(wf$nymphs),
                        sd = sd(wf$nymphs)),
        
        pois.d = dpois(mids,
                        lambda =poisson$estimate[[1]])
  )
}

#' Plot Probability Density over histogram to appreciate fit
#' @param x A \code{count.fit()} list
#' @return Plot.
#'
#' @examples
#' plot.count.fit(count.fit(wf.per.leaf))


plot.count.fit <- function(x,...){
  plot(x$hist, xlab = "Nypmh Count",...)
  lines(x$hist$mids,max(x$hist$counts) * x$nb.d / max(x$nb.d),
         col = gg_color_hue(4)[1], lwd=2, lty = 1)
  lines(x$hist$mids,max(x$hist$counts) * x$lnorm.d / max(x$lnorm.d),
        col = gg_color_hue(4)[2], lwd=2, lty = 2)
  lines(x$hist$mids,max(x$hist$counts) * x$pois.d / max(x$pois.d),
        col = gg_color_hue(4)[3], lwd=2, lty = 3)
  lines(x$hist$mids,max(x$hist$counts) * x$norm.d / max(x$norm.d),
        col = gg_color_hue(4)[4], lwd=2, lty = 4)
  
  # logspline aproximation, not good.
  # lines(wf.hist$mids,max(wf.hist$counts)*logspline.d/max(logspline.d))
  
  legend("topright",legend=factor(c("Negative Binomial", "Lognormal",
                                    "Poisson","Normal"),
                                  levels = c("Negative Binomial","Lognormal",
                                             "Poisson","Normal")),
         title = "Model Aproximation",
         col= gg_color_hue(4), lwd = 2, lty=1:4, bty ="n")
  
  # Q-Q plots
  
  # Innecesary for now
  # lnorm means lognormal
  

  # qqp(wf.plant$nymphs, ylab = "Nymph Counts")
  

  # qqp(wf.plant$nymphs + 1, "lnorm", ylab = "Nymph Counts")

  # qqp(wf.plant$nymphs, "pois", poisson$estimate, ylab = "Nymph Counts")
  
  

  # qqp(wf.plant$nymphs, "nbinom", 
  #     size = nbinom$estimate[[1]], 
  #     mu = nbinom$estimate[[2]],
  #     ylab = "Nymph Counts")
  

  # qqp(wf.plant$nymphs+1, "gamma", 
  #     shape = gammad$estimate[[1]],
  #     rate = gammad$estimate[[2]])
  recordPlot()
}

# Negative Binomial GLM fit ####################################################

#' Fit MASS Negative binomial General Linear Model for whitefly counts
#' @param formula formula object specifiying the GLM
#' @param wf A whitefly count dataframe. 
#' @return List with model fit, lsmeans and posthoc groups
#'
#' @examples
#' fit.nb.glm(nymphs ~ year + leaf, wf.per.leaf)

fit.nb.glm <-function(formula,wf){
  # TODO: validate formula for one predictor only
  response <- as.character(terms(formula)[[2]])
  predictor <- as.character(terms(formula)[[3]])
  fit <- glm.nb(formula, data = wf)
  lsm.fit <- lsmeans(fit, predictor,
                     adjust = "tuckey",
                     data = wf)
  ls.posthoc <- cld(lsm.fit, type = "response")
  #ls.posthoc$experiment <- gsub(" ","",ls.posthoc$.group)
  list( formula = formula,
        response = response,
        predictor = predictor,
        glm = fit, lsm = lsm.fit,
        posthoc = ls.posthoc,
        order = ls.posthoc[,predictor],
        group = as.factor(ls.posthoc$.group)
        )
}

#' Plot negative binomial model fit
#' @param wf A whitefly count dataframe.
#' @param np.fitted a \code{fit.nb.glm()} list.
#' @param type \code{"CI"} Confidence interval dotplot
#'             \code{"density"} probability density joyplot
#' @return Nothing
#'
#' @examples
#' plot.fit.nb.glm(wf.per.leaf,
#'                 fit.nb.glm(nymphs ~ year + leaf, wf.per.leaf),
#'                 type = "density")

plot.fit.nb.glm <- function(wf, nb.fitted, type = "CI", xmax =15000){
  group <- nb.fitted$group
  n.grp <- length(levels(group))
  grp.col <- gg_color_hue(n.grp)[group]
  response <- nb.fitted$response
  predictor <- nb.fitted$predictor
  order <- nb.fitted$order
  
    if(type=="CI"){
    plot(nb.fitted$posthoc , col = grp.col, 
         lwd =25, xlab = "Nymphs", 
         ylab = as.character(group),
         main = "Arithmetic Mean Confidence Interval \n from GLM")
      
    } else if(type=="density"){
    
    # change factor for joyplot to appear in the same order as posthoc plot
    wf[,predictor] <- factor(wf[,predictor], levels = order)
    ggplot(wf, aes_string(x = response, y = predictor, fill = predictor)) +
       xlim(0, xmax) +
       ggtitle("Density") + geom_joy() + 
       scale_fill_manual(values = grp.col) +
       theme(legend.position="none")
    
  }# TODO: else{stop("")}
}


#' Plot infestation X clone interaction as boxplot
#' @param wf A whitefly count dataframe.
#' @return Nothing
#'
#' @examples
#' plot.infestationXclone(wf.per.leaf)

plot.infestationXclone <-function(wf){
  # TODO:  validate that infestation is a column of wf
  wf.aggregate(wf, by.stat = "gmean")
  clone.order <- wf.aggregate(wf, by.stat = "gmean")$clone
  n.clones <-length(clone.order)
  
  wf$clone <- factor(wf$clone, levels = clone.order)
  
  par(mar=c(6,4,4,2)) 
  boxplot(nymphs+1 ~  infestation + clone, data=wf, frame = FALSE,
          border = "white", log="y", xaxt="n", yaxt="n")
  
  rect(seq(0.5,2*n.clones, by=4),
       rep(1,n.clones),
       seq(2.5,2*n.clones+0.5, by=4),
       rep(100000,n.clones),
       col="grey90",lty=0)
  
  boxplot(nymphs+1 ~  infestation + clone, las=2, log = "y",
          main = "Interaction: Infestation Regime x Clone",
          ylab = "Nymphs per Unit",xaxt="n", xlab = "",
          col= c("white", "grey"),
          data = wf, add=TRUE)
  
  axis(1, at=2*(1:n.clones)-0.5,labels=clone.order, las=2)
  
  mtext("Clone (sorted by geometric mean of nymphs)", side=1, line=5)
  
  legend(2*n.clones - 3,6, legend=c("High", "Low"),
         title = "Infestation",
         bg = "white",
         fill=c("grey", "white"))
}

#' Plot infestation X clone from GLM
#' @param cld \code{multcomp} common letter display \code{cld} object from
#'        posthoc comparison in \code{plot.fit.nb.glm()}
#' @return Nothing
#'
#' @examples
#' plot.gg.infestationXclone(nb.glm.fit$posthoc)

plot.gg.infestationXclone<-function(cld){
  title <- "Interaction: Infestation Regime x Clone"
  subtitle <- "(GLM Arithmetic Means + 95% Confidence Interval)"

  ggplot(cld,
         aes(x     = clone,
             y     = response,
             color = infestation,
             label = .group),
         log10="y") +
    geom_point(shape  = 15,
               size   = 4) +
    geom_text(aes( x= clone,
                   y = asymp.UCL,
                   label=.group),
                   nudge_y = 0.1) +
    geom_errorbar(aes(ymin  =  asymp.LCL,
                      ymax  =  asymp.UCL),
                  width =  0.2,
                  size  =  0.7) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5), 
          legend.position =  c(0.9, 0.2) ) +
    scale_y_continuous(trans=log10_trans(),
                       limits = c(1,10000),
                       breaks =c(1,10,100,1000,10000)) +
    ggtitle(paste(title,subtitle, sep = "\n"))
}

# Efficacy calculations ########################################################

#' Vertical to horizontal data reformatting and aggregation for whitefly counts
#' just a wrapper for \code{dcast} from \code{reshape2}
#' @param formula formula object specifiying the \code{dcast} reformatting
#' @param wf A whitefly count dataframe.
#' @param fun function for aggregation deafult \code{geometric.mean}
#' @return A wide formated dataframe
#'
#' @examples
#' clone.means <- wf.wide(nymphs ~ clone,fun = mean)

wf.wide <- function(formula, wf, fun = geometric.mean){
  dcast(wf,formula,
        value.var = "nymphs",
        fun.aggregate = fun)
}

# plot.wide.cor<-function(wf.wide, fun = "Geometric Means"){
#   pairs.panels(log10(wf.wide[,-1]), smooth = FALSE, density=FALSE,
#                ellipses=FALSE, rug = FALSE, pch = 1, lm= TRUE, stars = TRUE,
#                main = paste(cross,fun, "Correlation"))
# }


#' Efficacy calculations as ratio of count means, geometric means, or medians
#' @param wf A whitefly count dataframe.
#' @param control Control genotype string.
#'        Genotype Denominator in efficacy: e = u1/u0.
#' @return A dataframe with efficacy calculated per clone
#'
#' @examples
#' efficacy.df <- efficacy.aggregate(wf.per.plant ~ control = "COL1468")

efficacy.aggregate<- function (wf, control = "COL1468"){
  wf <-wf.count
  # This formula describes experimenl design
  form <- formula(nymphs ~ exp.cross + experiment + clone + group)
  wf.ag <- aggregate(form, data=wf,
                     FUN=function(x) c(n = length(x),
                                       mean  = mean(x, na.rm=TRUE), 
                                       geometric.mean = geometric.mean(x),
                                       median = median(x,na.rm=TRUE)))
  # fix aggregate column names
  wf.ag <- fix.ag.colnames(wf.ag)

  wf.x <- wf.ag[wf.ag$group == "infestation_check" & wf.ag$clone == control,
                !(colnames(wf.ag) %in% c("clone", "group","exp.cross"))]
  colnames(wf.x) <- gsub("nymphs.","infestation_",
                         colnames(wf.x),fixed=TRUE)
  colnames(wf.x) <- gsub("nymphs.","infestation_",
                         colnames(wf.x),fixed=TRUE)
  wf.y <- wf.ag[wf.ag$group != "infestation_check",
                colnames(wf.ag) != "group"]
  colnames(wf.y) <- gsub("nymphs.","clone_",
                         colnames(wf.y),fixed=TRUE) 
  wf.merge <- merge(wf.x, wf.y,
              by="experiment",
              all.y = TRUE)

  wf.y <- wf.ag[wf.ag$group != "infestation_check" & wf.ag$clone == control,
                !(colnames(wf.ag) %in% c("clone", "group","exp.cross"))]
  colnames(wf.y) <- gsub("nymphs.","control_",
                         colnames(wf.y),fixed=TRUE)
  
  wf.merge <- merge(wf.merge, wf.y,
                    by="experiment",
                    all.y = TRUE)
  
  wf.merge$mean.eff <- 1 - wf.merge$clone_mean/wf.merge$control_mean
  wf.merge$geometric.mean.eff <- 1 - wf.merge$clone_geometric.mean/wf.merge$control_geometric.mean
  wf.merge$median <- 1 - wf.merge$clone_median/wf.merge$control_median
  wf.merge
}

#' Vertical to horizontal data reformatting for efficacy calculations
#' just a wrapper for \code{dcast} from \code{reshape2}
#' @param ef efficacy dataframe from \code{efficacy.aggregate()} 
#' @return List of wide formated dataframes of efficacy by mean,
#'         geometric mean median
#'
#' @examples
#' ef.by.measure <- efficacy.list(efficacy.df)

efficacy.list <- function(ef){
list(mean = dcast(ef,clone ~ experiment,value.var = "mean.eff"),
     geometric.mean = dcast(ef,clone ~ experiment,value.var = "geometric.mean.eff"))
#   median = dcast(wf.merge,clone ~ experiment,value.var = "median"))

}


#' Plot of efficacy probability density
#' just a wrapper for \code{dcast} from \code{reshape2}
#' @param ef efficacy dataframe from \code{efficacy.aggregate()} 
#' @return List of wide formated dataframes of  efficacy by mean,
#'         geometric mean, median by clone for each experiment
#'
#' @examples
#' ef.by.measure <- efficacy.list(efficacy.df)

plot.ef.density <-function(ef,type, group){
  n.group <- length(unique(ef[,group]))
  ggplot(ef, 
         aes_string(x = type,
                    y = group,
                    fill = group)) +
    ggtitle("Efficacy Distribution") + 
    geom_joy() +  
    xlim(-2, max(ef[,type])) +
    scale_fill_cyclical(values = gg_color_hue(n.group)) +
    xlab("Efficacy")  
}

#' Emulate ggplot color palette.
#' @param n number of colors to generate
#' @return vector of hexadecimal color strings
#'
#' @examples
#' gg_color_hue(length(levels(wf$clone)))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# Correlation Analysis, reproducibility ########################################

#' Correlation panel for use un correlation matrix

panel.hist <- function(x, ...)
{ usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

#' Correlation panel for use un correlation matrix

panel.points <- function(x,y){
  points(x,y)
  abline(stats::lm(y ~ x), col = "red")
}

#' Correlation panel for use un correlation matrix

panel.cor <- function(x, y, digits=2, prefix="")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  cex <- 0.4/strwidth(txt)
  
  test <- cor.test(x,y)
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))
  
  text(0.5, 0.5, txt, cex = cex, col = "darkgrey")
  text(.65, .6, Signif, cex = cex , col = "darkgrey", pos= 4)
}

#' Plot wide formated dataframe correlations per clone
#' in log scale
#' @param wf.wide wide formatted data by clone
#' @return \code{pairs} plot of correlation
#'
#' @examples
#' plot.wide.cor(efficacy.list(ef)$mean)

plot.wide.cor<-function(wf.wide, main = "Correlation", ...){
  pairs(log10(wf.wide[,-1]),
        lower.panel = panel.points,
        diag.panel = panel.hist,
        upper.panel = panel.cor, 
        main = main)
}

#' Plot wide formated dataframe correlations per clone
#' with point labels, asumes plotted points to be ~ 10
#' more points would be extremely unreadable
#' @param check.wide wide formatted data by clone for checks only
#' @param main plot title 
#' @return \code{pairs} plot of correlation
#'
#' @examples
#' plot.check.cor(eplot.check.cor)
#' 
plot.check.cor <-function(check.wide,
                          main = "Correlation",
                          ...){
  pairs(check.wide[,-1],
        lower.panel = function(x,y){
          text(x,y, labels=check.wide$clone, cex= 1)
          abline(stats::lm(y ~ x), col = "red")
        },
        labels = gsub("_","\n", colnames(check.wide[,-1])),
        upper.panel = panel.cor, 
        gap = 0,
        main = main
        )
}

#' Plot wide formated dataframe correlations per clone
#' in arc sine scale
#' @param ef.par wide formatted data by clone
#' @param main plot title  
#' @return \code{pairs} plot of correlation
#'
#' @examples
#' plot.wide.cor(efficacy.list(ef)$mean)

plot.ef.cor<-function(ef.par, main = "Correlation", ...){
  pairs(asin(ef.par[,-1]),
        lower.panel = points.panel,
        diag.panel = panel.hist,
        upper.panel = panel.cor,
        main = main)
}

  