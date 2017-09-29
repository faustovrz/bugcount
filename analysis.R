#!/usr/local/bin/Rscript
source('bugcountGLM.R')

# Parse command line arguments #################################################

# Usage:
# Rscript anlysis.R wf_consolidated.tab resistance_analysis.pdf leaf

args <- commandArgs(TRUE)

if (is.na(args[1])) {
  args[1] <- 'leaf'
}

if(is.na(args[2])) {
    args[2] <- 'resistance_analysis'
}

if(is.na(args[3])) {
  args[3] <- 'wf_consolidated.tab'
}


# Whitefly analysis.

# Start pdf output device for all subsequent plots

prefix <- paste(args[1],args[2], sep = "_")

pdf(paste(prefix,"pdf", sep = "."), onefile=TRUE)

# 1. Read WF count data. #######################################################
#
# WF counts per picture. Each picture must have a unique combination of
# experiment + replicate + plant + pic
# experiment description consists of:
# exp.year + exp.cross + exp.propagation + exp.substrate
# **************************************************************************** #

wf.file <- args[3]

wf <- read.wf(file=wf.file, sep="\t",header=TRUE)

# 2. Counts per Plant / per Leaf. ##############################################

wf.leaf <- leaf.count(wf)
summary(wf.leaf$nymphs)

wf.plant <- plant.count(wf)
summary(wf.plant$nymphs)

# Pick analysis level plant or leaf
if(args[1]=="plant"){
wf.count <- wf.plant
} else if (args[1]=="leaf"){
wf.count <- wf.leaf
}

# Leave infestation check out?

# **************************************************************************** #
#
#  What indicator should we use as resistance phenotype?
#
#  Nymph counts are not normal 
#  Error is not distributed normally.
#  Nymph counts do not meet linear model assumptions (ANOVA, MapQTL).
#  Best fit is a negative Binomial Distribution
#
# **************************************************************************** #

# Probability density model fit to actual distribution

wf.count.fit <-  count.fit(wf.count)

# quartz()
plot.count.fit(wf.count.fit, main = "Nymphs")


# 3. Reproducibility ##########################################################

# Assume that nymphs per plant are distributed as negative binomial

# **************************************************************************** #
# 3.1 Whole population.


# nymphs ~ experiment GLM, posthoc and common letter difference 
by.experiment <- fit.nb.glm(nymphs ~ experiment,wf.count)

# quartz()
plot.fit.nb.glm(wf.count, by.experiment)


# quartz()
plot.fit.nb.glm(wf.count,by.experiment,type="density")

# add infestation regime as factor 
#    high,low: inferred from nyphs ~ experiment GLM

wf.count$infestation <- "low"
wf.count[wf.count$exp.year == 2017,]$infestation <- "high"
wf.count$infestation <- factor(wf.count$infestation, levels =c("low","high"))

# I need a single normally distributed number per cassava accesion
# for QTL analysis

# log transfromed geometric mean 
# Robust against extreme values 

# **************************************************************************** #
# Correlation  between experiments 

cor.title <- paste("Geometric Means Correlation","Log Scale", sep = "\n")

for (cross in levels(wf.count$exp.cross)) {
  
  wf.by.exp <- wf.wide(clone~experiment,
                       wf.count[wf.count$exp.cross==cross,])
  
#   # quartz()
  
  plot.wide.cor(wf.by.exp, main = paste(cross,cor.title))
  
}



# **************************************************************************** #
# 3.2 Checks.

# Leave infestation check out?
# wf.check <- wf.count[grepl("check", wf.count$group) &
wf.check <- wf.count[grepl("internal_check", wf.count$group) &
                    wf.count$clone != "Secundina",]

wf.check <- remove.levels(wf.check)
names(wf.check)
levels(wf.check$clone)

# Correlation  between experiments

check.by.exp <- wf.wide(clone ~ experiment, wf.check, 
                        fun = function(x) log10(geometric.mean(x)) )
levels(wf.check$clone)

# quartz()
plot.check.cor(check.by.exp, main = paste("Checks", cor.title))

head(wf.check)
# nymphs ~ clone GLM posthoc and common letter difference
fit <- fit.nb.glm(nymphs ~ clone, wf.check)

# quartz()
plot.fit.nb.glm(wf.check,fit)

# quartz()
plot.fit.nb.glm(wf.check,fit,type = "density")


# nymphs ~ infestation * clone GLM

fit <- glm.nb(nymphs ~  infestation * clone, data = wf.check)
# lsmeans(fit, ~ infestation * clone,
#         adjust = "tuckey", transform = "response") 
# lsmeans with transorm = "response" might give geometric means
# and geomeric means ratio tests and confidence intervals I need to confirm!!!

posthoc <- cld(lsmeans(fit, ~ infestation * clone, adjust = "tuckey"),
               type = "response")


# Plot post hoc

### Order the levels for printing
clone.order <- wf.aggregate(wf.check, by.stat = "gmean")$clone

posthoc$clone <- factor(posthoc$clone,
                      levels=clone.order)

posthoc$infestation = factor(posthoc$infestation,
                   levels=c("low", "high"))

###  Remove spaces in .group  

posthoc$.group=gsub(" ", "", posthoc$.group)

# quartz()
plot.gg.infestationXclone(posthoc)


# quartz()
plot.infestationXclone(wf.check)


# Nymph counts by Cross #########################################################

exp.levels <-levels(as.factor(wf.count$experiment))
exp.grp <-list(CM8996 = list(1,2,4,5,c(1,2),c(2,4),c(4,5),c(1,2,4), c(2,4,5)),
               GM8586 = list(3,6,c(3,6)))

# exp.grp <-list(CM8996 = list(1))


for (cross in levels(wf.count$exp.cross)){
  for (idx in 1:length(exp.grp[[cross]])){
  exp.allowed <- exp.levels[exp.grp[[cross]][[idx]]]
  
  wf.by.cross <- wf.count[wf.count$exp.cross==cross &
                            wf.count$experiment %in% exp.allowed,]
  
  exp.count <- aggregate(experiment ~ clone ,wf.by.cross,
                         FUN = function(x) length(unique(x)))
  selected.clones <- exp.count$clone[exp.count$experiment==length(exp.allowed)]
  wf.by.cross <- wf.by.cross[ wf.by.cross$clone %in%selected.clones,]
  wf.by.cross <- remove.levels(wf.by.cross)
  
  wf.by.cross$nymphs<-wf.by.cross$nymphs+1
  
  # # quartz()
  # hist(log10(wf.by.cross$nymphs))
  
  # by.cross.count.fit <-  count.fit(wf.by.cross[wf.by.cross$group == "offspring",])
  
  # # quartz()
  # plot.count.fit(by.cross.count.fit, 
  #                main = paste(c(exp.allowed,"Nymphs"), collapse ="\n"))
  
  
  # ANOVA ********************************************************************* #
  # Assuming Normal distribution
  #
  
  fit.lm <- lm(nymphs ~  clone, data = wf.by.cross)
  AIC(fit.lm)
  posthoc <- cld(lsmeans(fit.lm, "clone", adjust = "tuckey"))
  colnames(posthoc)
#   # quartz()
  plot.cross.anova(posthoc,exp.allowed)
  
  # GLM ********************************************************************** #
  # Assuming Negative Binomial distribution
  #
  
  fit.nb <- glm.nb(nymphs ~  clone, data = wf.by.cross)
  AIC(fit.nb)
  
  posthoc <- cld(lsmeans(fit.nb, ~ clone, adjust = "tuckey"),
                 type = "response")
  
  # Plot post hoc
#   # quartz()
  plot.cross.nb.fit(posthoc,exp.allowed)
  
  # Compare models
  anova(fit.lm,fit.nb)
  }
}



# 4. Ressistance Efficacy #####################################################

# Calculate ressitance efficacy per experiment
# and Correlation between experiments

control <- "COL1468"

ef<- efficacy.aggregate(wf.count, control = control)

# Clone Sample Size by Experiment

# quartz()
par(mar=c(12,4,4,2)) 
boxplot(clone_n ~ experiment, data = ef, 
        main = "Clone Sample Size by Experiment",
        las = 2, log ="y")

# quartz()
plot.ef.density(ef,"geometric.mean.eff","experiment")

# quartz()
par(mar=c(15,4,0,2)) 
boxplot(asin(geometric.mean.eff) ~ experiment, las = 2,
        ylab = "asin(Efficacy)", data=ef)

# Correlation  between experiments

cor.title <- paste("Efficacy Correlation","asin Transform", sep = "\n")

for (cross in levels(ef$exp.cross)) {

# remove outlying point for correlations
ef.by.cross <- ef[ef$exp.cross == cross & ef$clone != control,] 
ef.by.cross <-remove.levels(ef.by.cross)

ef.gmean <-efficacy.list(ef.by.cross)$geometric.mean

# quartz()
plot.ef.cor(ef.gmean,main = paste(cross,"/",control,cor.title))

}

graphics.off()

# 5. Sample size calculations #################################################

