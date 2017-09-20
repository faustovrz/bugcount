# Whitefly analysis.

source('bugcountGLM.R')

# Start pdf output device for all subsequent plots

# pdf.file <- 'resistance_analysis.pdf'
# pdf(pdf.file, onefile=TRUE)

# 1. Read WF count data. #######################################################
#
# WF counts per picture. Each picture must have a unique combination of
# experiment + replicate + plant + pic
# experiment description consists of:
# exp.year + exp.cross + exp.propagation + exp.substrate
# **************************************************************************** #

wf <- read.wf(file="WF_consolidated.tab", sep="\t",header=TRUE)

# 2. Counts per Plant / per Leaf. ##############################################

wf.leaf <- leaf.count(wf)
summary(wf.leaf$nymphs)

wf.plant <- plant.count(wf)
summary(wf.plant$nymphs)

# Pick analysis level plant or leaf

# wf.count <- wf.plant

wf.count <- wf.leaf

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

quartz()
plot.count.fit(wf.count.fit, main = "Nymphs")


# 3. Reproducibility ##########################################################

# Assume that nymphs per plant are distributed as negative binomial

# **************************************************************************** #
# 3.1 Whole population.


# nymphs ~ experiment GLM, posthoc and common letter difference 
by.experiment <- fit.nb.glm(nymphs ~ experiment,wf.count)

quartz()
plot.fit.nb.glm(wf.count, by.experiment)


quartz()
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
  
  quartz()
  
  plot.wide.cor(wf.by.exp, main = paste(cross,cor.title))
  
}



# **************************************************************************** #
# 3.2 Checks.

wf.check <- wf.count[grepl("check", wf.count$group) &
                    wf.count$clone != "Secundina",]

wf.check <- remove.levels(wf.check)
names(wf.check)
levels(wf.check$clone)

# Correlation  between experiments

check.by.exp <- wf.wide(clone ~ experiment, wf.check, 
                        fun = function(x) log10(geometric.mean(x)) )
levels(wf.check$clone)

quartz()
plot.check.cor(check.by.exp, main = paste("Checks", cor.title))

head(wf.check)
# nymphs ~ clone GLM posthoc and common letter difference
fit <- fit.nb.glm(nymphs ~ clone, wf.check)

quartz()
plot.fit.nb.glm(wf.check,fit)

quartz()
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

quartz()
plot.gg.infestationXclone(posthoc)


quartz()
plot.infestationXclone(wf.check)


# Nymph counts by Cross #########################################################
# Select just the clones where there are more than 2 experiments
cross <-"CM8996"
exp.levels <-levels(as.factor(wf.count$experiment))
exp.grp <-list(1,2,4,5,c(1,2),c(2,4),c(4,5),c(1,2,4), c(2,4,5))
for (idx in 1:length(exp.grp)){
exp.allowed <- exp.levels[exp.grp[[idx]]]

wf.by.cross <- wf.count[wf.count$exp.cross==cross &
                          wf.count$experiment %in% exp.allowed,]

exp.count <- aggregate(experiment ~ clone ,wf.by.cross,
                       FUN = function(x) length(unique(x)))
selected.clones <- exp.count$clone[exp.count$experiment==length(exp.allowed)]
wf.by.cross <- wf.by.cross[ wf.by.cross$clone %in%selected.clones,]
wf.by.cross <- remove.levels(wf.by.cross)

wf.by.cross$nymphs<-wf.by.cross$nymphs+1

# quartz()
# hist(log10(wf.by.cross$nymphs))

by.cross.count.fit <-  count.fit(wf.by.cross[wf.by.cross$group == "offspring",])

# quartz()
# plot.count.fit(by.cross.count.fit, 
#                main = paste(c(exp.allowed,"Nymphs"), collapse ="\n"))


# ANOVA ********************************************************************* #
# Assuming Normal distribution
#
tail(wf.by.cross)
fit.lm <- lm(nymphs ~  clone, data = wf.by.cross)
AIC(fit.lm)
posthoc <- cld(lsmeans(fit.lm, "clone", adjust = "tuckey"))
colnames(posthoc)

### Order the levels for printing
# clone.order <- wf.aggregate(wf.by.cross, by.stat = "gmean")$clone
clone.order <- posthoc$clone
posthoc$clone <- factor(posthoc$clone,
                      levels=clone.order)
###  Remove spaces in .group  

posthoc$.group=gsub(" ", "", posthoc$.group)

lm.post <- posthoc

# quartz()
# print(ggplot(lm.post,
#          aes(x     = clone,
#              y     = lsmean,
#              color = .group,
#              label = .group),
#          log10="y") +
#     geom_point(shape  = 15,
#                size   = 4) +
#     geom_errorbar(aes(ymin  =  lower.CL, 
#                       ymax  =  upper.CL),
#                   width =  0.2,
#                   size  =  0.7) +
#     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#     #      legend.position =  c(0.9, 0.2) ) +
#     # scale_y_continuous(trans=log10_trans(), limits = c(1,10000) ) +
#     ggtitle(paste(c(exp.allowed,"ANOVA, Mean + 95% Confidence Interval"),
#                   collapse ="\n")))


# GLM ********************************************************************** #
# Assuming Negative Binomial distribution
#

fit.nb <- glm.nb(nymphs ~  clone, data = wf.by.cross)
AIC(fit.nb)
anova(fit.lm,fit.nb)

posthoc <- cld(lsmeans(fit.nb, ~ clone, adjust = "tuckey"),
               type = "response")
colnames(posthoc)

# Plot post hoc
# plot(posthoc)

### Order the levels for printing
# clone.order <- wf.aggregate(wf.by.cross, by.stat = "gmean")$clone
clone.order <- posthoc$clone
posthoc$clone <- factor(posthoc$clone,
                      levels=clone.order)

###  Remove spaces in .group  

posthoc$.group=gsub(" ", "", posthoc$.group)

nb.glm.post <- posthoc

quartz()
print(ggplot(nb.glm.post,
         aes(x     = clone,
             y     = response,
             color = .group,
             label = .group),
         log10="y") +
    geom_point(shape  = 15,
               size   = 4) +
    geom_errorbar(aes(ymin  =  asymp.LCL, 
                      ymax  =  asymp.UCL),
                  width =  0.2,
                  size  =  0.7) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          legend.position="none") +
    scale_y_continuous(trans=log10_trans(),limits = c(0.001,100000)) +
    ggtitle(paste(c(exp.allowed,
                  "Negative Binomial GLM, Mean + 95% Confidence Interval"),
                  collapse ="\n")))
}

p <- merge(lm.post, nb.glm.post, by="clone")

quartz()
plot(response~lsmean, data=p)



# 4. Ressistance Efficacy #####################################################

# Calculate ressitance efficacy per experiment
# and Correlation between experiments

control <- "COL1468"

ef<- efficacy.aggregate(wf.count, control = control)

# Clone Sample Size by Experiment

quartz()
par(mar=c(12,4,4,2)) 
boxplot(clone_n ~ experiment, data = ef, 
        main = "Clone Sample Size by Experiment",
        las = 2, log ="y")

quartz()
plot.ef.density(ef,"geometric.mean.eff","experiment")

quartz()
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

quartz()
plot.ef.cor(ef.gmean,main = paste(cross,"/",control,cor.title))

}

# graphics.off()

# 5. Sample size calculations #################################################

