#!/usr/local/bin/Rscript
source('bugcountGLM.R')

# Parse command line arguments #################################################

# Usage:
# Rscript blup_heritability.R WF_consolidated.tab blup_analysis leaf

args <- commandArgs(TRUE)

if (is.na(args[1])) {
  args[1] <- 'WF_consolidated.tab'
}

wf.file <- args[1]

wf <- read.wf(file=wf.file, sep="\t",header=TRUE)

if(is.na(args[2])) {
  args[2] <- 'blup_analysis'
}

if(is.na(args[3])) {
  args[3] <- 'leaf'
}

wf.leaf <- leaf.count(wf)
summary(wf.leaf$nymphs)

wf.plant <- plant.count(wf)
summary(wf.plant$nymphs)
# args[3] <- "plant"

# Pick analysis level plant or leaf
if(args[3]=="plant"){
  wf.count <- wf.plant
} else if (args[3]=="leaf"){
  wf.count <- wf.leaf
}

# I need to make 4 sets of data per experimental unit (leaf or plant):


# Define dataset 1 #############################################################
# internal_check & offspring data
# infestation_check plants excluded excluded from the analysis.
# infestation_check nymph count instead should be used 
# as a measurement of infestation level

dataset1 <- wf.count[wf.count$group %in% c("internal_check","offspring"), ]
dataset1 <- remove.levels(dataset1)

print(unique(dataset1$clone[!grepl("CM|GM",dataset1$clone)]))

# Define dataset 2 #############################################################
# internal_check & offspring data excluding '2013_CM8996_invitro_poor'

dataset2 <- dataset1[dataset1$experiment != '2013_CM8996_invitro_poor', ]
dataset2 <- remove.levels(dataset2)

print(unique(dataset2$clone[!grepl("CM|GM",dataset2$clone)]))

# Define dataset 3 #############################################################
# Complete internal_check & offspring data including 2013_CM8996_invitro_poor
# ~ 40 CM8996 individuals 

dataset3 <- dataset1[dataset1$clone %in% select.complete.clones(dataset1), ]
dataset3 <- remove.levels(dataset3)

print(unique(dataset3$clone[!grepl("CM|GM",dataset3$clone)]))
unique(grep("CM",dataset3$clone, value = TRUE))
# Define dataset 4 #############################################################
# Complete internal_check & offspring data excluding 2013_CM8996_invitro_poor
# ~  CM8996 individuals 

dataset4 <- dataset2[dataset2$clone %in% select.complete.clones(dataset2), ]
dataset4 <- remove.levels(dataset4)

print(unique(dataset4$clone[!grepl("CM|GM",dataset4$clone)]))

# Dataset analysis ###########################################################

dataset <- list(dataset1, dataset2, dataset3, dataset4)

for (idx in 1:length(dataset)){
  # Initialize output per dataset ##############################################

  dataset.str <- paste("dataset", idx, sep = "")
  prefix <- paste(args[3],args[2],dataset.str, sep = "_")
  pdf.file <- paste(prefix,"pdf",sep = ".")
  pdf(file = pdf.file)
  cat(paste("\n", dataset.str, strrep("<",60), "\n\n"))
  
  # Negative Binomial GLM with random effects  #################################
  
  nb.model<- glmer.nb(nymphs ~ (1|experiment) + (1|clone) +
                                      (1|experiment:clone),
                             data = dataset[[idx]],
                             verbose=TRUE)

  # Most models will not converge to this point. 
  # Starting from this model and adding 200k (2e5) iterations to the optimizer 
  # usually makes the model converge
  
  start.pars <- getME(nb.model,c("theta","fixef"))
  nb.model <- update(nb.model, start=start.pars ,
                     control=glmerControl(optCtrl=list(maxfun=2e5)))
  
  # Print output in order to check for non-convergence
  print(summary(nb.model))
  print(max(nb.model@optinfo$derivs$gradient))
  
  cat(paste("\n",strrep(">",10), "glmer.nb warnings:\n"), sep = "")
  print(nb.model@optinfo$warnings)
  
  # check residuals with Dharma package
  # simulationOutput <- simulateResiduals( fittedModel = nb.model, n = 500, use.u = T)
  
  # # quartz()
  # plotSimulatedResiduals(simulationOutput = simulationOutput, quantreg = F)
  # Exits with error because there are no fixed effects. Plot interactively.
  
  # Calculate total BLUPs ######################################################
  # (including all experiments in the dataset)
  
  blup <- ranef(nb.model)
  blup.df <- data.frame(clone = rownames(blup$clone), blup = blup$clone[,1]) 
  blup.df <-blup.df[order(-blup.df$blup),]
  
  print(blup.df[!grepl("CM|GM",blup.df$clone),])
  
  
  # Calculate BLUPs by experiment ##############################################
  
  by.exp <- str_split_fixed(rownames(blup$`experiment:clone`), ":", 2)
  by.exp <- as.data.frame(by.exp)
  by.exp$blup <- blup$`experiment:clone`[[1]]
  colnames(by.exp) <- c("experiment","clone","blup")
  blup.by.exp<-dcast(by.exp, clone ~ experiment, value.var = "blup")
  # # quartz()
  pairs.panels(blup.by.exp[,-1])
  # # quartz()
  hist(blup$clone[,1], main = prefix)
  
  # Make output of BLUPs as phenotypes for MapQTL  #############################
  # combining total BLUP and BLUP by experiment
  
  out <- merge(blup.df, blup.by.exp, by="clone", all.x = TRUE)
  out <- merge(x = data.frame(clone = levels(dataset1$clone)),
               y = out, by="clone", all.x = TRUE)
  write.table(out, file = paste(prefix,"tab",sep = "."), 
              row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Comparing offspring BLUP with parentals  ###################################
  
  gm <-blup.df[grep("GM8586",blup.df$clone),]
  nrow(gm)
  gm[ gm$blup >= blup.df$blup[blup.df$clone == "TMS60444"],] 
  nrow(gm[ gm$blup >= blup.df$blup[blup.df$clone == "TMS60444"], ])
  gm[ gm$blup <= blup.df$blup[blup.df$clone == "ECU72"], ]
  nrow(gm[ gm$blup <= blup.df$blup[blup.df$clone == "ECU72"], ])
  
  cm <- blup.df[grep("CM8996",blup.df$clone),]
  nrow(cm)
  cm[ cm$blup >= blup.df$blup[blup.df$clone == "COL1468"], ]
  nrow(cm[ cm$blup >= blup.df$blup[blup.df$clone == "COL1468"], ])
  cm[ cm$blup <= blup.df$blup[blup.df$clone == "ECU72"], ]
  nrow(cm[ cm$blup <= blup.df$blup[blup.df$clone == "ECU72"], ])

  # Plotting distribution by cross  ########################################  

  # quartz()
  hist(cm$blup, main = paste("CM8996",prefix))
  
  # quartz()
  hist(gm$blup, main = paste("GM8586",prefix))

  # Checking fit of BLUP with nymph counts #####################################
  
  nymph.means <- aggregate(nymphs~clone, 
                           data = dataset1, FUN=mean, na.rm=TRUE)
  
  nymph.merge  <- merge(nymph.means,blup.df, by="clone")
  
  cross.type <- 1 * grepl("CM8996", nymph.merge$clone) +
                2 * grepl("GM8586", nymph.merge$clone) +
                3 * !grepl("CM8996|GM8586", nymph.merge$clone)
  nymph.merge$cross <- factor(c("CM8996","GM8586", "check")[cross.type],
                                 levels = c("CM8996","GM8586", "check"))
# quartz()
  print( ggplot(nymph.merge) + ggtitle(prefix) +
    geom_point(aes(x=blup, y=nymphs, color=cross), shape=1) +
    geom_text(data=subset(nymph.merge, cross == "check"),
            aes(x=blup, y=nymphs,label=clone, color=cross),
            hjust = 0, nudge_x = 0.005,
            show.legend = FALSE))

#  quartz()
  print( ggplot(nymph.merge) + ggtitle(prefix) +
        geom_point(aes(x=blup, y=nymphs+1, color=cross), shape=1) +
        geom_text(data=subset(nymph.merge, cross == "check"),
                    aes(x=blup, y=nymphs,label=clone, color=cross),
                    hjust = 0, nudge_x = 0.005,
                    show.legend = FALSE) +
        scale_y_continuous(trans=log2_trans()))


  print(summary(lm(nymphs ~ blup,data=nymph.merge)))
  

  graphics.off()
  
}





