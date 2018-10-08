---
title: "BIRM-CAM pilot"
author: "Steven Kiddle"
date: "8 October 2018"
output: html_document
---

```{r, }

# Load R packages, may need to install first (using install.packages())
library(data.table)
library(pheatmap)
library(grid)
library(PReMiuM)



# load data
load('../output/proc_data/CV.Rdata')

print('data loaded')

# extract 10,000 patients, and the columns (here 16:51) listing presence or absence of conditions you are studying
small <- copd_train[1:10000,c(6,16:51)]

# outcome column (five-year mortality) must be renamed outcome
colnames(small)[1] <- 'outcome'


set.seed(2)



# Format as required for PReMiuM:

covNames                <- names(small)[-1]

# 2: Put data and modelling options into an "inputs" data frame
inputs                  <- generateSampleDataFile(clusSummaryBernoulliDiscrete()) # This is just to initialise the "inputs" variable, so that it has all of the right fieldnames
inputs$inputData        <- small   #Put our data in the "inputData" slot
inputs$covNames         <- covNames
inputs$nCovariates      <- length(covNames)
inputs$fixedEffectNames <- NULL
inputs$xModel='Discrete'
inputs$discreteCovs <- covNames[2:36]
inputs$continuousCovs <- NULL


print('running')

# uncomment the runInfoObj lines to run clustering from scratch

# 3: Run PReMiuM!
# For the unsupervised case, we set "excludeY = TRUE"
# Note: nSweeps in the number of Gibbs iterations - might want to set to be a small number initially
#runInfoObj  <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=20000, nClusInit=100, nBurn=1, data=inputs$inputData, output="outcome_guided_20000", covNames = inputs$covNames, reportBurnIn = TRUE,excludeY=FALSE,discreteCovs = inputs$discreteCovs)

#runInfoObj2   <- profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=20000, nClusInit=100, nBurn=1, data=inputs$inputData, output="unsupervised_20000", covNames = inputs$covNames, reportBurnIn = TRUE,excludeY=TRUE,discreteCovs = inputs$discreteCovs)



# Load outcome-guided MCMC chain
z <- fread('outcome_guided_20000_z.txt')

# Load the posterior probability for each iteration of chain
logPost <- fread('outcome_guided_20000_logPost.txt')

# Identify Maximum A Posteriori (MAP) solution, i.e. the most probable one given the prior and data 
map1 <- which.max(logPost$V1)

# Extract the MAP solution, and count patients in each
map_s <- as.numeric(z[map1,])
table(map_s)

# Load multimorbidity data, and associated co-variates, looking only at first 10,000 patients
small <- copd_train[1:10000,]

# combine MAP clustering solution with this data
res <- cbind(small,map_s)

# for each cluster extract the average profile, age, five year mortality and number of co-morbidities
cluster_means <- matrix(NA,8,38)

for (i in 0:7){
  
  print(i)
  
  ind <- which(map_s == i)
  
  cluster_means[i+1,] <- c(mean(as.numeric(res[ind,'dead5']))-1,apply(res[ind,c(16:51)]=='present',2,'mean'),median(apply(res[ind,c(16:51)]=='present',1,'sum')))
  
}

# human readable versions of labels
var_labels <- c('Five_Year_Mortality','Chronic_Kidney_Disease','Alcohol_problems','Anorexia/bulimia','Atrial_fibrillation','Blindness/low_vision','Bronciectasis','Chronic_liver_disease','Sinusitis','Coronary_heart_disease','Dementia','Diabetes','Diverticular_disease','Hearing_loss','Heart_failure','Hypertension','Inflammatory_bowel_disease','Learning_disability','Multiple_sclerosis','Peripheral_vascular_disorder','Parkinsons','Prostate_disorders','Substance_abuse','Connective_tissue_disorders','Stroke','Thyroid_disorders','Anxiety','Asthma','Cancer','Constipation','Depression','Epilepsy','Painful_condition','Irritable_bowel_syndrome','Migrane','Psoriasis/eczema','Pyschosis/bipolar','Normalised_MM_burden')

# sort clusters based on risk of five year mortality
ord <- order(cluster_means[,1])

# make MM matrix for plotting
heatData <- t(cluster_means[ord,c(2:37)]*100)
rownames(heatData)<-var_labels[c(2:37)]
colnames(heatData) <- paste('O',1:8,sep='')

# find median age for each cluster
res <- as.data.table(res)
ages <- res[,list(median=median(age)),by='map_s']
setkey(ages,map_s)

# prepare extra data at top of heatmap, ordering clusters based on five year mortality
anno_col <- data.frame(Percent_patients=as.numeric(table(map_s)/100),Median_age=ages$median,cluster_means[,c(38,1)])[ord,]
anno_col[,2] <- anno_col[,2] + 67.13123
anno_col[,4] <- anno_col[,4] * 100
colnames(anno_col)[3:4] <- c('Median # co-morbidities',var_labels[1])
rownames(anno_col) <- paste('O',1:8,sep='')

# prepare colours for plot
myBlues             <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(100)
ann_colors = list(Five_Year_Mortality = rev(heat.colors(100)),'Median # co-morbidities' = cm.colors(100))

# Outcome-guided heatmap, will later use hierarchical cluster of conditions here to order them
a <- pheatmap::pheatmap(heatData,cluster_rows = T, cluster_cols = F, cellwidth = 12, 
                   cellheight = 8, fontsize = 8, color = myBlues, 
                  annotation_col = anno_col,annotation_colors = ann_colors)



# load unsupervised results, extracting MAP solution
z2 <- fread('unsupervised_20000_z.txt')
logPost2 <- fread('unsupervised_20000_logPost.txt')
map2 <- which.max(logPost2$V1)
map_u <- as.numeric(z2[map2,])
table(map_u)

# put solution with the data
res_u <- cbind(small,map_u)

# extract average MM profile, age, five year mortality and number of co-morbidities
cluster_means_u <- matrix(NA,7,38)

for (i in 0:6){
  
  print(i)
  
  ind <- which(map_u == i)
  
  cluster_means_u[i+1,] <- c(mean(as.numeric(res_u[ind,'dead5']))-1,apply(res_u[ind,c(16:51)]=='present',2,'mean'),median(apply(res_u[ind,c(16:51)]=='present',1,'sum')))
  
}

# the order of clusters five year mortality
ord_u <- order(cluster_means_u[,1])

# prepare for plotting
heatData_u <- t(cluster_means_u[ord_u,c(2:37)]*100)
rownames(heatData_u)<-var_labels[c(2:37)]
colnames(heatData_u) <- paste('U',1:7,sep='')

# get median ages
res_u <- as.data.table(res_u)
ages_u <- res_u[,list(median=median(age)),by='map_u']
setkey(ages_u,map_u)

# get data for top bars
anno_col_u <- data.frame(Percent_patients=as.numeric(table(map_u)/100),Median_age=ages_u$median,cluster_means_u[,c(38,1)])[ord_u,]
anno_col_u[,2] <- anno_col_u[,2] + 67.13123
anno_col_u[,4] <- anno_col_u[,4] * 100
colnames(anno_col_u)[3:4] <- c('Median # co-morbidities',var_labels[1])
rownames(anno_col_u) <- paste('U',1:7,sep='')

# add COPD as all patients have it, order the rest based on order from outcome-guided only plot
combinedPlusCOPD <- rbind(100,cbind(heatData,heatData_u)[a$tree_row$order,])
rownames(combinedPlusCOPD)[1] <- 'COPD'

pheatmap::pheatmap(combinedPlusCOPD,cluster_rows = F, cluster_cols = F, cellwidth = 12, 
         cellheight = 8, fontsize = 8, color = myBlues, 
         annotation_col = rbind(anno_col,anno_col_u),annotation_colors = ann_colors)

# produce plot, hierachical tree was added from the outcome-guided plot earlier, 
# powerpoint was used to add numbers and tidy up key colour bars
pdf('heatmap_pilot.pdf')
pheatmap::pheatmap(combinedPlusCOPD,cluster_rows = F, cluster_cols = F, cellwidth = 12, 
         cellheight = 8, fontsize = 8, color = myBlues, 
         annotation_col = rbind(anno_col,anno_col_u),annotation_colors = ann_colors)
dev.off()
```
