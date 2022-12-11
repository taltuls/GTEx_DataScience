load("mat.brain.RData")
load("pheno.brain.RData")
if (!require("devtools")) install.packages("devtools")
library(devtools)
library(ggplot2)
install.packages("ggfortify")
library(ggfortify)
# categorize brain tissues = SMTSD

# PCA with mat brain
brainPCA = prcomp(t(mat.brain))
autoplot(brainPCA, data = pheno.brain, colour = 'SMTSD',main = "PCA for all brain samples")

# non relevant brain tissues = cerebellum, cerebellar hemisphere
brainFilter = which(pheno.brain$SMTSD %in% c('Brain - Cerebellar Hemisphere', 'Brain - Cerebellum'))
filtered.pheno.brain = pheno.brain[-brainFilter,] #285
filtered.brain.samples = as.vector(pheno.brain[-brainFilter,'SAMPID'])
filtered.mat.brain = mat.brain[,filtered.brain.samples] #285

# plot pca with color categories = SMTSD
brainPCA = prcomp(t(filtered.mat.brain))
autoplot(brainPCA, data = filtered.pheno.brain, colour = 'SMTSD', 
         main = "PCA w/o Cerebellar Hemisphere and Cerrebellum")

# plot pca with color categories = SMRIN
brainPCA = prcomp(t(filtered.mat.brain))
autoplot(brainPCA, data = filtered.pheno.brain, colour = 'SMRIN', 
         main = "Samples Quality")

# plot pca with color categories = DTHHRDY
autoplot(brainPCA, data = filtered.pheno.brain, colour = 'DTHHRDY', 
         main = "Cause of death")

# plot pca with color categories = GENDER
autoplot(brainPCA, data = filtered.pheno.brain, colour = 'GENDER', 
         main = "Gender")

# plot pca with color categories = AGE
autoplot(brainPCA, data = filtered.pheno.brain, colour = 'AGE', 
         main = "Age")
