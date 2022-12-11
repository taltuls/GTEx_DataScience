load("mat.colon.RData")
load("pheno.colon.RData")
if (!require("devtools")) install.packages("devtools")
library(devtools)
#library(ggplot2)
#install.packages("ggfortify")
library(ggfortify)
# categorize colon tissues = SMTSD

# PCA with mat colon
colonPCA = prcomp(t(mat.colon))
autoplot(colonPCA, data = pheno.colon, colour = 'SMTSD',main = "PCA for all colon samples")

# non relevant colon tissues = cerebellum, cerebellar hemisphere
colonFilter = which(pheno.colon$SMTSD %in% c('Colon - Cerebellar Hemisphere', 'Colon - Cerebellum'))
filtered.pheno.colon = pheno.colon[-colonFilter,] #285
filtered.colon.samples = as.vector(pheno.colon[-colonFilter,'SAMPID'])
filtered.mat.colon = mat.colon[,filtered.colon.samples] #285

# plot pca with color categories = SMRIN
colonPCA = prcomp(t(filtered.mat.colon))
autoplot(colonPCA, data = pheno.colon, colour = 'SMRIN', 
         main = "Samples Quality")

# plot pca with color categories = DTHHRDY
autoplot(colonPCA, data = pheno.colon, colour = 'DTHHRDY', 
         main = "Cause of death", labels = colnames(mat.colon))

# plot pca with color categories = GENDER
autoplot(colonPCA, data = filtered.pheno.colon, colour = 'GENDER', 
         main = "Gender")

# plot pca with color categories = AGE
autoplot(colonPCA, data = filtered.pheno.colon, colour = 'AGE', 
         main = "Age")
