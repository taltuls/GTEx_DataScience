if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("preprocessCore")
}

library(preprocessCore)

#rows are genes, columns are samples
quantile.normalize.raw.gtex <- function(edata.mat)
{
  norm_edata = normalize.quantiles(as.matrix(edata.mat))
  rownames(norm_edata) = rownames(edata.mat)
  colnames(norm_edata) = colnames(edata.mat)
  return(norm_edata)
}

############# COLON #############
# initial box plot
boxplot(filtered.mat.colon, main = "Colon boxplot before Quantile Normalization")
# quantile normalization
qn.colon = quantile.normalize.raw.gtex(filtered.mat.colon)
boxplot(qn.colon[], main = "Colon boxplot after Quantile Normalization")

############# BRAIN #############
# initial box plot
boxplot(filtered.mat.brain, main = "Brain boxplot before Quantile Normalization")
# quantile normalization
qn.brain = quantile.normalize.raw.gtex(filtered.mat.brain)
boxplot(qn.brain[], main = "Brain boxplot after Quantile Normalization")
