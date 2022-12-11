# hirarchy clustering with brain mat data
sampleTree = hclust(dist(t(mat.brain)), method = "average")
par(cex = 0.3)
par(mar = c(0,4,2,0))

# plot the cluster
plot(sampleTree, main = "Brain Clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# calculate cut average for k = 2, k is 2 because we know there are two different groups 
# of samples in the brain, as we can see from the dendogram.
cut_avg <- cutree(sampleTree, k = 2)
plot(sampleTree)
rect.hclust(sampleTree , k = 2, border = 2:6)

# Look at cluster results
ind1 = which(cut_avg == 1)
ind2 = which(cut_avg == 2)

# Get samples in cluster 1 and filter brain and mat data accordingly
brain.clust1.indxs = which(pheno.brain$SAMPID %in% names(ind1))
filtered.pheno.brain = pheno.brain[brain.clust1.indxs,] # 285
filtered.mat.brain = mat.brain[,brain.clust1.indxs]

# get samples in cluster 2 and ignore it - we don't need these samples
# these samples are outliers
brain.clust2.indxs = which(pheno.brain$SAMPID %in% names(ind2))
out.pheno.brain = pheno.brain[brain.clust2.indxs,] # 68
