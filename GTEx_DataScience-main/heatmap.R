
#tracking which subjects correspond to which columns in the brain dataframe
brainsubtrack=rep(NA, dim(qn.brain.res)[2])
for (i in 1:dim(qn.brain.res)[2]){
  brainsubtrack[i]=strsplit(colnames(qn.brain.res)[i], "-")[[1]][2]
}

#tracking which subjects correspond to which columns in the colon dataframe
colonsubtrack=rep(NaN, dim(qn.colon.res)[2])
for (i in 1:dim(qn.colon.res)[2]){
  colonsubtrack[i]=strsplit(colnames(qn.colon.res)[i], "-")[[1]][2]
}

#tracking which subjects in the brain appear in the colon dataframe
intrack=rep(NA, length(brainsubtrack))
for (i in 1:length(brainsubtrack)){
  intrack[i]=brainsubtrack[i] %in% colonsubtrack
}

#removing subjects which do not appear
qn.brain.res.new=qn.brain.res[,intrack]

#renewing brainsubtracking according to new dataframe with removed columns
brainsubtrack=rep(NA, dim(qn.brain.res.new)[2])
for (i in 1:dim(qn.brain.res.new)[2]){
  brainsubtrack[i]=strsplit(colnames(qn.brain.res.new)[i], "-")[[1]][2]
}

#preparing a colon dataframe to correspond to the brain dataframe
qn.colon.res.braincomp=qn.colon.res[,which(colonsubtrack==brainsubtrack[1])[1]]
for (i in 2:length(brainsubtrack)){
  qn.colon.res.braincomp=cbind(qn.colon.res.braincomp,qn.colon.res[,which(colonsubtrack==brainsubtrack[i])[1]])
}

#calculating correlation matrix
colon_brain_heatmap=cor(t(qn.colon.res.braincomp), t(qn.brain.res.new))

#calculating the number of correlations whose absolute value is greater than 4 for each gene in brain sample
cols_genes_sum =colSums(abs(colon_brain_heatmap)>.4)

#sorting
cols_genes_sum_sorted=sort(cols_genes_sum, decreasing=TRUE)

#building gene rank dataframe according to the genes which correlate the most to colon samples
generankdf=data.frame(rep(NA, length(cols_genes_sum)))
colnames(generankdf)=c("rank")
generankdf$genedescription=rep(NA, length(cols_genes_sum))
generankdf$numcor=NA

for (i in 1:length(cols_genes_sum)){
  generankdf$rank[i]=i
  generankdf$genedescription[i]=gene.f[,"Description"][gene.f[,"Name"]==row.names(qn.brain.res.new)[order(cols_genes_sum, decreasing=TRUE)[i]]]
  generankdf$genename[i]=row.names(qn.brain.res.new)[order(cols_genes_sum, decreasing=TRUE)[i]]
}
generankdf$numcor=cols_genes_sum_sorted

# 10 target genes in the brain to be predicted
yframe=t(qn.brain.res.new[generankdf$genename[1:10],])

# genes in the colon - the features of the ML model 
modelX=t(qn.colon.res.braincomp)

# list all the genes in the colon that are correlated >|0.4| with the 10 target
# genes in the brain
connectgeneslists=list()
for (i in 1:10){
  connectgeneslists[[i]]=c(row.names(qn.colon.res.braincomp)[abs(colon_brain_heatmap[,generankdf$genename[i]])>.4])
}
