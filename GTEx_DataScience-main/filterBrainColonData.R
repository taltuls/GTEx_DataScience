# Look for mat.bs object, if doesn't exist -> run StartPreprocessData.R to generate it.
if(!exists("mat.bc")){
  source("StartPreprocessData.R")
}

# create a copy of mat.bs
reads.src1 = mat.bc

# transpose reads.src1
t.reads.src = t(reads.src1)

#delete genes with low values - with 80% of expression is 0.1.
vec.1 = apply(reads.src1 , 1, function(x) length(which( x > log(0.1+1, 2) )))
row.index = which(vec.1 > (0.8*(ncol(reads.src1))))

# leave just rows with expression at least 80% of the samples
src.reads = reads.src1 [row.index, ] #13227 genes

brain.cols = which(pheno.bc$SMTS %in% "Brain")
colon.cols = which(pheno.bc$SMTS %in% "Colon")

mat.brain = src.reads[,brain.cols]
mat.colon = src.reads[,colon.cols]

# Save processed data locally
save(mat.brain, file="mat.brain.Rdata")
save(mat.colon, file="mat.colon.Rdata")
# this is the filtered mat.bc 
save(src.reads, file="src.reads.Rdata")
save(mat.bc, file="mat.bc.Rdata")
save(pheno.bc, file="pheno.bc.Rdata")
save(pheno.colon, file="pheno.colon.Rdata")
save(pheno.brain, file="pheno.brain.Rdata")
#################################################################

#delete genes with variance = 0
var.data <- apply(src.reads, 1, var) #generate variance of each row - gene
low.var.indxs = which(var.data == 0)
if(length(low.var.indxs) > 0){
  data.free = src.reads
  #now we get smaller matrix, with no genes with variance 0
  src.reads <- data.free[-low.var.indxs,]
}
