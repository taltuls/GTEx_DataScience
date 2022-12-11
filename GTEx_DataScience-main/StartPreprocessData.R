# Load the data
load("gene.f.RData")
load("gene.f.with.entrez.RData")
load("mat.f.coding.RData")
load("pheno.f.RData")

# Get brain and colon data from pheno
pheno.brain <- pheno.f[pheno.f['SMTS'] == "Brain",] # 1259
pheno.colon <- pheno.f[pheno.f['SMTS'] == "Colon",] # 345

#initialize empty lists to insert into the pairs
samplesBrain <- c()
samplesColon <- c()
SubjId <- c()

# Get SampleIDs for brain and colon
sampBrain <- as.vector(pheno.brain$SAMPID) # 1259 sampleID's
sampColon <- as.vector(pheno.colon$SAMPID) # 345 sampleID's

# Get SubjectIDs for brain and colon
subjBrain <- as.vector(pheno.brain$SUBJID) # 1259 sujectID's
subjColon <- as.vector(pheno.colon$SUBJID) # 345 subjectID's

# Loop to find pairs and insert them into the samplesBrain, samplesColon, SubjId lists
for (i in 1:nrow(pheno.brain)) {
  for (j in 1:nrow(pheno.colon)) {
    if (subjBrain[i] == subjColon[j]) {
      samplesBrain <- append(samplesBrain, sampBrain[i])
      samplesColon <- append(samplesColon, sampColon[j])
      SubjId <- append(SubjId, subjColon[j])
    }
  }
}

# Get unique subject ids in the pairs
SubjUnique <- unique(SubjId)
bcSampsID <- append(unique(samplesBrain), unique(samplesColon))

# Create a subset of colon and brain pheno
pheno.bc <- pheno.f[pheno.f$SAMPID %in% bcSampsID,]
# Create a subset of brain pheno
pheno.brain <- pheno.brain[pheno.brain$SAMPID %in% unique(samplesBrain),] # 353
# Create a subset of colon pheno
pheno.colon <- pheno.colon[pheno.colon$SAMPID %in% unique(samplesColon),] # 66

# Clear unused data, helps working with big files. 
gc()
memory.limit(10000)

# Create a subset of colon and brain mat
mat.bc = mat.f.coding[,bcSampsID]
# Create a subset of brain mat
mat.brain = mat.f.coding[,unique(samplesBrain)]
# Create a subset of colon mat
mat.colon = mat.f.coding[,unique(samplesColon)]

rm(subjBrain, subjColon, mat.f.coding)
gc()
