#XGBOOST
#install.packages("drat", repos="https://cran.rstudio.com")
#drat:::addRepo("dmlc"):
#  install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")
install.packages("InformationValue")#install.packages("caret")

library(caret)
library(xgboost)
library(InformationValue)
library(glmnet)

xgboost_model<-function(modelX.wgcna, gene_index){
  # Attach relevant target to features, according to gene index [1..10]
  dataML = data.frame(cbind(modelX.wgcna, yframe.wgcna[,gene_index]))
  
  set.seed(100)  # For reproducibility
  # Create index for testing and training data
  inTrain <- createDataPartition(y = dataML[,length(dataML)], p = 0.8, list = FALSE)
  # subset train data to training
  training <- dataML[inTrain,]
  # subset the rest to test
  testing <- dataML[-inTrain,]
  
  X_train = as.matrix(training[,-length(training)])
  y_train = training[,length(training)]
  X_test = as.matrix(testing[,-length(testing)])
  y_test = testing[,length(testing)]
  
  # train control to have cross validation with 5 fold
  xgb_trcontrol = trainControl(
    method = "cv",
    number = 5,  
    allowParallel = TRUE,
    verboseIter = FALSE,
    returnData = FALSE
  )
  
  # The hyperparameters to optimize 
  xgbGrid <- expand.grid(nrounds = c(15, 20, 25),  
                         max_depth = c(10, 15, 20),
                         colsample_bytree = seq(0.5, 0.9, length.out = 3),
                         eta = 0.1,
                         gamma=0,
                         min_child_weight = 1,
                         subsample = 1
  )
  
  # train the model
  set.seed(0) 
  xgb_model = train(
    X_train, y_train,  
    trControl = xgb_trcontrol,
    tuneGrid = xgbGrid,
    method = "xgbTree",
    verbose = FALSE,
    verbosity = 0
  )
  
  # print best model
  print(xgb_model$bestTune)
  
  #predict test target with test subset
  predicted = predict(xgb_model, X_test)
  
  # calculate errors: MSE R squared
  residuals = y_test - predicted
  MSE = mean(residuals^2)
  cat('The mean square error of the test data is ', round(MSE,3),'\n')
  
  y_test_mean = mean(y_test)
  # Calculate total sum of squares
  tss =  sum((y_test - y_test_mean)^2)
  # Calculate residual sum of squares
  rss =  sum(residuals^2)
  # Calculate R-squared
  rsq  =  1 - (rss/tss)
  cat('The R-square of the test data is ', round(rsq,3), '\n')
  
  return(list("mse" = MSE, "cor" = rsq))
}

# Generate empty list to populate in each iteration in the following for loop
# These lists will be a part of a dataframe we'll use to compare models
mselist=rep(NA,3)
corrlist=rep(NA,3)

lasso_mselist=rep(NA,3)
lasso_corrlist=rep(NA,3)

xgbcorr_mselist=rep(NA,3)
xgbcorr_corrlist=rep(NA,3)

# loop over all the 10 target genes in the brain and execute XGBoost with the following data:
# 1. All the genes in the colon
# 2. Genes selected by Lasso feature selection 
# 3. Highly correlated genes

for(i in 2:3){
  
  # XGBoost with all the genes in the colon
  mse_cor = xgboost_model(modelX.wgcna, i)
  mselist[i]=mse_cor$mse
  corrlist[i]=mse_cor$cor
  
  
  #lasso
  lambdas <- 10^seq(2, -2, by = -.1)
  
  # Setting alpha = 1 implements lasso regression
  lasso_reg <- cv.glmnet(modelX.wgcna, yframe.wgcna[,i], alpha = 1, standardize = TRUE, nfolds = 5)
  df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
  coef_genes = names(df_coef[which(abs(df_coef)>0),])[-1]
  # XGBoost with Lasso selected genes
  mse_cor_lasso = xgboost_model(modelX.wgcna[,coef_genes], i)
  lasso_mselist[i]=mse_cor_lasso$mse
  lasso_corrlist[i]=mse_cor_lasso$cor
  
  #xgb with correlated genes
  mse_cor_corr = xgboost_model(modelX.wgcna[,connectgeneslists.wgcna[[i]]], i)
  xgbcorr_mselist[i]=mse_cor_corr$mse
  xgbcorr_corrlist[i]=mse_cor_corr$cor
}

# translate results into dataframe
resultsdf1.wgcna=data.frame(generankdf.wgcna$genename[1:3])
rownames(resultsdf1.wgcna)=colnames(yframe.wgcna[1:3])
colnames(resultsdf1.wgcna)=c("description")
resultsdf1.wgcna$mse=mselist
resultsdf1.wgcna$rsq=corrlist
resultsdf1.wgcna$after_lasso_mse=lasso_mselist
resultsdf1.wgcna$after_lasso_rsq=lasso_corrlist
resultsdf1.wgcna$xgb_corr_mse=xgbcorr_mselist
resultsdf1.wgcna$xgb_corr_rsq=xgbcorr_corrlist
write.csv(resultsdf1.wgcna, "resultsdf_WGCNA.csv")

# Get Lasso selected genes for MEdarkolivegreen gene in the brain
gene.f_df_MEdarkolivegreen = as.data.frame(gene.f)
genes_MEdarkolivegreen = gene.f_df_MEdarkolivegreen[gene.f_df_MEdarkolivegreen$Name %in% coef_genes,]
write.csv(genes_MEdarkolivegreen, "MEdarkolivegreen_lasso_genes.csv")

# Get Lasso selected genes for MEivory gene in the brain
gene.f_df_MEivory = as.data.frame(gene.f)
genes_MEivory = gene.f_df_MEivory[gene.f_df_MEivory$Name %in% coef_genes,]
write.csv(genes_MEivory, "MEivory_lasso_genes.csv")

# get genes in olive & ivory modules
ivory_genes = names(datExpr)[moduleColors=="ivory"]
darkOliveGreen_genes = names(datExpr)[moduleColors=="darkolivegreen"]

# genes in ivory module
MEivory_genes_in_module = as.data.frame(gene.f)
genes_in_MEivory = MEivory_genes_in_module[MEivory_genes_in_module$Name %in% ivory_genes,]
write.csv(genes_in_MEivory, "MEivory_module_genes.csv")

# genes in olive module
MEdarkolivegreen_genes_in_module = as.data.frame(gene.f)
genes_in_MEdarkolivegreen = MEdarkolivegreen_genes_in_module[MEdarkolivegreen_genes_in_module$Name %in% darkOliveGreen_genes,]
write.csv(genes_in_MEdarkolivegreen, "MEdarkolivegreen_module_genes.csv")
