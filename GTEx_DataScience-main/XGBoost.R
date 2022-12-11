#XGBOOST
#install.packages("drat", repos="https://cran.rstudio.com")
#drat:::addRepo("dmlc"):
#  install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")
install.packages("InformationValue")#install.packages("caret")

library(caret)
library(xgboost)
library(InformationValue)
library(glmnet)

xgboost_model<-function(modelX, gene_index){
  # Attach relevant target to features, according to gene index [1..10]
  dataML = data.frame(cbind(modelX, yframe[,gene_index]))
  
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
mselist=rep(NA,10)
corrlist=rep(NA,10)

lasso_mselist=rep(NA,10)
lasso_corrlist=rep(NA,10)

xgbcorr_mselist=rep(NA,10)
xgbcorr_corrlist=rep(NA,10)

# loop over all the 10 target genes in the brain and execute XGBoost with the following data:
# 1. All the genes in the colon
# 2. Genes selected by Lasso feature selection 
# 3. Highly correlated genes

for(i in 1:10){

  # XGBoost with all the genes in the colon
  mse_cor = xgboost_model(modelX, i)
  mselist[i]=mse_cor$mse
  corrlist[i]=mse_cor$cor
  
  
  #lasso
  lambdas <- 10^seq(2, -3, by = -.1)
  
  # Setting alpha = 1 implements lasso regression
  lasso_reg <- cv.glmnet(modelX, yframe[,i], alpha = 1, standardize = TRUE, nfolds = 5)
  df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
  coef_genes = names(df_coef[which(abs(df_coef)>0),])[-1]
  # XGBoost with Lasso selected genes
  xgboost_model(modelX[,coef_genes], i)
  mse_cor_lasso = xgboost_model(modelX[,coef_genes], i)
  lasso_mselist[i]=mse_cor_lasso$mse
  lasso_corrlist[i]=mse_cor_lasso$cor
  
  #xgb with correlated genes
  mse_cor_corr = xgboost_model(modelX[,connectgeneslists[[i]]], i)
  xgbcorr_mselist[i]=mse_cor_corr$mse
  xgbcorr_corrlist[i]=mse_cor_corr$cor
}

# translate results into dataframe
resultsdf1=data.frame(generankdf$genedescription[1:10])
rownames(resultsdf1)=colnames(yframe)
colnames(resultsdf1)=c("description")
resultsdf1$mse=mselist
resultsdf1$rsq=corrlist
resultsdf1$after_lasso_mse=lasso_mselist
resultsdf1$after_lasso_rsq=lasso_corrlist
resultsdf1$xgb_corr_mse=xgbcorr_mselist
resultsdf1$xgb_corr_rsq=xgbcorr_corrlist
write.csv(resultsdf1, "resultsdf_fixed.csv")

# Get Lasso selected genes for FKBP5 gene in the brain
gene.f_df_fkb = as.data.frame(gene.f)
genes_fkbp = gene.f_df_fkb[gene.f_df_fkb$Name %in% coef_genes,]
write.csv(genes_fkbp, "FKBP5_lasso_genes.csv")

# Get Lasso selected genes for SLCO4A1 gene in the brain
gene.f_df_sloc = as.data.frame(gene.f)
genes_sloc = gene.f_df_sloc[gene.f_df_sloc$Name %in% coef_genes,]
write.csv(genes_sloc, "SLCO4A1_lasso_genes.csv")


# which(pheno.f$Name %in% coef_genes)
# gene.f[,1]
# gene.f[gene.f$]
# 
# xcorr = train[,connectgeneslists[[i]]]
# mse_cor = run_xgboost(train[,connectgeneslists[[i]]],test[,connectgeneslists[[i]]],i)
# xgbcorr_mselist[i]=mse_cor$mse
# xgbcorr_corrlist[i]=mse_cor$cor













# 
# 
# 
# 
# 
# set.seed(100)
# max_depth = 4
# num_rounds = 8
# 
# run_xgboost <- function(train, test, i){                   data               label                                                                                  "binary:logistic"
#   bstSparse <- xgboost(data = as.matrix(train[,-ncol(train)]), label = as.numeric(train[,ncol(train)]), max.depth = max_depth, eta = 1.5, nthread = 2, nrounds = num_rounds, objective = "reg:squarederror")
#   pred <- predict(bstSparse, as.matrix(test[,-ncol(test)]))
#   
#   pred
#   #                      test labels
#   err <- mean(pred - test[,ncol(test)])^2
#   #corrlist[i]=cor(pred,test[,ncol(test)])
#   
#   print(paste("mse=", err))
#   #mselist[i]=err
#   
#   return(list("mse" = err, "cor" = cor(pred,test[,ncol(test)])))
# }
# 
# 
# 
# mselist=rep(NA,10)
# corrlist=rep(NA,10)
# 
# lasso_mselist=rep(NA,10)
# lasso_corrlist=rep(NA,10)
# 
# xgbcorr_mselist=rep(NA,10)
# xgbcorr_corrlist=rep(NA,10)
# 
# split1 <- sample(c(rep(0, 0.7 * nrow(newdf)), rep(1, 0.3 * nrow(newdf))))
# 
# for(i in 1:10){
#   newdf = data.frame(cbind(modelX, yframe[,i]))
#   colnames(newdf)[length(colnames(newdf))]="label"
#   train = newdf[split1==0,]
#   test = newdf[split1==1,]
#   
#   mse_cor = run_xgboost(train,test,i)
#   mselist[i]=mse_cor$mse
#   corrlist[i]=mse_cor$cor
#   
#   #lasso
#   lambdas <- 10^seq(2, -3, by = -.1)
#   
#   # Setting alpha = 1 implements lasso regression
#   lasso_reg <- cv.glmnet(modelX, yframe[,i], alpha = 1, standardize = TRUE, nfolds = 5)
#   #cat('Min Lambda: ', lasso_reg$lambda.min, '\n 1Sd Lambda: ', lasso_reg$lambda.1se)
#   df_coef <- round(as.matrix(coef(lasso_reg, s=lasso_reg$lambda.min)), 2)
#   coef_genes = names(df_coef[which(abs(df_coef)>0),])[-1]
#   x = train[,coef_genes]
#   #xgboost 2
#   mse_cor = run_xgboost(train[,coef_genes],test[,coef_genes],i)
#   lasso_mselist[i]=mse_cor$mse
#   lasso_corrlist[i]=mse_cor$cor
#   
#   #xgb with correlated genes
#   xcorr = train[,connectgeneslists[[i]]]
#   mse_cor = run_xgboost(train[,connectgeneslists[[i]]],test[,connectgeneslists[[i]]],i)
#   xgbcorr_mselist[i]=mse_cor$mse
#   xgbcorr_corrlist[i]=mse_cor$cor
# }
# resultsdf=data.frame(generankdf$genedescription[1:10])
# rownames(resultsdf)=colnames(yframe)
# colnames(resultsdf)=c("description")
# resultsdf$mse=mselist
# resultsdf$rsq=corrlist*corrlist
# resultsdf$after_lasso_mse=lasso_mselist
# resultsdf$after_lasso_rsq=lasso_corrlist*lasso_corrlist
# resultsdf$xgb_corr_mse=xgbcorr_mselist
# resultsdf$xgb_corr_rsq=xgbcorr_corrlist*xgbcorr_corrlist
