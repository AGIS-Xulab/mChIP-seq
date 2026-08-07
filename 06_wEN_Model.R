################################################################################
################################################################################
#####
###set up directory
setwd("path to/05-combine_bin")

required.packages <- c("tidyverse",
                       "data.table",
                       "glmnet",
                       "caret",
                       "limma",
					   "pROC")
installed.packages <- rownames(installed.packages())
if (all(required.packages %in% installed.packages)) {
  for (package in required.packages) {
    suppressMessages(library(package, character.only=TRUE))
  }
} else {
  required.packages[!required.packages %in% installed.packages]
  stop(paste("Install these packages: ", paste(required.packages[!required.packages %in% installed.packages], collapse=",")))
}


#generate variables

nfolds <- 5
k <- 2
q <- 0.6
balanced <- TRUE
set.seed(12345)

#define function for splitting samples into random folds
balanced.folds <- function(y, nfolds=min(min(table(y)), 10)) {
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)     
  nfolds <- max(nfolds, 2) # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)  # nice we to get the ids in a list, split by class
  
  # make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds)*nfolds, length(totals))
  for (i in seq(totals)) {
    if (length(yids[[i]])>1) { bigmat[seq(totals[i]), i] <- sample(yids[[i]]) }
    if (length(yids[[i]])==1) { bigmat[seq(totals[i]), i] <- yids[[i]] }
  }
  # reshape the matrix
  smallmat <- matrix(bigmat, nrow=nfolds)
  
  #now do a clever sort to mix up the NAs
  #smallmat <- permute.rows(t(smallmat))
  smallmat <- t(smallmat)
  smallmat <- smallmat[sample(nrow(smallmat)),]
  
  #now a clever unlisting
  #note the "clever" unlist doesn't work when there are no NAs
  apply(smallmat, 2, function(x) x[!is.na(x)])
  res <- vector("list", nfolds)
  for(j in 1:nfolds) {
    jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j]
  }
  return(res)
}

if (balanced) {
  folds=balanced.folds(ytrain, nfolds)
  foldss=rep(NA, length(ytrain))
  for (i in 1:nfolds) { foldss[folds[[i]]]=i }
  folds=foldss
} else {
  folds=sample(1:nfolds, size=length(ytrain), replace=TRUE)
}

#create output dir
outdir <- file.path(dirname(fr_file),"caret")
dir.create(outdir, showWarnings=F, recursive=T)
###file to get frag info
fr_file <- "binfrag_summary.rds"
sampleInfo<-"Path to sampleInfo.txt"

#read in data
fr <- readRDS(fr_file)
meta <- read.table(sampleInfo,header=TRUE)
meta <-meta |> rename(id=Sample)

fr.ratio <- fr %>% select("id","bin","ratio1","ratio1.corrected") %>% inner_join(meta, by="id") 

####reshape to wider
fr.wider <- fr.ratio %>% 
            mutate( cancer = recode(Type,
                         `Tum` = "Cancer",                    
                         `Hea` = "No cancer")
			) %>% 
			pivot_wider(names_from = bin,values_from = c(ratio1,ratio1.corrected))
			
colnames(fr.wider)[1]<-"sample"


##prepared df for tranning and validation
fr.df<-fr.wider %>% select(sample, cancer, ratio1.corrected_1:ratio1.corrected_492)


#split samples into train and test sets
inTrain <- createDataPartition( y = fr.wider$cancer,
  ## the outcome data are needed
  p = .5,
  ## The percentage of data in the
  ## training set
  list = FALSE
)

training=fr.df[inTrain,]
xtrain<-training %>% select(-sample,-cancer)
ytrain<-training$cancer
testing=fr.df[-inTrain,]
xval=testing %>% select(-sample,-cancer)
yval=testing$cancer

#center and scale all features
xtrainpp <- preProcess(xtrain, method=c("center", "scale"))
xtrain <- predict(xtrainpp, newdata=xtrain)

############
#####TrainEN
#establish hyperparameter tuning grid
#alphas <- 1
alphas <- seq(0, 1, length=11)
len <- 10
init <- glmnet(xtrain, ytrain, family="binomial", alpha=0.5, nlambda=len)
lambdas <- unique(init$lambda)
lambdas <- lambdas[-c(1, length(lambdas))]
lambdas <- lambdas[1:min(length(lambdas), len)]
fitGrid <- expand.grid(alpha=alphas, lambda=lambdas) # creates grid of desired alphas and lambdas
fitControl <- trainControl("cv", number=5, classProbs=TRUE, summaryFunction=twoClassSummary) # runs 5-fold CV for every combo of alpha/lambda in tuneGrid


#train model 
sample_pred <- list()
feature_sel <- list()
metrics <- list()
rocobj <- list()
if (nfolds > 1) {
  for (i in 1:nfolds) {
    print(paste0("Fold ", i))
    
    xout <- xtrain[folds!=i,,drop=FALSE]
    yout <- ytrain[folds!=i]
    xin <- xtrain[folds==i,,drop=FALSE]
    yin <- ytrain[folds==i]
    
    model <- train(xout, make.names(yout), method="glmnet", family="binomial", metric="ROC", 
                   trControl=fitControl, tuneGrid=fitGrid)

    # select best model (i.e. highest accuracy)
    alpha <- model$bestTune$alpha
    lambda <- model$bestTune$lambda
    sample_pred[[i]] <- data.frame(fold=i,
                                     ID=rownames(xin),
                                     true_class=yin,
                                     response=predict(model, xin, type="prob", s=lambda)[,2]) 
  }
}



# generate final model using FULL training set
# above used to estimate the performance of this model
print("All folds")

model <- train(xtrain, make.names(ytrain), method="glmnet", family="binomial", metric="ROC", 
               trControl=fitControl, tuneGrid=fitGrid)

# select best model (i.e. highest accuracy)
alpha <- model$bestTune$alpha
lambda <- model$bestTune$lambda


# extract info for samples & assess performance
sample_pred[[length(sample_pred)+1]] <- data.frame(fold="full",
                                                   ID=rownames(xtrain),
                                                   true_class=ytrain,
                                                   response=predict(model, xtrain, type="prob", s=lambda)[,2]) 






#  write out model performance metrics
sample_pred <- do.call(rbind, sample_pred)
write.table(sample_pred %>% filter(fold=="full"), file=file.path(outdir, "model_sample_predictions.txt"), quote=F, row.names=F, sep="\t")

trainroc <- roc(sample_pred$true_class[sample_pred$fold=="full"], sample_pred$response[sample_pred$fold=="full"])
cvroc <- roc(sample_pred$true_class[sample_pred$fold!="full"], sample_pred$response[sample_pred$fold!="full"])

metrics <- data.frame(train=c("Train", "CV"), auc=c(trainroc$auc, cvroc$auc))
write.table(metrics, file=file.path(outdir, "model_metrics.txt"), quote=F, row.names=F, sep="\t")


# write out final model
saveRDS(model, file=file.path(outdir, "H3K9me3_model.rds"))

####################
##########validateEN
# pre-process val data
pp<-xtrainpp
xval <- predict(pp, newdata=xval)


#read in final model
###model <- readRDS(mdl_file)


# make predictions
sample_pred <- data.frame(ID=rownames(xval), 
                          true_class=yval,
                          response=predict(model, newdata=xval, type="prob", s="lambda.min")[,2])
write.table(sample_pred, file.path(outdir, "model_sample_predictions_validation.txt"), row.names=FALSE, quote=FALSE, sep="\t")


val_roc <- roc(sample_pred$true_class, sample_pred$response) ## 0.9625
aucs <- data.frame(train="Validation", 
                   auc=val_roc$auc)

pdf("validation_roc.pdf")
plot(val_roc,type="S",col="red")
dev.off()
