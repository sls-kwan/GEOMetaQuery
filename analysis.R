library(ROCR)
library(caret)
library(e1071)    #SVM
###
#  For two class (binary) classification problem
###
source("GEOgrab.R")
source("metamapextract.R")
source("freetext.R")

load("~/COMP4930/Results/GSM/trainingset.RData");
load("~/COMP4930/Results/GSE/metamapscoring.RData");
load("~/COMP4930/Results/GSM/scoringmatrix.RData");

#outmatrix <- cbind(metamapscoring, scoring.matrix, presmetamatrix, presfreematrix)
outmatrix <- cbind(scoring.matrix, metamapscoring)

outmatrix <- scoring.matrix
#outmatrix <- metamapscoring
input <- trainingset 
#Classification between perturbation and case/control
class0.ind<-which(input[,2]==1)
class1.ind<-which(input[,2]!=1)

class0<-outmatrix[class0.ind,]
class1<-outmatrix[class1.ind,]

label<-c(rep(0,length(class0.ind)),rep(1,length(class1.ind)) )
dd<-rbind(class0,class1);
control <- which(trainingset$V2 == 0)
treatment <- which(trainingset$V2 == 1)
modelpredic <- function (z, y){
  ### Feature selection
  flds <- y
  cvflds <- dd[-flds[[z]],]
  #Choose Features that have AUC >0.6
  auc.single<-array(0,length(cvflds[1,]));
  elements <- c(1:length(dd[1,]))
  cvelements <- setdiff(elements, -flds[[z]])
  for(i in cvelements){
    feat<-dd[,i];
    pred<-prediction(feat,label,label.ordering=c(0,1));
    #perf<-performance(pred,"tpr","fpr");
    auc.single[i]<-attr(performance(pred,"auc"),"y.value")[[1]];
    if(auc.single[i]<0.5){
      auc.single[i]<-1-auc.single[i]
    }
  }
  ff.good<-which(auc.single>0.6)
  #training set
  dd.ff<-dd[-flds[[z]], ff.good]
  #test set
  dd.ts <- dd[flds[[z]], ff.good]
  classifier<-svm(dd.ff,label[-flds[[z]]], kernel="radial", cost=1,  probability=TRUE);
  deVal<-attr(predict(classifier,dd.ts, decision.values=TRUE),"decision.values")
  return(list(auc.single, ff.good, deVal))
}

AUCf <- function (x){
  #4 fold cross validation
  nooOFolds <- as.numeric(5)
  #creationofFolds
  
  flds <- createFolds(c(1: nrow(dd)), k=nooOFolds, list=TRUE, returnTrain=FALSE)
  
  enumatc <- c(1:nooOFolds)
  featuresnpred <- lapply(enumatc, modelpredic, flds)
  aucall <- unlist(lapply(enumatc, function(x){return(unlist(featuresnpred[[x]][1]))}))
  features <- unlist(lapply(enumatc, function(x){return(unlist(featuresnpred[[x]][2]))}))
  test <- unlist(lapply(enumatc, function(x){return(unlist(featuresnpred[[x]][3]))}))
  # Get AUROC
  #pred<-prediction(unlist(test), label[as.vector(unlist(flds))],label.ordering=c(0,1));
  #AUC <- attr(performance(pred,"auc"),"y.value")[[1]];
  #perf <- performance(pred, "tpr", "fpr" )
  #return(perf)
  #return(list(AUC, perf, features, aucall))
  return(list(test, label[as.vector(unlist(flds))], features, aucall))
}

getMissCl <- function(fItems){
  p <- cumsum(rle(fItems)$lengths) + 1
  return(p[-length(p)])
}

onetohund <- c(1:20)
AUCnperf <- lapply(onetohund, AUCf)
predictions <- lapply(onetohund, function(x){return(unlist(AUCnperf[[x]][1]))})
labels <- lapply(onetohund, function(x){return(unlist(AUCnperf[[x]][2]))})
#AUCreturn <- unlist(lapply(onetohund, function(x){return(unlist(AUCnperf[[x]][1][[1]]))}))
#perfreturn <- unlist(lapply(onetohund, function(x){return(unlist(AUCnperf[[x]][2][[1]]))}))
featreturn <- unlist(lapply(onetohund, function(x){return(unlist(AUCnperf[[x]][3]))}))
featreturn <- unique(featreturn)
feataucval <- lapply(onetohund, function(x){return(unlist(AUCnperf[[x]][4]))})
#plot(AUCreturn[[55]])
#boxplot(AUCreturn)
pred <- prediction(predictions, labels)
perf <- performance(pred,"tpr","fpr")
AUC <- attr(performance(pred,"auc"),"y.value")
AUC <- unlist(AUC)
mean(AUC)
fPred <- lapply(as.vector(attr(pred, "fp")), getMissCl)

fNeg <- lapply(as.vector(attr(pred, "fn")), getMissCl)

alwaysfP <- gse[as.factor(names(which(table(unlist(fPred)) == 20)))]
alwaysfN <- gse[as.factor(names(which(table(unlist(fNeg)) == 20)))]
par(pty="s")
plot(perf,col="grey82",lty=3, main=paste("ROC Curve with Standard Deviations, AUC =", round(mean(AUC), 3),sep=" "), asp=1) 
plot(perf,lwd=3,avg="vertical",spread.estimate="stddev",add=TRUE)
abline(a=0, b=1, col = "gray60", lty=2)
#### CREATING FEATURE MATRIX HEATMAP#########
colorfile <- paste(getwd(), "colorcoding.csv", sep="/")
colorinput <- read.table(colorfile, sep=",", quote="\"", na.strings = "NA", stringsAsFactors=FALSE)
trainingset$Col <- colorinput$V2[match(trainingset$V3, colorinput$V3)]
resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}
save(featreturn, file="featreturn.RData")
resetPar()
par()
outmatrixfeat <- outmatrix[,featreturn] 
heatmap.2(outmatrix, RowSideColors=as.character(trainingset$Col), margin = c(12,12), key = TRUE, key.xlab=expression('Log'[2]*' Values'), trace="none", cexRow = 0.7, cexCol = 0.4, col = c(colorpanel(length(unique(unlist(as.list(outmatrixfeat)))), "#e5f5f9", "#99d8c9", "#2ca25f")))
heatmap.2(outmatrixfeat, RowSideColors=as.character(trainingset$Col), margin = c(12,12), key = TRUE, trace="none", cexRow = 0.5, cexCol = 0.7, col = c(colorpanel(length(unique(unlist(as.list(outmatrixfeat)))), "#e5f5f9", "#99d8c9", "#2ca25f")))
colorinput[is.na(colorinput$V1),]$V1 <- "NA"
par(mar=c(0, 0, 0, 0))
legend4heat <- legend("bottomleft", legend=c("Perturbation", "Not Perturbation"), cex=1.0, bty="n",col=as.vector(colorinput$V2),pch=19)
legend4heat <- legend("bottomleft", legend=c("Control", "Treatment"), cex=1.0, bty="n",col=as.vector(colorinput$V2),pch=19)
cnames <- colnames(outmatrix)
newold <-cnames[featreturn]
