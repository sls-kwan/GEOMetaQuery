library(ROCR)
library(caret)
library(e1071)    #SVM
###
#  For two class (binary) classification problem
###

source("metamapextract.R")
source("freetext.R")

load("trainingset.RData");
load("metamapscoring.RData");
load("scoringmatrix.RData");

outmatrix <- cbind(metamapscoring, scoring.matrix)


#outmatrix <- scoring.matrix
#outmatrix <- metamapscoring
input <- trainingset 
#Classification between perturbation and case/control
class0.ind<-which(input[,3]==1)
class1.ind<-which(input[,3]!=1)

class0<-outmatrix[class0.ind,]
class1<-outmatrix[class1.ind,]

label<-c(rep(0,length(class0.ind)),rep(1,length(class1.ind)) )
dd<-cbind(t(class0),t(class1));


AUCf <- function (x){
  #4 fold cross validation
  nooOFolds <- as.numeric(4)
  #creationofFolds
  
  flds <- createFolds(c(1: ncol(dd)), k=nooOFolds, list=TRUE, returnTrain=FALSE)
  
  enumatc <- c(1:nooOFolds)
  cvfeatureselection<- function (z){
    ### Feature selection
    cvflds <- dd[,flds[[z]]]
    #Choose Features that have AUC >0.6
    auc.single<-array(0,length(cvflds[,1]));
    elements <- c(1:length(dd[,1]))
    cvelements <- setdiff(elements, flds[[z]])
    for(i in cvelements){
      feat<-dd[i,];
      pred<-prediction(feat,label,label.ordering=c(0,1));
      #perf<-performance(pred,"tpr","fpr");
      auc.single[i]<-attr(performance(pred,"auc"),"y.value")[[1]];
      if(auc.single[i]<0.5){
        auc.single[i]<-1-auc.single[i]
      }
    }
    
    ff.good<-which(auc.single>0.6)
    return(ff.good)
  }
  cvfeatall <- lapply(enumatc, cvfeatureselection)
  cvfeat <- Reduce(intersect, cvfeatall)
  ### Building binary classification
  modelpredic <- function (z){
    #training set
    dd.ff<-dd[cvfeat,-flds[[z]]]
    #test set
    dd.ts <- dd[cvfeat, flds[[z]]]
    classifier<-svm(t(dd.ff),label[-flds[[z]]], kernel="radial", cost=1,  probability=TRUE);
    deVal<-attr(predict(classifier,t(dd.ts), decision.values=TRUE),"decision.values")
    return(deVal)
  }
  test <- unlist(lapply(enumatc, modelpredic))
  # Get AUROC
  pred<-prediction(unlist(test), label[as.vector(unlist(flds))],label.ordering=c(0,1));
  AUC <- attr(performance(pred,"auc"),"y.value")[[1]];
  #perf <- performance(pred, "tpr", "fpr" )
  #return(perf)
  return(cvfeat)
}
onetohund <- c(1:100)
AUCreturn <- lapply(onetohund, AUCf)
#plot(AUCreturn[[55]])
boxplot(unlist(AUCreturn))


#### CREATING FEATURE MATRIX#########
colorfile <- paste(getwd(), "colorcoding.csv", sep="/")
colorinput <- read.table(colorfile, sep=",", quote="\"", na.strings = "NA", stringsAsFactors=FALSE)
trainingset$Col <- colorinput$V2[match(trainingset$V5, colorinput$V3)]

outmatrix <- outmatrix[,featuregood <- Reduce(intersect, AUCreturn)] 

heatmap.2(outmatrix, RowSideColors=as.character(trainingset$Col), trace="none", margins = c(12,12))

colorinput[is.na(colorinput$V1),]$V1 <- "NA"
legend4heat <- legend("topright", legend=colorinput$V1, cex=1.0, bty="n",col=as.vector(colorinput$V2),pch=19)

