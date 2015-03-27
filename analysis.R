
###
#  For two class (binary) classification problem
###
load("outmatrix.RData");
load("input.RData");

class0.ind<-which(input[,5]==1)
class1.ind<-which(input[,5]==2)

class0<-outmatrix[class0.ind,]
class1<-outmatrix[class1.ind,]

label<-c(rep(0,length(class0.ind)),rep(1,length(class1.ind)) )
dd<-cbind(t(class0),t(class1));

heatmap(dd,scale="none")
### Cross-validation
##Split up 10% test, 90% training. recursively 10 times. 
#predict after
library(caret)
nooOFolds <- as.numeric(4)
flds <- createFolds(c(1: ncol(dd)), k=nooOFolds, list=TRUE, returnTrain=FALSE)



### Feature selection
library(ROCR)

auc.single<-array(0,length(dd[,1]));
for(i in 1:length(dd[,1])){
   feat<-dd[i,];
   pred<-prediction(feat,label,label.ordering=c(0,1));
   #perf<-performance(pred,"tpr","fpr");
   auc.single[i]<-attr(performance(pred,"auc"),"y.value")[[1]];
   if(auc.single[i]<0.5){
       auc.single[i]<-1-auc.single[i]
   }
}

ff.good<-which(auc.single>0.6)

pdf(file="AUROC_single_C_P.pdf", width=8, height=7)
par(mai=c(1.5, 0.9, 1,0.1)); #c(bottom, left, top, right)
barplot(sort(auc.single,decreasing=T)-0.5,las=2, main="AUROC of single feature", horiz=F,offset=0.5,ylim=c(0.5,1));
abline(h=0.9,col="red",lty=2,lwd=1.5);
dev.off();

round(sort(auc.single,decreasing=T),3)




### Building binary classification
library(e1071)    #SVM
enumatc <- c(1:nooOFolds)
modelpredic <- function (z){
  dd.ff<-dd[ff.good,-flds[[z]]]
  dd.ts <- dd[ff.good, flds[[z]]]
  classifier<-svm(t(dd.ff),label[-flds[[z]]], kernel="radial", cost=1,  probability=TRUE);
  
  deVal<-attr(predict(classifier,t(dd.ts), decision.values=TRUE),"decision.values")
  barplot(t(deVal),las=2)
  abline(h=0,col="red")
  return(deVal)
}
test <- lapply(enumatc, modelpredic)
## Get AUROC
pred<-prediction(unlist(test), label[as.vector(unlist(flds))],label.ordering=c(0,1));
attr(performance(pred,"auc"),"y.value")[[1]];
#[1] 0.7746711        ##cross validation

##randomForest implementation


