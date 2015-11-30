library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(textir)
library(gdata)
#####
filenames = c("trainingoutp.csv",
              "feats.txt")
#####
egfiles <- paste(getwd(), filenames[1], sep="/") ##Input has GSE
featfile <- paste(getwd(), filenames[2], sep="/") ##Input has GSE
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')

trainingset <- read.csv(egfiles, header=FALSE)
save(trainingset, file="~/COMP4930/Results/GSE/trainingset.RData")
featset <- read.table(featfile, sep= ",", quote="\"")
output <- list()
gse <- as.vector.factor(trainingset$V1)
gselist <- as.list(gse)
#gsmplatforms <- lapply(GSMList(gselist),function(x) {Meta(x)$platform})
trainingset$row.names <- NULL
#Function that grabs all corresponding GSE metadata
#GSEQuery requires GEOgrab.R
source("GEOgrab.R")
graballinfo <- function(gseid) {
  gsemeta <- unlist(dbGetQuery(metadata,paste(GSEQuery,
                                              " from gse",
                                              " where gse.gse =\"", gseid,"\"", "\n", sep="")))
  gsetogsm<- lapply(gseid, grabgsm)
  titlenchargsm <- unlist(lapply(as.list(unlist(gsetogsm)), grabtitle))
  allinfo <- list(gsemeta,titlenchargsm)
  return(allinfo)
}



metadata.all <- lapply(as.list(gse), graballinfo)
save(metadata.all, file="~/COMP4930/Results/GSE/metadatafreetxt.RData")

#Fowards each feature to be scored against GSE Metadata
parsetoscore <- function (features, metatable){
  metatable.score <- mapply(scoring, elementmeta = metatable, feature = features)
  return(unlist(metatable.score))
}

#Scoring function for GSE Metadata
scoring <- function(elementmeta, feature){
  featCount <- sum(str_count(elementmeta, regex(feature, ignore_case = TRUE)))
  if(featCount> 0){
    freq <- featCount
    return(as.numeric(freq))  ##Using Filename
  }
  
  return(as.numeric(0))
  
}

#Scoring function for GSE Metadata

scoring.list <- lapply(as.list(as.character(featset$V1)), parsetoscore, metatable = metadata.all )
scoring.matrix <- do.call(cbind, scoring.list)
scoring.matrix <- tfidf(scoring.matrix)
scoring.matrix[is.na(scoring.matrix)] <- 0
rownames(scoring.matrix) <- gse

colnames(scoring.matrix) <- as.vector(featset$V1)

save(scoring.matrix , file="~/COMP4930/Results/GSE/scoringmatrix.RData")

dbDisconnect(metadata)

##TESTS
if(max(scoring.matrix) > 0){
  print("Scoring function is running\n")
} 
if(nrow(trainingset) != nrow(scoring.matrix)){
  print("Scoring matrix is incorrect, wrong GSE")
} else if (length(featset$V1) != ncol(scoring.matrix)){
  print("Incorrect matrix, wrong features")
} else {
  print("Correct matrix")
}



