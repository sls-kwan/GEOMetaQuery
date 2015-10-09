library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(textir)
library(gdata)
#####
filenames = c("sample_types.csv",
              "features_samples.txt")
#####
egfiles <- paste(getwd(), filenames[1], sep="/") ##Input has GSE
featfile <- paste(getwd(), filenames[2], sep="/") ##Input has GSE
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
  metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')

GEOmetadata <- read.csv("~/COMP4930/GEOmetaScrape/GEOmetadata.txt", header=FALSE)
GEOmetadata$V9 <- NULL
trainingset <- read.csv(egfiles, header=FALSE)
save(trainingset, file="Results/GSM/trainingset.RData")
featset <- read.table(featfile, sep= ",", quote="\"")
output <- list()
gsm <- as.vector.factor(trainingset$V1)
gsmlist <- as.list(gsm)
rownames(GEOmetadata) <- gsm
#gsmplatforms <- lapply(GSMList(gselist),function(x) {Meta(x)$platform})
trainingset$row.names <- NULL
#Function that grabs all corresponding GSE metadata
graballinfo <- function(gsmid) {
  #gsemeta <- unlist(dbGetQuery(metadata,paste(GSEQuery,
  #                                            " from gse",
  #                                            " where gse.gse =\"", gseid,"\"", "\n", sep="")))
  #gsetogsm<- lapply(gseid, grabgsm)
  GSMrow <- which(rownames(GEOmetadata) == gsmid)
  titlenchargsm <- GEOmetadata[GSMrow,]
  titlenchargsm <- do.call(paste, titlenchargsm)
  return(titlenchargsm)
}



metadata.all <- lapply(as.list(gsm), graballinfo)
save(metadata.all, file="Results/GSM/metadatafreetxt.RData")

#Fowards each feature to be scored against GSE Metadata
parsetoscore <- function (features, metatable){
  metatable.score <- mapply(scoring, elementmeta = metatable, feature = features, metaindex = seq_along(metatable))
  return(unlist(metatable.score))
}

#Scoring function for GSE Metadata
scoring <- function(elementmeta, feature, metaindex){
  featCount <- sum(str_count(elementmeta, regex(feature, ignore_case = TRUE)))
  if(featCount> 0){
    freq <- featCount
    return(as.numeric(freq))
  }
  
  return(as.numeric(0))
  
}

scoring.list <- lapply(as.list(as.character(featset$V1)), parsetoscore, metatable = metadata.all )
scoring.matrix <- do.call(cbind, scoring.list)
scoring.matrix <- tfidf(scoring.matrix)
scoring.matrix[is.na(scoring.matrix)] <- 0
rownames(scoring.matrix) <- gsm

colnames(scoring.matrix) <- as.vector(featset$V1)

save(scoring.matrix , file="Results/GSM/scoringmatrix.RData")

#save occurence matrix
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



