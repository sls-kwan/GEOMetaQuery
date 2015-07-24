library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(gdata)
#####
filenames = c("trainingoutput.csv",
              "features.txt")
#####
egfiles <- paste(getwd(), filenames[1], sep="/") ##Input has GSE
featfile <- paste(getwd(), filenames[2], sep="/") ##Input has GSE
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')

trainingset <- read.csv(egfiles, header=FALSE)
save(trainingset, file="trainingset.RData")
featset <- read.table(featfile, sep= ",", quote="\"")
output <- list()
gse <- as.vector.factor(trainingset$V1)
gselist <- as.list(gse)
#gsmplatforms <- lapply(GSMList(gselist),function(x) {Meta(x)$platform})
trainingset$row.names <- NULL
#Function that grabs all corresponding GSE metadata
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
save(metadata.all, file="metadatafreetxt.RData")

#Fowards each feature to be scored against GSE Metadata
parsetoscore <- function (features, metatable){
  metatable.score <- mapply(scoring, elementmeta = metatable, feature = features, metaindex = seq_along(metatable))
  return(unlist(metatable.score))
}

#Scoring function for GSE Metadata
scoring <- function(elementmeta, feature, metaindex){
  if(length(which(grepl(feature, elementmeta[[2]], ignore.case = TRUE))) > 0){
    freq <- length(which(grepl(feature, elementmeta[[2]], ignore.case = TRUE)))
    #adjfreq <- (1+log2(as.numeric(freq)))*log2(as.numeric(totalrow)/as.numeric(category))
    return(as.numeric(freq))  ##Using Filename
  } else if(length(which(grepl(feature, elementmeta[[1]], ignore.case = TRUE))) > 0){
    freq <- length(which(grepl(feature, elementmeta[[1]], ignore.case = TRUE)))
    #adjfreq <- (1+log2(as.numeric(freq)))*log2(as.numeric(totalrow)/as.numeric(category))
    return(as.numeric(freq))
  }
  
  return(as.numeric(0))
  
}
parsetomat <- function (metatable){
  tfmax <- max(metatable)
  #print(tfmax)
  metatable.score <- lapply(metatable, FUN=scoringnorm, tfmax = tfmax)
  return(unlist(metatable.score))
  #return(unlist(metatable.score))
  #return(tfmax)
}

#Scoring function for GSE Metadata
scoringnorm <- function(elementmeta, tfmax){
  tfmax <- tfmax
  if(as.numeric(elementmeta) > 0){
    freq <- elementmeta
    adjfreq <- 0.4+(1-0.4)*(freq/tfmax)
    return(as.numeric(adjfreq))  ##Using Filename
  }
  return(as.numeric(0))
  
}

scoring.list <- lapply(as.list(as.character(featset$V1)), parsetoscore, metatable = metadata.all )
adj.scoring.list <- lapply(scoring.list, parsetomat)
scoring.matrix <- do.call(cbind, adj.scoring.list)
rownames(scoring.matrix) <- gse

colnames(scoring.matrix) <- as.vector(featset$V1)

save(scoring.matrix , file="scoringmatrix.RData")

presfreematrix <- scoring.matrix
presfreematrix[presfreematrix > 0] <- 1.0
save(presfreematrix, file="presfreematrix.RData")
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



