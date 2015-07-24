library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(gdata)
library(RCytoscape)
#####
setwd("~/COMP4930")
filenames = c("CHDseries_w_class.csv")
#####
egfiles <- paste(getwd(), filenames[1], sep="/") ##Input has GSE
featfile <- paste(getwd(), "SemanticTypes_2013AA.txt", sep="/") ##Input has GSE

if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')
trainingset <- read.csv(egfiles, header=FALSE)

featset <- read.table(featfile, sep= "|", quote="\"")
output <- list()
gse <- as.vector.factor(trainingset$V1)
trainingset$row.names <- NULL

#Function that grabs all corresponding GSE metadata
graballinfo <- function(gseid) {
  gsemeta <- unlist(dbGetQuery(metadata,paste(GSEQuery,
                             " from gse",
                             " where gse.gse =\"", gseid,"\"", "\n", sep="")))
  #gsetogsm<- lapply(gseid, grabgsm)
  #titlenchargsm <- unlist(lapply(as.list(unlist(gsetogsm)), grabtitle))
  #gseplusgsmmeta <- paste(gsemeta, titlenchargsm, collapse="\n")
  gseplusgsmmeta <- paste(gsemeta, collapse = "\n")
  metmap <- paste ("cd ~ && echo \"", unlist(gseplusgsmmeta), "\" | ./Downloads/public_mm/bin/metamap --XMLf > COMP4930/MetaMapOut.xml")
  system(metmap)
  semtype<- system("./xml.pl", intern = TRUE)
  return(semtype)
}
metadata.all <- lapply(as.list(gse), graballinfo)
metadata.all.metamap <- metadata.all
save(metadata.all.metamap, file="CHD/metamapmetadata.RData")
load("metamapmetadata.RData")

#Fowards each feature to be scored against GSE Metadata
parsetoscore <- function (features, metatable){
  metatable.score <- mapply(elementmeta = metatable, scoring, feature = features)
  return(unlist(metatable.score))
}

#Scoring function for GSE Metadata
scoring <- function(elementmeta, feature){
  if(length(which(grepl(feature, elementmeta, fixed = TRUE))) > 0){
    freq <- length(which(grepl(feature, elementmeta, fixed = TRUE)))
    return(list(as.numeric(freq)))  ##Using Filename
  }
  return(list(as.numeric(0)))
  
}

parsetomat <- function (metatable){
  tfmax <- max(metatable)
  #print(tfmax)
  metatable.score <- lapply(as.list(metatable), FUN=scoringnorm, tfmax)
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

scoring.list <- lapply(as.list(as.character(featset$V3)), parsetoscore, metadata.all.metamap )
adj.scoring.list <- lapply(as.list(scoring.list), parsetomat)
scoring.matrix <- do.call(cbind, adj.scoring.list)


######


rownames(scoring.matrix) <- gse
colnames(scoring.matrix) <- as.vector(featset$V3)
metamapscoring <- scoring.matrix
save(metamapscoring, file ="metamapscoring.RData")
presmetamatrix <- metamapscoring
presmetamatrix[presmetamatrix > 0] <- 1.0

save(presmetamatrix, file="presmetamatrix.RData")
dbDisconnect(metadata)

##TESTS
if(max(scoring.matrix) > 0){
  print("MetaMap is running\n")
} 
if(nrow(trainingset) != nrow(scoring.matrix)){
  print("Scoring matrix is incorrect, wrong GSE")
} else if (length(featset$V3) != ncol(scoring.matrix)){
  print("Incorrect matrix, wrong features")
} else {
  print("Correct matrix")
}

