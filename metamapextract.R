library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(gdata)
#####
filenames = c("trainingoutput.csv")
#####
egfiles <- paste(getwd(), filenames[1], sep="/") ##Input has GSE
featfile <- paste(getwd(), "SemanticTypes_2013AA.txt", sep="/") ##Input has GSE

if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')
trainingset <- read.table(egfiles, sep=",", quote="\"")

featset <- read.table(featfile, sep= "|", quote="\"")
output <- list()
gse <- as.vector.factor(trainingset$V1)


trainingset$row.names <- NULL
GSESQLQuery <- c("select gse.gse", "gse.summary", "gse.title", "gse.overall_design",
                             "gse.pubmed_id","gse.title", "gse.supplementary_file", "gse.repeats", "gse.repeats_sample_list",
                             "gse.variable", "gse.variable_description")
#Full select Query for GSE
GSEQuery <- paste(as.vector(GSESQLQuery), collapse=",")

#Grabs all GSM ids given GSE
grabgsm <- function(convgse){
  gsm <- dbGetQuery(metadata,paste("select gse_gsm.gsm",
                                   " from gse_gsm",
                                   " where gse_gsm.gse =\"", convgse,"\"", "\n", sep=""))
  finalgsm <- as.vector(gsm)
  return(finalgsm$gsm)
}

#Grabs Title and Charactersitics for GSM
grabtitle <- function (gsmid){
  titlenchar <- unlist(dbGetQuery(metadata,paste("select gsm.characteristics_ch1, gsm.title",
                                " from gsm",
                                " where gsm.gsm =\"", gsmid,"\"", "\n", sep="")))
  return(titlenchar)
}

#Function that grabs all corresponding GSE metadata
graballinfo <- function(gseid) {
  gsemeta <- unlist(dbGetQuery(metadata,paste(GSEQuery,
                             " from gse",
                             " where gse.gse =\"", gseid,"\"", "\n", sep="")))
  gsetogsm<- lapply(gseid, grabgsm)
  titlenchargsm <- unlist(lapply(as.list(unlist(gsetogsm)), grabtitle))
  gseplusgsmmeta <- paste(gsemeta, titlenchargsm)
  metmap <- paste ("cd ~ && echo \"", gseplusgsmmeta, "\" | ./Downloads/public_mm/bin/metamap --XMLf > COMP4930/MetaMapOut.xml")
  system(metmap)
  semtype<- system("./xml.pl", intern = TRUE)
  return(semtype)
}
metadata.all <- lapply(as.list(gse), graballinfo)
metadata.all.metamap <- metadata.all
save(metadata.all.metamap, file="metamapmetadata.RData")

#Fowards each feature to be scored against GSE Metadata
parsetoscore <- function (features, metatable){
  metatable.score <- mapply(scoring, elementmeta = metatable, feature = features, metaindex = seq_along(metatable))
  return(unlist(metatable.score))
}

#Scoring function for GSE Metadata
scoring <- function(elementmeta, feature, metaindex){
  metarow <- trainingset[metaindex,]
  totalrow <- nrow(trainingset)
  nacheck <- as.numeric(is.na(metarow$V5))
  if(nacheck == 1){
    category <- length(which(as.numeric(is.na(trainingset$V5)) == 1))
  } else {
    category <- length(which(trainingset$V5 == metarow$V5))
  }
  if(length(which(grepl(feature, elementmeta, ignore.case = TRUE))) > 0){
    freq <- length(which(grepl(feature, elementmeta, ignore.case = TRUE)))
    adjfreq <- (1+log2(as.numeric(freq)))*log2(as.numeric(totalrow)/as.numeric(category))
    return(as.numeric(adjfreq))  ##Using Filename
  }
  return(as.numeric(0))
  
}


scoring.list <- lapply(as.list(as.character(featset$V3)), parsetoscore, metatable = metadata.all.metamap )
scoring.matrix <- do.call(cbind, scoring.list)

rownames(scoring.matrix) <- gse
colnames(scoring.matrix) <- as.vector(featset$V3)
metamapscoring <- scoring.matrix
save(metamapscoring, file ="metamapscoring.RData")

dbDisconnect(metadata)
