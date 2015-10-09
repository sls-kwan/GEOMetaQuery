library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(gdata)

if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
metadata <- dbConnect(SQLite(),'GEOmetadb.sqlite')
dbListTables(metadata)
dbListFields(metadata, 'gsm')
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
grabgse <- function(convgsm){
  gse<- dbGetQuery(metadata,paste("select gsm.series_id",
                                   " from gsm",
                                   " where gsm.gsm =\"", convgsm,"\"", "\n", sep=""))
  finalgse <- as.vector(gse)
  return(finalgse$series_id)
}
#Grabs Title and Charactersitics for GSM
grabtitle <- function (gsmid){
  titlenchar <- unlist(dbGetQuery(metadata,paste("select gsm.characteristics_ch1, gsm.title, gsm.description, gsm.source_name_ch1",
                                                 " from gsm",
                                                 " where gsm.gsm =\"", gsmid,"\"", "\n", sep="")))
  return(titlenchar)
}



grabgpl <- function (gsmid){
  gplid <- unlist(dbGetQuery(metadata,paste("select gsm.gpl",
                                                 " from gsm",
                                                 " where gsm.gsm =\"", gsmid,"\"", "\n", sep="")))
  return(gplid)
}