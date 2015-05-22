library(GEOmetadb)
library(MASS)
library(gplots)
library(RColorBrewer)
library("parallel")
library(gtools)
library(gdata)
wdin <- paste(getwd(), "trainingoutput.csv", sep="/") ##Input has GSE
fdin <- paste(getwd(), "allv3.txt", sep="/") ##Input has GSE
##sdin <- paste(getwd(), "testfeatt.txt", sep="/")
wout <- paste(getwd(), "output.csv", sep="/")
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
input <- read.table(wdin, sep=",", quote="\"")
#class0.ind<-which(input[,5]==2)
#input<-input[class0.ind,]
finput <- read.table(fdin, sep= ",", quote="\"")
output <- list()
gse <- as.vector.factor(input$V1)

#frequency barplot
#freqty <- table(input$V5)
#rownames(freqty) <- c("C", "P", "T", "C+P", "C+T", "P+T")
#barplot(freqty, ylab="Frequencies")
#colnames(outmatrix) <- 
#cout <- colnames(outmatrix) 
#rownames(outmatrix) <- as.vector.factor(input$V1)
#rout <- rownames(outmatrix)
input$row.names <- NULL
SQLQuery <- c("select gse.gse", "gse.summary", "gse.title", "gse.overall_design",
                             "gse.pubmed_id","gse.title", "gse.supplementary_file", "gse.repeats", "gse.repeats_sample_list",
                             "gse.variable", "gse.variable_description")
StringQuery <- paste(as.vector(SQLQuery), collapse=",")
grabgsm <- function(t){
  gsm <- dbGetQuery(con,paste("select gse_gsm.gsm",
                                   " from gse_gsm",
                                   " where gse_gsm.gse =\"", t,"\"", "\n", sep=""))
  finalgsm <- as.vector(gsm)
  return(finalgsm$gsm)
}
grabtitle <- function (z){
  titlenspec <- unlist(dbGetQuery(con,paste("select gsm.characteristics_ch1, gsm.title",
                                " from gsm",
                                " where gsm.gsm =\"", z,"\"", "\n", sep="")))
  #finf <- unlist(dbGetQuery(con, paste("select type from gds_subset where sample_id = \"", z, "\"", "\n", sep="")))
  #titlenspec <- c("")
  #fullgsminfo <- c(titlenspec, finf)
  return(titlenspec)
}
grab <- function(i) {
  j <- i
  rs <- unlist(dbGetQuery(con,paste(StringQuery,
                             " from gse",
                             " where gse.gse =\"", i,"\"", "\n", sep="")))
  conversiongsm <- lapply(j, grabgsm)
  titlegse <- unlist(lapply(as.list(unlist(conversiongsm)), grabtitle))
  modifvec <- paste(rs, titlegse)
  command <- paste ("cd ~ && echo \"", modifvec, "\" | ./Downloads/public_mm/bin/metamap --XMLf > COMP4930/test.xml")
  system(command)
  rsa <- system.time(system("./xml.pl", intern = TRUE))
  rs <- paste(rs, rsa)
  #rs <- c("")
  fullinfo <- list(rs,titlegse)
  ##print(fullinfo)
  return(fullinfo)
}

check <- function (text, sr){
  outmatrix.c <- lapply(sr, clarify, tex = text)
  return(unlist(outmatrix.c))
}
clarify <- function(rss, tex){
  #make case-insensitive
  # Is this needed? Better AUROC without this part
  if(length(which(grepl(tex, rss[[2]]))) > 0){
    return(as.numeric(1*length(which(grepl(tex, rss[[2]])))))  ##Using Filename
  } else if(length(which(grepl(tex,rss[[1]]))) > 0){
    return(as.numeric(0.9*length(which(grepl(tex, rss[[1]])))))  #using other data
  }

  return(as.numeric(0))
  
}
 
outmatrix.grab <- lapply(as.list(gse), grab)


#commatrix <- gtools::permutations(n=nrow(finput), r=2, v=as.vector(finput$V1))
#tempcomp <- paste(commatrix[,1], commatrix[,2], sep=".*")
#outmatrix.c <- lapply(as.list((tempcomp)), check, sr = outmatrix.grab )
outmatrix.clar <- lapply(as.list(as.character(finput$V1)), check, sr = outmatrix.grab )
outmatrix <- do.call(cbind, outmatrix.clar)
#outmatrix <- do.call(cbind , c(outmatrix.c, outmatrix.clar))
rownames(outmatrix) <- gse

#heatmap.2() library(gplots)
#rowsidecolors 


save(outmatrix , file="outmatrix.RData")
colnames(outmatrix) <- as.vector(finput$V1)
#rownames(outmatrix) <- rout
#print(outmatrix)
outmatrixs <- as.matrix(sapply(outmatrix, as.numeric)) 
sumsin <- .rowSums(outmatrix, nrow(outmatrix), ncol(outmatrix))
summ <- matrix(ncol=1, nrow=length(input$V1))
summs <- matrix(ncol=1, nrow=length(input$V1))
summ <- do.call(rbind, as.list(sumsin))
rownames(summ) <- input$V1
colnames(summ) <- "Score"
#print(summ)
ddin <- paste(getwd(), "colorcoding.csv", sep="/")
ddinput <- read.table(ddin, sep=",", quote="\"")
input$Col <- ddinput$V2[match(input$V5, ddinput$V3)]
save(input, file="input.RData")

##save(ddinput, file="ddinput.RData")
#need to cluster them together
#heatmap.2(outmatrix, RowSideColors=as.character(input$Col), trace="none", margins = c(12,12))

outmatrix2 <- outmatrix
outmatrix2[outmatrix2 == 0] <- 0.001
heatmap.2(log2(outmatrix2), RowSideColors=as.character(input$Col), trace="none", margins = c(12,12))
outmatrix <- outmatrix2
save(outmatrix, file="outmatrix.RData")

legend4heat <- legend("topright", legend=ddinput$V1, cex=1.0, bty="n",col=as.vector(ddinput$V2),pch=19)
#save(legend4heat, file="legend.RData")
#write.table(output, file = wout, sep = ",", append = FALSE, eol = "\n", na = "NA", row.names = TRUE, col.names= TRUE);
dbDisconnect(con)
