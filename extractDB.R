library(GEOmetadb)
library(MASS)
wdin <- paste(getwd(), "input.txt", sep="/") ##Input has GSE
fdin <- paste(getwd(), "testfeat.txt", sep="/") ##Input has GSE
wout <- paste(getwd(), "output.csv", sep="/")
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
input <- read.table(wdin, quote="\"")
finput <- read.table(fdin, quote="\"")
output <- list()
gse <- as.vector.factor(input$V1)
outmatrix <- matrix(ncol=length(finput$V1), nrow=length(input$V1))
rownames(outmatrix) <- gse
rout <- rownames(outmatrix)
colnames(outmatrix) <- as.vector.factor(finput$V1)
cout <- colnames(outmatrix) 
unique(input)  #Make sure there are no duplicates
grab <- function(i) {
  rs <- dbGetQuery(con,paste("select gse.gse, gse.summary, gse.title, gse.overall_design,",
                             "gse.pubmed_id,gse.title, gse.supplementary_file, gse.repeats, gse.repeats_sample_list,",
                             "gse.variable, gse.variable_description",
                             " from gse",
                             " where gse.gse =\"", i,"\"", "\n", sep=""))
  return(rs)
}
check <- function (sr, text){
  outmatrix.c <- lapply(sr, clarify, tex = text)
}
clarify <- function(rss, tex){
  if(length(which(grepl(tex, rss))) > 0){
    return(as.numeric(1))
  }
    return(as.numeric(0))
  
}

  outmatrix.grab <- lapply(as.list(rownames(outmatrix)), grab)
  outmatrix.clar <- lapply(as.list(colnames(outmatrix)), check, sr = outmatrix.grab )
  outmatrix <-do.call(cbind, outmatrix.clar)
  colnames(outmatrix) <- cout
  rownames(outmatrix) <- rout
  print(outmatrix)
#write.table(output, file = wout, sep = ",", append = FALSE, eol = "\n", na = "NA", row.names = TRUE, col.names= TRUE);
dbDisconnect(con)
