library(GEOmetadb)
library(MASS)
wdin <- paste(getwd(), "input.txt", sep="/") ##Input has GSE
wout <- paste(getwd(), "output.csv", sep="/")
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
input <- read.table(wdin, quote="\"")
output <- list()
outmatrix <- matrix()
unique(input)  #Make sure there are no duplicates
for (i in input$V1){
  rs <- dbGetQuery(con,paste("select gse.gse, gse.summary, gse.title, gse.overall_design,",
                           "gse.pubmed_id,gse.title, gse.supplementary_file, gse.repeats, gse.repeats_sample_list,",
                           "gse.variable, gse.variable_description",
                           " from gse",
                           " where gse.gse =\"", i,"\"", "\n", sep=""))
  ##output <- unlist(rs, use.names = FALSE)
  ##print(output)
  output <- rbind(output, rs)
  ##outmatrix <- cbind(outmatrix, output)
  ##print(outmatrix)
 
}
write.table(output, file = wout, sep = ",", append = FALSE, eol = "\n", na = "NA", row.names = FALSE, col.names= FALSE);
dbDisconnect(con)
