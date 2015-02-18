library(GEOmetadb)
wdin <- paste(getwd(), "input.txt", sep="/") ##Input has GSE
wout <- paste(getwd(), "output.txt", sep="/")
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
input <- read.table(wdin, quote="\"")
unique(input)  #Make sure there are no duplicates
for (i in input$V1){
  rs <- dbGetQuery(con,paste("select gse.gse, gse.summary,",
                           "gse.pubmed_id,gse.type, gse.supplementary_file",
                           " from gse",
                           " where gse.gse =\"", i,"\"", "\n", sep=""))
  write.table(rs, file=wout, row.names = FALSE, col.names = FALSE, sep="\n", append = TRUE) 
}
dbDisconnect(con)
