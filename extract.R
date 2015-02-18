library(GEOquery)
wdin <- paste(getwd(), "input.txt", sep="/") ##Input has GSE
wout <- paste(getwd(), "output.txt", sep="/")
input <- read.table(wdin, quote="\"")
unique(input)  #Make sure there are no duplicates
for (i in input$V1){
  gse <- getGEO(i,GSEMatrix=FALSE, getGPL=F)
  a <- vector()
  for (j in (Meta(gse)$summary)){
    b <- paste(a, j);
    a <- b
  }
  t <- paste(i, " Type", Meta(gse)$type, a)
  write(t, file=wout,append = TRUE)  #appends to file, does not overwrite
}




