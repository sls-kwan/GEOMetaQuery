wdin <- paste(getwd(), "input.txt", sep="/")
wout <- paste(getwd(), "output.txt", sep="/")
input <- read.table(wdin, quote="\"")
unique(input)
for (i in input$V1){
  gse <- getGEO(i,GSEMatrix=FALSE, getGPL=F)
  a <- vector()
  for (j in (Meta(gse)$summary)){
    b <- paste(a, j);
    a <- b
  }
  t <- paste(i, " Type", Meta(gse)$type, a)
  write(t, file=wout,append = TRUE)
}




