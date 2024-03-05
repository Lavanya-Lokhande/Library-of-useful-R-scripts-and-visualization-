

setwd("C:/Users/annas/Documents/R") #change to your wd (and make sure you have the input data there)

dat <- read.delim("Example_Data.txt", header = TRUE,sep = "\t") #change file type to suit your format

SampleNames <- as.character(dat[,1]) #change so it fits your data
GeneNames <- colnames(dat[,-c(1:2)]) #change so it fits your data
groups <- dat[,2] #change to the column number where you define your groups

#transpose and define gene data matrix.
genedat <- t(dat[,-c(1:2)]) #change the 2 to the number of annotation columns in your data

#define groups to compare (change to your annotations)
group1<-"Low"
group2<-"High" 

#define output
analysis <- paste(group1," versus ",group2, sep="")
outputfilename <- paste("P-values of ",analysis,".txt", sep="")


subset1 <- is.element(groups , strsplit(group1,",")[[1]])
subset2 <- is.element(groups , strsplit(group2,",")[[1]])
N <- length(GeneNames)
wil<-rep(0,N)
tt<-rep(0,N)

for(k in 1:N){
  t<-t.test(genedat[k,subset1 ],genedat[k,subset2 ])
  w<-wilcox.test(genedat[k,subset1 ],genedat[k,subset2 ])
  wil[k]<-w$p.value
  tt[k]<- t$p.value
}
o<-order(wil)

# Write table to txt-fil
analysis <- paste(group1,"versus",group2)
tabrow<-list("Analyte"=GeneNames[o],"Wilcoxon"=wil[o],"T"=tt[o])
write(analysis ,file=outputfilename )
write("" , file=outputfilename , append=TRUE)
write.table(tabrow,file=outputfilename , sep="\t" , quote=FALSE, row.names=FALSE, append=TRUE)

#output table will end up in your working directory!

