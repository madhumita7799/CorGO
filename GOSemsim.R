######for obtaining semantic similarity matrix from gene list
setwd("current working directory")

####input log2 RSEM normalized gene expression data
rna <- read.table('CESC.uncv2.mRNAseq_RSEM_normalized_log2.txt', header=T,row.names=1,sep='\t', check.names = F)

###filtering rows with >60% NA
rna<-rna[-which(rowMeans(is.na(rna)) > 0.6), ]
rna[is.na(rna)] <- 0
gene<- rownames(rna)
library(org.Hs.eg.db)
library(GOSemSim)

hs <- org.Hs.eg.db
my.symbols <- as.character(genes)
a<- select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")

###Semantic Similarity on the basis of involvement in the same biological process is computed
d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
simmat<- mgeneSim(a$ENTREZID, semData = d, measure = "Wang", drop = "IEA", combine = "BMA")

###Output semantic similarity file
write.csv(simmat, file="simmat")




