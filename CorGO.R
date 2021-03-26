setwd("<working directory>")

###Input gene expression data file
data1<- read.table("CESC.uncv2.mRNAseq_RSEM_normalized_log2.txt", header=T,row.names=1,sep='\t', check.names = F)

###filter rows with >60% NA
data1<-data1[-which(rowMeans(is.na(data1)) > 0.6), ]
data1[is.na(data1)] <- 0
#genes<- rownames(data1)

##transpose
datat<- t(data1)

##create correlation matrix
rho<-cor(datat,datat)
rownames(rho)<- rownames(data1)
colnames(rho)<- rownames(data1)

###Input semantic similarity matrix file
simmat<- read.table("simmat", header=T, sep=",")

###from entrez id to gene symbol
library(annotate)
genesymbol<- getSYMBOL(rownames(simmat_cc), data='org.Hs.eg')
genesymbol<- as.data.frame(genesymbol)
entrezID<- rownames(genesymbol)
entreztosymbol<- as.data.frame(cbind(entrezID, genesymbol))

####keeping only those genes whose semantic similarity is available
entreztosymbol<- entreztosymbol[!duplicated(entreztosymbol$genesymbol), ]
#entreztosymbol<- entreztosymbol[-278,]
i<- intersect(rownames(data1), entreztosymbol$entrezID)
data1<- data1[i,]
rownames(data1)<- entreztosymbol$genesymbol
i<- intersect(entreztosymbol$entrezID, (rownames(simmat)))
entreztosymbol<- entreztosymbol[i,]
rownames(simmat)<- entreztosymbol$genesymbol
colnames(simmat)<- entreztosymbol$genesymbol
i<- intersect(rownames(rho), rownames(simmat))
rho2<- rho[i,i]
#########13348 genes

##taking absolute value of correlation coefficient
rho2<- abs(rho2)

##set thereshold for delta and delta2
delta<- 0.9
delta2<- 0.7
alpha<- c(0.5, 0.6, 0.7, 0.8, 0.9)

for (k in 1: length(alpha)){
######integrated matrix with combined score from correlation as well as semantic similarity
int_mat<- matrix(c(0), nrow = nrow(rho2), ncol = ncol(rho2))
for (i in 1: nrow(rho2)){
  for (j in 1: ncol(rho2)){
    int_mat[i,j]<- alpha[k]*rho2[i,j]+(1-alpha[k])*simmat_breast2[i,j]
  }
}
rownames(int_mat)<- rownames(rho2)
colnames(int_mat)<- colnames(rho2)

####assigning N score

score<- matrix(data=c(0),nrow= nrow(rho2), ncol= 1)

for ( i in 1: nrow(rho2)){
  for (j in 1: ncol(rho2)){
    if (int_mat[i,j] > delta)
      score[i]= score[i]+1
  }
}

score<- as.data.frame(score)
score_list<- as.data.frame(cbind(rownames(rho2), score$V1))
colnames(score_list)<- c("mRNA", "score")
sc_mat = as.matrix(score_list)

#########sorting in descending order
sort_final = sc_mat[rev(order(as.numeric(sc_mat[,2]))),]
genefinalscore<- sort_final

#######defining threshold for further filtering of genes

avgthr<- mean(as.numeric(genefinalscore[,2]))
tr_filtered<- matrix(nrow = nrow(genefinalscore), ncol = 2)
for ( i in 1: nrow(genefinalscore)){
  if (as.numeric(genefinalscore[i,2])> avgthr)
    tr_filtered[i,]<- genefinalscore[i,]
}

colnames(tr_filtered)<- colnames(genefinalscore)

tr_filtered<- tr_filtered[!(apply(tr_filtered, 1, function(y) any(y == 0))),]
list_score<- tr_filtered[complete.cases(tr_filtered),]


######clustering taking high-scoring genes as initial centroids

int_mat2<- int_mat 
clust_all<- NULL
clust_count<-0
for (i in 1: nrow(list_score)){
   centroid<- list_score[i,1]
for (j in i: nrow(list_score)){ 
  if(int_mat2[which(rownames(int_mat2)== list_score[i,1]), which(colnames(int_mat2)== list_score[j,1])] > delta2){
  centroid <- rbind(centroid, list_score[j,1])
  int_mat2[,which(colnames(int_mat2)== list_score[j,1])]<- 0
  }
}
   
 
  if(length(centroid)>5){
  clust_count<- clust_count+1
  write.csv(centroid, file = paste("Cluster_",i,"-", delta,"-", delta2,"-", alpha[k]))
  centroid<- cbind(centroid, c(clust_count))
  clust_all<- rbind(clust_all, centroid)

  }
   
}

#####print clusters
write.csv(clust_all, file = paste("all_clusters",delta,"-",delta2,"-",alpha[k], clust_count))

############Computing Cluster Validity Indices
########requires: Combined clusters file (generated previously) and gene expression matrix

dist_mat<- as.matrix(dist(data1))

library(clValid)
library( clusterSim)
library(cluster)

########input combined cluster file
u<- read.csv(paste("all_clusters",delta,"-",delta2,"-",alpha[k], clust_count), header = T)
u<- u[-1]
######taking unique genes(removing centroid from first line)
u<- u[!duplicated(u$mRNA), ]

i<- intersect(as.character(u$mRNA), rownames(dist_mat))
clust_data<- dist_mat[i,i]

#dim(clust_data)

clust_data<- as.matrix(clust_data)
#######Silhouette Index
sil<- silhouette(as.numeric(u$X.1), dist = clust_data)
plot(sil, col=c("red","green"))
silhouette_Index<- summary(sil)$avg.width

#######Dunn Index

Dunn_Index<- dunn(clusters= u$X.1, method = "euclidean", Data = clust_data)


######DB Index
db<-index.DB(clust_data, u$X.1, d=NULL, centrotypes="centroids", p=2, q=2)
DB_Index<- db$DB

valid_result<- cbind(silhouette_Index, Dunn_Index, DB_Index)

write.csv(valid_result, paste("cluster_validation",delta,"-",delta2,"-",alpha[k], clust_count))


rm(tr_filtered, sort_final, list_score, int_mat, score, score_list, sc_mat, genefinalscore, avgthr, valid_result)
}

