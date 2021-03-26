# CorGO
An Integrated Method for Clustering Functionally Similar Genes

# Related article
Pant, N., Madhumita, M. & Paul, S. CorGO: An Integrated Method for Clustering Functionally Similar Genes. 
Interdiscip Sciences Computational Life Sciences (2021). https://doi.org/10.1007/s12539-021-00424-9.

# Abstract
Identification of groups of co-expressed or co-regulated genes is critical for exploring the underlying mechanism
behind a particular disease like cancer. Condition-specific (disease-specific) gene-expression profiles acquired
from different platforms are widely utilized by researchers to get insight into the regulatory mechanism of the disease. 
Several clustering algorithms are developed using gene expression profiles to identify the group of similar genes.
These algorithms are computationally efficient but are not able to capture the functional similarity present between
the genes, which is very important from a biological perspective. In this study, an algorithm named CorGO is introduced,
that specifically deals with the identification of functionally similar gene-clusters. Two types of relationships are
calculated for this purpose. Firstly, the Correlation (Cor) between the genes are captured from the gene-expression data,
which helps in deciphering the relationship between genes based on its expression across several diseased samples.
Secondly, Gene Ontology (GO)-based semantic similarity information available for the genes is utilized, that helps
in adding up biological relevance to the identified gene-clusters. A similarity measure is defined by integrating
these two components that help in the identification of homogeneous and functionally similar groups of genes. 
CorGO is applied to four different types of gene expression profiles of different types of cancer. Gene-clusters
identified by CorGO, are further validated by pathway enrichment, disease enrichment, and network analysis. 
These biological analyses demonstrated significant connectivity and functional relatedness within the 
genes of the same cluster. A comparative study with commonly used clustering algorithms is also performed
to show the efficacy of the proposed method.
