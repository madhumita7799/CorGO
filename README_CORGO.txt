This document is created to help the user to implement CorGO algorithm on their data.

INPUT:
This code requires:
1. Gene expression data (RSEM and log2 normalized) in a tab separated .txt file. File should be in the format given below:

gene	TCGA-2W-A8YY-01	TCGA-4J-AA1J-01	TCGA-BI-A0VR-01	TCGA-BI-A0VS-01	TCGA-BI-A20A-01
ZXDA|7789	5.02443050143841	5.43108803384847	5.81234676007428	5.66225108171523	5.8884315444003
ZXDB|158586	8.22571934014447	8.70361223726215	9.12922221637279	8.01524888895441	8.84652246292395
ZXDC|79364	10.2947127522261	9.85901134565728	11.2683355520681	10.6650096748703	10.665241215006
ZYG11A|440590	-1.48317797655422	3.94279597935065	7.47902436863684	6.59700667921068	5.00452835824537
ZYG11B|79699	9.53274747472069	9.54518960697191	9.62186249908162	8.81241383033939	9.17775488200431

#########################################################################################################################################################################################################

PROCESS:

The algorithm operates in two steps:

1. Run GOSemSim.R to generate semantic similarity matrix from the gene expression data
2. Run CorGO.R with the semantic similarity matrix obtained from Step 1. to generate clusters

#########################################################################################################################################################################################################

OUTPUT:

1. "all_clusters" file: Cluster file with clusters for all the possible values of alpha (with values of alpha, delta and delta2 mentioned in file name)
2. "cluster_validation" file : Cluster Validity Indices (Silhoutte, Dunn and DB Index) file for all the clusters generated at different values of alpha (with values of alpha, delta and delta2 mentioned in file name)
