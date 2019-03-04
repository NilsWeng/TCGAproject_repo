#Pipeline for finding significant changes in geneexpression
#for samples from highly mutated , per cluster basis.
# TLDr - Underlying mechanisms of mutational patterns

rm(list=ls())


# Library -----------------------------------------------------------------
library(dplyr)



# Read data ---------------------------------------------------------------

#High-mutation samples that have been clustered
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations(3s)")
clustered_samples <- read.table("Sample_In_Cluster_complete_5.txt",header=TRUE,stringsAsFactors = FALSE)
clustered_samples$Sample_id <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",clustered_samples$Sample_id)



#mRNA expression data
setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads")
mRNA_exp <- read.table("PANCAN_mRNA_expression.txt",stringsAsFactors = FALSE,header=TRUE)
colnames(mRNA_exp) <- gsub("\\.","-",colnames(mRNA_exp))
colnames(mRNA_exp) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(mRNA_exp))


#Add row with cluster data to later split by
cluster_vector <- as.data.frame(colnames(mRNA_exp)[-1])
colnames(cluster_vector) <- "Sample_id"
cluster_vector <- merge(cluster_vector,clustered_samples,by="Sample_id",all=TRUE)
#cluster_vector$Cluster[is.na(cluster_vector$Cluster)] <- "0"

row.names(mRNA_exp) <- mRNA_exp$`Hybridization-REF`

mRNA_exp <- select(mRNA_exp, -one_of("Hybridization-REF"))



# Split mRNA into cluster groups ------------------------------------------
mRNA_exp1 <- rbind("Cluster"=cluster_vector$Cluster,mRNA_exp)
mRNA_exp1 <- t(mRNA_exp1)
mRNA_exp1 <- as.data.frame(mRNA_exp1)

#Split table of all high-mut samples in each cluster
mRNA_cluster <- split(mRNA_exp1, f = mRNA_exp1$Cluster)

# Significant - Over all samples/cancers --------------------------------------------------
library(matrixStats)
mRNA_exp <- as.matrix(mRNA_exp)
std_row <- rowMads(mRNA_exp)


