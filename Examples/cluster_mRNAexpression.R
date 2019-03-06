#Pipeline for finding significant changes in geneexpression
#for samples from highly mutated , per cluster basis.
# TLDr - Underlying mechanisms of mutational patterns

rm(list=ls())


# Library -----------------------------------------------------------------
library(dplyr)
library(ggplot2)


# Read data ---------------------------------------------------------------

#High-mutation samples that have been clustered
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations(3s)")
clustered_samples <- read.table("Sample_In_Cluster_complete_7.txt",header=TRUE,stringsAsFactors = FALSE)
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

#rm("mRNA_exp1")

# Significant - Over all samples/cancers --------------------------------------------------
library(matrixStats)
mRNA_exp <- as.matrix(mRNA_exp)
std_row <- rowMads(mRNA_exp)



get_Tstat_for_gene <- function(gene_name){
  #Function that returns p-value of statistic test for each cluster given a gene name.
  
  gene_vector <- vector(mode="numeric", length=length(mRNA_cluster))
  pop_row <- mRNA_exp[rownames(mRNA_exp) == gene_name ,]
  
  for (i in 1:length(mRNA_cluster)){
    
    sub_pop <- mRNA_cluster[[i]]
    sub_pop <- select(sub_pop, -one_of("Cluster"))
    sub_pop <- t(sub_pop)
    
    sub_pop_row <- sub_pop[rownames(sub_pop) == gene_name ,]
    
    #statistic <- t.test(sub_pop_row,pop_row)$p.value
    statistic <- wilcox.test(sub_pop_row,pop_row)$p.value
    
    gene_vector[i] <- statistic
    #qqnorm(sub_pop_row,main="QQ plot of sub_pop",pch=19)
    #qqline(sub_pop_row)
    
    
  }
  
  #Only works when naming 5 clusters
  names(gene_vector) <- paste("cluster",c(1:length(mRNA_cluster)),sep="")
  return(gene_vector)
      
  
}

#Returns p-value of wilcox.test for each cluster and gene
Statistic_DF <- as.data.frame(sapply(rownames(mRNA_exp),get_Tstat_for_gene))


# More advanced statistic (looking within each cancer) --------------------

# Select genes with significant p-value from Statistic_DF -----------------

get_significant_genes <- function(Cluster){
  
  #significant_genes <- colnames(Statistic_DF)[Statistic_DF[i, ] < 0.05]
  
  significant_genes <- colnames(Statistic_DF)[Statistic_DF[(rownames(Statistic_DF) == Cluster) ,] < 0.05]
  return(significant_genes)
  
}

sig_Genes_In_Cluster <- lapply(rownames(Statistic_DF),get_significant_genes)




#Just get T-stat for each cluster given a gene-name
Stat_Gene <- function(gene_name){
  
  pop_row <- mRNA_exp[rownames(mRNA_exp) == gene_name ,]
 
  for (i in 1:length(mRNA_cluster)){
    
    print(paste("Cluster",i,gene_name,sep=" "))
    sub_pop <- mRNA_cluster[[i]]
    sub_pop <- select(sub_pop, -one_of("Cluster"))
    sub_pop <- t(sub_pop)
    
    sub_pop_row <- sub_pop[rownames(sub_pop) == gene_name ,]
    
    #if p<0.05 we cant really assume normal distrubution, ie cant use t-test
    print("Population")
    #Cant take more than 5k samples
    #print(shapiro.test(pop_row))
    #ks.test(pop_row,"pnorm")
    
    print("subpopulation")
    print(shapiro.test(sub_pop_row))
    
    print(t.test(sub_pop_row,pop_row))
    
    #qqnorm(pop_row,main="Normal Q-Q plot",pch=19)
    #print(wilcox.test(sub_pop_row,pop_row))
  }

}

sig_Genes_In_Cluster
Stat_Gene("ITPA")




# Plot distribution of significant genes ----------------------------------

  
plot_distribution <- function(gene_name,cluster){
  
  DF <- mRNA_exp1 %>% select(gene_name,"Cluster")
  DF$Cluster[is.na(DF$Cluster)] <- "0"
  DF$Cluster[DF$Cluster != cluster] <- "Other"
  colnames(DF) <- c("expression","Cluster")
  
  print(ggplot(DF, aes(expression, fill = Cluster)) + 
          geom_density(alpha = 0.2) + xlim(0,6000) +
          ggtitle(paste("Density plot of gene expression",gene_name,sep=" - "))) 
  
  print(qqnorm(DF$expression,
               main=paste("QQ normal plot (All samples)",gene_name,sep=" - "),pch=19))
  print(qqline(as.numeric(DF$expression)))
  
}


sig_Genes_In_Cluster  
genes <- sig_Genes_In_Cluster[[3]]
cluster <- "3"


setwd("C:/Users/Nils_/Downloads/Firebrowse/Pictures")
pdf_name <- paste("Normality_plot_cluster_",cluster,".pdf",sep="")
pdf(pdf_name)

for (gene in genes){
  
  print(plot_distribution(gene,cluster))
  
}

dev.off()

