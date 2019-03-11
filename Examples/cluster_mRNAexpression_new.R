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


sampleDF <- data.frame("Sample_id" = colnames(mRNA_exp)[2:length(colnames(mRNA_exp))])

#Get cancer abb for each sample
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
TSS2Study <- read.table("TSS2Studyabb.txt",header=TRUE,stringsAsFactors = FALSE)
code <- data.frame("TSS.Code"=gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*$","",sampleDF$Sample_id)))
code <- left_join(code,TSS2Study,by="TSS.Code")
sampleDF$cancer <- code$Study.Abbreviation

sampleDF <- left_join(sampleDF,clustered_samples,by="Sample_id")
colnames(sampleDF) <- c("sample","cancer","cluster")
sampleDF[is.na(sampleDF$cluster) ,3] <- "0"

rm(code,TSS2Study)






# More advanced statistic (looking within each cancer) --------------------

#Add cancer_type_data
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cohort <- read.table("cohort.txt",stringsAsFactors = FALSE)


#Loop over each cluster

get_stat <- function(gene_name){
  
  
  gene_exp_tot <- mRNA_exp %>% filter(`Hybridization-REF` == gene_name) %>% select(-one_of("Hybridization-REF"))
  
  return_vector <- c()
  name_vector <- c()
  for (i in 1:length(unique(clustered_samples$Cluster))){
    
    
    samples_in_cluster <- sampleDF %>% filter(cluster==i)
    
    
    for (cancer_type in unique(samples_in_cluster$cancer)){
      
      
      samples_in_cancer <- sampleDF %>% filter(cancer==cancer_type) %>% select(sample)
      samples_in_cancer_cluster <- samples_in_cluster %>% filter(cancer == cancer_type) %>% select(sample)
      #filter away samples in cluster from cancer pop
      samples_in_cancer <- data.frame("sample"=samples_in_cancer[!(samples_in_cancer$sample %in% samples_in_cancer_cluster$sample) ,])
      
      cancer_exp <- as.numeric(gene_exp_tot[, colnames(gene_exp_tot) %in% samples_in_cancer$sample])
      cluster_exp <- as.numeric(gene_exp_tot[, colnames(gene_exp_tot) %in% samples_in_cancer_cluster$sample])
      
      stat <- wilcox.test(cancer_exp,cluster_exp)$p.value
      
      name <- paste(i,cancer_type,sep="_")
      return_vector <- c(return_vector,stat)
      name_vector <- c(name_vector,name)
      
    }
    
  }
  print(paste("Done with",gene_name,sep=" "))
  names(return_vector) <- name_vector
  return(return_vector)
  
}




Statistic_DF <- as.data.frame(sapply(mRNA_exp$`Hybridization-REF`, get_stat))

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
write.table(Statistic_DF,"Statistic_DF.txt")









# Select genes with significant p-value from Statistic_DF -----------------

get_significant_genes <- function(Cluster){
  

  
  significant_genes <- colnames(Statistic_DF)[Statistic_DF[(rownames(Statistic_DF) == Cluster) ,] < 0.05]
  return(significant_genes)
  
}

sig_Genes_In_Cluster <- lapply(rownames(Statistic_DF),get_significant_genes)


get_all_sig_genes <- function(Cluster){
  
  hit <- paste(Cluster,"_",sep="")
  genes <- Statistic_DF[grep(hit,rownames(Statistic_DF)),]
  significant_genes <- colnames(genes[(colMeans(genes) < 0.25)])
  return(significant_genes)
}


all_sig_Genes_In_Cluster <- lapply(c(1:7),get_all_sig_genes)
all_sig_Genes_In_Cluster



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

