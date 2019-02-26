#Correlation betwen gistic and mRNA-exp data

rm(list=ls())

setwd("C:/Users/Nils_/Downloads/Firebrowse/ACC")


#Load used packages
library(data.table)
library(dplyr)

mRNA_exp <- fread("ACC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",stringsAsFactors = FALSE)
GISTIC <- fread("all_data_by_genes.txt",stringsAsFactors = FALSE)


#Sample_id common for both tables
samples_mRNA <- colnames(mRNA_exp)[2:length(colnames(mRNA_exp))]
samples_mRNA <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples_mRNA)

samples_GISTIC <- colnames(GISTIC)[4:length(colnames(GISTIC))]
samples_GISTIC <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples_GISTIC)

common_samples <- intersect(samples_GISTIC,samples_mRNA)


#Shorten TCGA-barcode
colnames(mRNA_exp)[2:length(colnames(mRNA_exp))] <- samples_mRNA
colnames(GISTIC)[4:length(colnames(GISTIC))] <- samples_GISTIC

#Trim away samples not in common_samples
mRNA_exp <- select(mRNA_exp,c("Hybridization REF",common_samples))
GISTIC <- select(GISTIC,c("Gene Symbol","Locus ID","Cytoband",common_samples))



#Gene_id common for both tables
gene_GISTIC <- GISTIC$`Gene Symbol`
#Rename samples with gene_name|21039123 to gene_name
gene_GISTIC <- gsub("\\|.*","",gene_GISTIC)
#Some samples have "small characters" -Problem?



gene_mRNA <- mRNA_exp$`Hybridization REF`
gene_mRNA <- gsub("\\|.*","",gene_mRNA)


common_genes <- intersect(gene_mRNA,gene_GISTIC)

#Rewrite gene-id:s
mRNA_exp$`Hybridization REF` <- gene_mRNA
GISTIC$`Gene Symbol` <- gene_GISTIC


#Trim away genes not in common_genes
GISTIC <- GISTIC[GISTIC$`Gene Symbol` %in% common_genes ,]
mRNA_exp <- mRNA_exp[mRNA_exp$`Hybridization REF` %in% common_genes ,]

#Some genes in GISTIC have same gene_sympol but several locust positions about 300 casesof 18k. Remove these
GISTIC <- GISTIC[!duplicated(GISTIC$`Gene Symbol`) ,] 
mRNA_exp <- mRNA_exp[!duplicated(mRNA_exp$`Hybridization REF`) ,]


#GISTIC cytoband and locust id not needed 
GISTIC <- GISTIC[, -(2:3)]


correlation_DF <- data.frame()


gene_vector <- as.character(GISTIC$`Gene Symbol`)

get_correlation <- function(gene_name){
  
  CN <- as.numeric(GISTIC[GISTIC$`Gene Symbol` == gene_name ,][, -1])
  EXP <- as.numeric(mRNA_exp[mRNA_exp$`Hybridization REF` == gene_name ,][, -1])
  correlation_list <- cor.test(CN,EXP,method="pearson")
  
  return_vector <- c(correlation_list$estimate,correlation_list$p.value)
  names(return_vector) <- c("correlation","p-value")
  return(return_vector)
  
}


#### Returns same value for all !!!
test <- sapply(gene_vector[1:100],get_correlation)



for (gene in GISTIC$`Gene Symbol`){
  
  print(gene)
  CN <- as.numeric(GISTIC[GISTIC$`Gene Symbol` == gene ,][, -1])
  EXP <- as.numeric(mRNA_exp[mRNA_exp$`Hybridization REF` == gene ,][, -1])
  
  
  correlation_list <- cor.test(CN,EXP,method="pearson")
  
  
  correlation_DF <-
  
  
  break()
}
  
