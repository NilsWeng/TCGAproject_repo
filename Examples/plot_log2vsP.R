#Plot log2(cn/2) against expression/avg expression in that cancer
rm(list=ls())


#Needed packages
library(dplyr)
library(plyr)
library(GenomicRanges)

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data");

CNV_table <- read.table('Xena-GDC-PANCAN.masked_cnv.tsv',header = TRUE);


#Add cancer abberation to CNV_table
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3");
TSS2Study <- read.table("TSS2Studyabb.txt",header = TRUE);
cancer_abb_vector <- as.data.frame(gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*$","",CNV_table$Sample)))
colnames(cancer_abb_vector) <- "TSS.Code"
cancer_abb_vector <- join(cancer_abb_vector,TSS2Study,by="TSS.Code")
cancer_abb_vector <- as.vector(cancer_abb_vector[, 2])

CNV_table$Study_abb <- cancer_abb_vector






#Make subset of CNV_table for certain cancertype
selected_cancer <- "ACC"
CNV_table <- CNV_table[CNV_table$Study_abb %in% selected_cancer ,]



#Make a GRange object 
df2Grange <- function(dfName){
  suppressPackageStartupMessages(library("GenomicRanges"))
  
  GRange_object <- GRanges(seqnames =  dfName$Chromosome,
                           ranges=IRanges(dfName$Start,dfName$End),
                           strand='*',
                           Segment_Mean=dfName$Segment_Mean,
                           Sample=dfName$Sample,
                           Study_abb = dfName$Study_abb)
  
  return (GRange_object)
  
  
}


CNV_Grange<-df2Grange(CNV_table)

rm(CNV_table,TSS2Study,cancer_abb_vector)
#Make a postional list over all gene regions
GeneRegions <- function (genelist) {
  
  
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(BSgenome))
  suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE))
  
  #From Malin, get all genes and the Granges positions stored in 'genelistgranges'
  mart <- "ensembl"
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") # OBS !!! What version was used for the annotation? Default GRCh38  
  ensemblgenes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol',"entrezgene", 'chromosome_name',
                                     'start_position','end_position'), mart = ensembl)
  
  xidx <- which(ensemblgenes[,3]=="X")
  yidx <- which(ensemblgenes[,3]=="Y")
  ensemblgenes[xidx, 2] <- 23
  ensemblgenes[yidx, 2] <- 24
  
  genelistgranges <- GRanges(seqnames=ensemblgenes$chromosome_name,
                             ranges=IRanges(ensemblgenes$start_position, ensemblgenes$end_position),
                             hgnc_symbol=ensemblgenes$hgnc_symbol,
                             entrezgene = ensemblgenes$entrezgene
                             ,ensemble_gene_id = ensemblgenes$ensembl_gene_id)
  
  
  
  keeplist <- c(1:24)
  tokeep<- keeplist[which(keeplist %in% levels(factor(seqnames(genelistgranges))))]
  genelistgranges<- keepSeqlevels(genelistgranges,tokeep, pruning.mode="coarse")
  
}
GenelistGrange <- GeneRegions()




#Find all genes within regions of CNVs segements

hits <- findOverlaps(GenelistGrange,CNV_Grange, type="within")  


CNV_Grange_Genes <- CNV_Grange[subjectHits(hits)]
CNV_Grange_Genes@elementMetadata$DeletedGene <- GenelistGrange[queryHits(hits)]@elementMetadata$ensemble_gene_id 
CNV_Grange_Genes@elementMetadata$DeletedGeneHGNC <- GenelistGrange[queryHits(hits)]@elementMetadata$hgnc_symbol
CNV_Grange_Genes@elementMetadata$DeletedGeneEntrez <- GenelistGrange[queryHits(hits)]@elementMetadata$entrezgene


#A <- as.vector(CNV_Grange_Genes@elementMetadata$DeletedGene)
# <- as.vector(CNV_Grange_Genes@elementMetadata$DeletedGeneHGNC)
#C <- as.vector(CNV_Grange_Genes@elementMetadata$DeletedGeneEntrez)

#test <- data.frame(A,B,C)

# Gene expression ---------------------------------------------------------


library(readr)

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")

EXP_table <- as.data.frame(read_tsv('mRNA.geneEXP.tsv'))


get_hgcn <- function(name){
  
  name <-  strsplit(name, "\\|")[[1]][1]
  return (name)
}

hgcn_Name <- sapply(EXP_table$gene_id,get_hgcn);
unname(hgcn_Name);
EXP_table$gene_id <- hgcn_Name;

#Shorter TCGA-barcode
colnames(EXP_table) <- c("gene_id",gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(EXP_table)[2:length(colnames(EXP_table))]))


#Select only genes in changed segments
EXP_hits <- EXP_table[EXP_table$gene_id %in% CNV_Grange_Genes$DeletedGeneHGNC ,];

#Add rownames as gene_id and remove gene_id as column
row.names(EXP_hits) <- EXP_hits$gene_id
EXP_hits <- EXP_hits[, !colnames(EXP_hits) == "gene_id"]


#translate TSS into cancer type
TSS_vector <- colnames(EXP_hits)
TSS_vector <- gsub("-[A-Z0-9]*-[A-Z0-9]*$","",TSS_vector)
TSS_vector <- gsub("TCGA-","",TSS_vector)
TSS_vector <- as.data.frame(TSS_vector)
colnames(TSS_vector) <- "TSS.Code"

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3");
TSS2Study <- read.table("TSS2Studyabb.txt",header = TRUE)
#Create vector of studytype matching EXP_hits
library(dplyr)
library(plyr)



TSS_vector <- join(TSS_vector,TSS2Study, by = "TSS.Code")
TSS_vector <- as.vector(TSS_vector[,2])

#Select only samples from selected_cancer type.
EXP_hits <- EXP_hits[TSS_vector %in% selected_cancer] 


#EXP hits now only contain relevant rows and colums





#Get gene expression given a sample_id and Gene_id
get_gene_expression <- function(Gene_id,Sample_id){
  
  
  #Sample_id ="TCGA-OR-A5J1-01A"
  #Gene_id = "A1BG"
  

  selected_rows <- EXP_hits[(rownames(EXP_hits) %in% Gene_id) ,]
  
  #Get sample expression
  expression <- as.numeric(selected_rows[, Sample_id])
  
  #Get avg expression for that cancertype
  
  
  avg_expression_cancer <- as.numeric(rowMeans(selected_rows,na.rm = TRUE))
  
  #avg_expression_cancer <- as.numeric(AVG_expression_cancer[row.names(AVG_expression_cancer) == Gene_id ,])
  
  
  P <- as.numeric(expression/avg_expression_cancer)
  return_vector <- c(P)
  names(return_vector) <- c("P")
  return(return_vector)
  
}




sample_id <- as.vector(CNV_Grange_Genes$Sample);
gene_id <- CNV_Grange_Genes$DeletedGeneHGNC;
segment_mean <- CNV_Grange_Genes$Segment_Mean

DF <- data.frame(sample_id,gene_id,segment_mean);

#Since matching on hgnc lots of rows are missing hgnc symbol(but have ensembl id)
DF <-  DF[!DF$gene_id == "" , ]

#Crashes because not all samples are present in both. Checked two samples , these didnt have gene expression data (from GBM)
DF <- DF[DF$sample_id %in% colnames(EXP_hits) ,];
DF <- DF[DF$gene_id %in% rownames(EXP_hits) ,];

DF <- na.omit(DF)

#get gene expression for these samples.
#DF_B <- DF
#Test reducing DF size

DF <- DF[ (DF$segment_mean < -0) ,]




#Make a random sampling from all genes
DF <- sample_n(DF, 50000)
#DF <- DF[1:100 ,]




#result_Matrix <- as.data.frame(mapply(get_gene_expression,DF$gene_id,as.vector(DF$sample_id)))
DF$P <- mapply(get_gene_expression,DF$gene_id,as.vector(DF$sample_id))
#NAN values when there is no expression in any of the samples (since 0/0)


#Save DF
#setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/log2vsP")
#save(DF,)
#DF$CN <- (2^DF$segment_mean)* 2

#Plot log2 vs P
plot(-DF$segment_mean,DF$P,
     xlab = "-log2(CN/2)",
     ylab = "Expression/Avg Expression")


#Correlation
cor.test(DF$segment_mean,DF$P)

#colnames(result_Matrix) <- paste(as.vector(DF$sample_id),DF$gene_id, sep= "<-->")
#remove duplicated results
#result_Matrix <- result_Matrix[, !duplicated(colnames(result_Matrix))]
#result_Matrix

library(ggplot2)

scatter_plot <- ggplot(DF, aes(segment_mean, P))
scatter_plot + geom_point() + labs(x = "- log2(CN/2)", y = "Expression/AVG expression (within cancer)") + 
  ggtitle("ACC 50k genes")
  #+  geom_smooth(method="auto")


scatter_plot <- ggplot(DF, aes(segment_mean, log(P)))
scatter_plot + geom_point() + labs(x = "- log2(CN/2)", y = "log(Expression/AVG expression (within cancer))") 
