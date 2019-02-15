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



#########PROBLEM WITH LOTS OF NA IN SELECTED CNV_TABLE


#Make subset of CNV_table for certain cancertype
CNV_table1 <- CNV_table[CNV_table$Study_abb == "ACC" ,]


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
  suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE))
  
  #From Malin, get all genes and the Granges positions stored in 'genelistgranges'
  mart <- "ensembl"
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") # OBS !!! What version was used for the annotation? Default GRCh38  
  ensemblgenes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name',
                                     'start_position','end_position'), mart = ensembl)
  
  xidx <- which(ensemblgenes[,3]=="X")
  yidx <- which(ensemblgenes[,3]=="Y")
  ensemblgenes[xidx, 2] <- 23
  ensemblgenes[yidx, 2] <- 24
  
  genelistgranges <- GRanges(seqnames=ensemblgenes$chromosome_name,
                             ranges=IRanges(ensemblgenes$start_position, ensemblgenes$end_position),
                             hgnc_symbol=ensemblgenes$hgnc_symbol
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


