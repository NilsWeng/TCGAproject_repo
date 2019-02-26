#Firebrowse

rm(list=ls())
library(GenomicRanges)
library(plyr)
#library(nilsTCGA);


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/Firebrowse")


CNV_table <- read.table("CHOL.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"
                        ,header = TRUE,stringsAsFactors = FALSE)


#This command makes duplications on 5 from 8k samples - NO it doesnt?
CNV_table$Sample <- gsub("[A-Z]-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",CNV_table$Sample)

#Remove duplicated rows 
CNV_table <- CNV_table %>% distinct



df2Grange <- function(dfName){
  suppressPackageStartupMessages(library("GenomicRanges"))
  
  GRange_object <- GRanges(seqnames =  dfName$Chromosome,
                           ranges=IRanges(dfName$Start,dfName$End),
                           strand='*',
                           Segment_Mean=dfName$Segment_Mean,
                           Sample=dfName$Sample)
  
  return (GRange_object)
  
  
}


CNV_Grange<-df2Grange(CNV_table)



GeneRegions <- function (genelist) {
  
  
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(BSgenome))
  suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE))
  
  #From Malin, get all genes and the Granges positions stored in 'genelistgranges'
  mart <- "ensembl"
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",GRCh = 37) # OBS !!! What version was used for the annotation? Default GRCh38  
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


hits <- findOverlaps(GenelistGrange,CNV_Grange, type="within")

CNV_Grange_Genes <- CNV_Grange[subjectHits(hits)]
CNV_Grange_Genes@elementMetadata$DeletedGene <- GenelistGrange[queryHits(hits)]@elementMetadata$ensemble_gene_id 
CNV_Grange_Genes@elementMetadata$DeletedGeneHGNC <- GenelistGrange[queryHits(hits)]@elementMetadata$hgnc_symbol


#Remove duplicated rows
DF <- data.frame(CNV_Grange_Genes$Sample,CNV_Grange_Genes$Segment_Mean,CNV_Grange_Genes$DeletedGene,CNV_Grange_Genes$DeletedGeneHGNC)


#Overlapping regions seems to add some duplicated samples 187303 to be exact(CHOL)
DF <- DF %>% distinct
DF <- DF[!DF$CNV_Grange_Genes.DeletedGeneHGNC =="" ,]

colnames(DF) <- c("sample_id","segment_mean","gene_id_ensemble","gene_id_HGNC")

#remove all genes with ensembl id but lacking HGNC symbol
CNV_Grange_Genes<- CNV_Grange_Genes[!CNV_Grange_Genes$DeletedGeneHGNC == "" ,]





# Expression in genes -----------------------------------------------------
EXP_table <- read.table("CHOL.uncv2.mRNAseq_RSEM_normalized_log2.txt",header = TRUE,stringsAsFactors = FALSE)

#remove Entrez gene id from $gene
get_hgcn <- function(name){
  
  name <-  strsplit(name, "\\|")[[1]][1]
  return (name)
}

hgcn_Name <- sapply(EXP_table$gene,get_hgcn);
EXP_table$gene <- hgcn_Name


#Rewrite TCGA-Barcode from TCGA.3X.AAV9.01 -> TCGA-3X-AAV9-01
colnames(EXP_table) <- gsub("\\.","-",colnames(EXP_table))


#Select only genes in changed segments
EXP_hits <- EXP_table[EXP_table$gene %in% DF$gene_id_HGNC ,]


get_expression <- function(gene,id){
  
  
  
  #Select row for right gene
  gene_row <- EXP_hits[EXP_hits$gene == gene ,]
  
  expression <- as.numeric(gene_row[colnames(gene_row) == id])
  
  gene_row <- gene_row[2:length(gene_row)]
  avg_expression  <- as.numeric(rowMeans(gene_row,na.rm = TRUE))
  
  #p <- as.numeric(expression/avg_expression)
  
  return_vector <- c(expression,avg_expression)
  names(return_vector) <- c("Expression","Avg_expression")
  
  #break()
  return(return_vector)
  
}

#gene1 <- "A1BG"
#id1 <- "TCGA-3X-AAVB-01"

#hej <- get_expression(gene1,id1)





#Select on samples present in the mRNA expression data
DF <- DF[DF$sample_id %in% colnames(EXP_hits) ,]
DF <- DF[DF$gene_id_HGNC %in% EXP_hits$gene ,]


DF1 <- sample_n(DF, 50000)
#DF <- DF[1:100 ,]


TEST <- mapply(get_expression,DF1$gene_id_HGNC,as.vector(DF1$sample_id))

DF1$Expression <- as.vector(TEST[1 ,])
DF1$AVG_expression <- as.vector(TEST[2 ,])
DF1$P <- DF1$Expression/DF1$AVG_expression


library(ggplot2)

scatter_plot <- ggplot(DF1, aes(segment_mean, P))
scatter_plot + geom_point() + labs(x = "log2(CN/2)", y = "Expression/AVG expression (within cancer)") + 
  ggtitle("Firebrowse-CHOL")
