#Firebrowse

rm(list=ls())
library(GenomicRanges)
library(plyr)
library(dplyr)
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
#EXP_table1 <- read.table("CHOL.uncv2.mRNAseq_RSEM_normalized_log2.txt",header = TRUE,stringsAsFactors = FALSE)





library(readr)
library(data.table)
library(ggplot2)

EXP_table <- fread("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",stringsAsFactors = FALSE)
EXP_table <- as.data.frame(EXP_table)
#EXP_table <- droplevels.data.frame(EXP_table)



#Some modification to make script work for level3 RSEM genes
colnames(EXP_table)<- c("gene",colnames(EXP_table)[2:length(colnames(EXP_table))])
colnames(EXP_table) <- gsub("[A-Z]-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(EXP_table))
colnames(EXP_table) <- as.vector(gsub("-",".",colnames(EXP_table)))
EXP_table <- EXP_table[2:nrow(EXP_table) ,]

#plot Expression data

if (FALSE){
  test <- as.vector(EXP_table[, 2:ncol(EXP_table)])
  test1 <- c()
  for (i in ncol(test)){
    
    test1 <- c(test1,test[, i])
    
  }
  test1 <- as.numeric(test1)
  test1<- data.frame(test1)
  colnames(test1) <- "Expression"
  #Density
  ggplot(test1, aes(Expression)) + geom_density() + xlim(0, 4000)
  #Histogram
  ggplot(test1, aes(Expression)) + geom_histogram() + xlim(0, 4000)
  
  
  
}

#remove Entrez gene id from $gene
get_hgcn <- function(name){
  
  name <-  strsplit(name, "\\|")[[1]][1]
  return (name)
}

hgcn_Name <- sapply(as.character(EXP_table$gene),get_hgcn);
EXP_table$gene <- as.vector(hgcn_Name)


#Rewrite TCGA-Barcode from TCGA.3X.AAV9.01 -> TCGA-3X-AAV9-01
colnames(EXP_table) <- gsub("\\.","-",colnames(EXP_table))


#Select only genes in changed segments
EXP_hits <- EXP_table[EXP_table$gene %in% DF$gene_id_HGNC ,]


get_expression <- function(gene,id){
  
  
 
  #Select row for right gene
  gene_row <-(EXP_hits[EXP_hits$gene == gene ,])
  
  expression <- as.numeric((gene_row[colnames(gene_row) == id]))
  
  gene_row <- gene_row[2:length(gene_row)]
 
  avg_expression  <- as.numeric(mean(as.numeric(gene_row),na.rm = TRUE))
  
  #p <- as.numeric(expression/avg_expression)
  
  return_vector <- c(expression,avg_expression)
  names(return_vector) <- c("Expression","Avg_expression")
  
  #break()
  return(return_vector)
  
}

#gene <- "A1BG"
#id <- "TCGA-3X-AAVB-01"

#hej <- get_expression(gene1,id1)





#Select on samples present in the mRNA expression data
DF <- DF[DF$sample_id %in% colnames(EXP_hits) ,]
DF <- DF[DF$gene_id_HGNC %in% EXP_hits$gene ,]

#DF1 <- DF
DF1 <- sample_n(DF, 100000)
#DF <- DF[1:100 ,]


TEST <- mapply(get_expression,as.vector(DF1$gene_id_HGNC),as.vector(DF1$sample_id))

#TEST<- mapply(get_expression,gene,id)


DF1$Expression <- as.numeric(as.vector(TEST[1 ,]))
DF1$AVG_expression <- as.numeric(as.vector(TEST[2 ,]))
DF1$P <- as.numeric(as.numeric(DF1$Expression)/as.numeric(DF1$AVG_expression))


library(ggplot2)

scatter_plot <- ggplot(DF1, aes(segment_mean, Expression))
scatter_plot + geom_point() + labs(x = "log2(CN/2)", y = "Expression/AVG expression (within cancer)") + 
  ggtitle("Firebrowse-CHOL")



#Density plot of Expression
ggplot(DF1, aes(Expression)) + geom_density() + xlim(0, 4000)


#Find if there is any correlation between CN and EXP for each gene separatley


#clear workspace
rm(list=ls()[! ls() %in% c("DF","DF1")])





genes <- unique(DF1$gene_id_HGNC)

# reduce loop time by only selecting genes with more than 6 samples
Number <- count(DF1,DF1$gene_id_HGNC)
Number <- Number[Number$n > 10 ,]
#DF_subset <- split(DF1,DF1$gene_id_HGNC) Just crashes the entire program
genes <- as.character(genes[genes %in% Number$`DF1$gene_id_HGNC`])

corr_DF <- data.frame()


for (gene in genes){
  
  
  
  DF_subset <- DF1[DF1$gene_id_HGNC == gene ,]
  

  
  N <- nrow(DF_subset)
  correlation <- cor(DF_subset$segment_mean,DF_subset$Expression)
  gene_name <- gene
  
  return_DF <- data.frame(gene_name,N,correlation)
  
  corr_DF <- rbind(corr_DF,return_DF)

  
  
}

corr_DF


#Print CN vs EXP for some

corr_DF <- corr_DF[!is.na(corr_DF$correlation),]
dev.off()

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/Firebrowse/Pictures")
pdf("Correlation_CHOL.pdf")

for (i in 1:nrow(corr_DF)){
  
  subset <- DF1[as.character(DF1$gene_id_HGNC) == as.character(corr_DF$gene_name[i]) ,]
  gene_name <- as.character(unique(subset$gene_id_HGNC))
  
  print(ggplot(subset, aes(segment_mean,Expression)) + geom_point() +
    ggtitle(gene_name) +
    xlab("Segment_mean") + ylab("Expression"))
       
  
  #print(plot(subset$segment_mean,subset$Expression))
  
}

dev.off()
