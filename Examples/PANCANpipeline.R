#Test of handling pancan TCGA data
rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")

library(nilsTCGA)


df2Grange <- function(dfName){
  suppressPackageStartupMessages(library("GenomicRanges"))
  
  GRange_object <- GRanges(seqnames =  dfName$Chromosome,
                         ranges=IRanges(dfName$Start,dfName$End),
                         strand='*',
                         
                         Num_Probes=dfName$Num_Probes,
                         Segment_Mean=dfName$Segment_Mean,
                         Sample=dfName$Sample)
  
  return (GRange_object)
  
  
}


# try to compare the different data tables with my original
#CNV_table[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-01$", "", CNV_table[,1])
#B<-CNV_table[CNV_table$Sample == "TCGA-OR-A5KQ-01A",]

#readCNVs('PANCAN_SNP.seg')
CNV_table <- read.table('PANCAN_CNV.seg',header = TRUE)

# Extract what sample types are in the SNP_table
types <- gsub("TCGA-[A-Z0-9]*-[A-Z0-9]*-", "", (gsub("[A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CNV_table[,1])))
unique(types)



# What types match to
#normal <- c("10A", "10B", "11A", "11B")
#tumor <- c("01A", "01B", "06A")
tumor <- c("01","06","02","05")
normal <- c("10","11","12","14")

CN_tumor <- CNV_table[which(types %in% tumor),];
CN_normal <- CNV_table[which(types %in% normal),];


# Change ID to match with other data-types
CN_tumor[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN_tumor[,1])
CN_normal[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN_normal[,1])


# Convert CN_tumor into GRange object

CN_tumor_Grange<-df2Grange(CN_tumor)

# Find all segments where the Segment_Mean < -2
potential_loss <- CN_tumor_Grange[CN_tumor_Grange@elementMetadata$Segment_Mean < -2,]


#Find all genes that are expressed in these regions (Script from Malin)

GeneRegions <- function (genelist) {
  

  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(biomaRt))
  suppressPackageStartupMessages(library(BSgenome))
  suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE))
  
  #From Malin, get all genes and the Granges positions stored in 'genelistgranges'
  mart <- "ensembl"
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") # OBS !!! What version was used for the annotation? Default GRCh38  
  ensemblgenes <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name','start_position','end_position'), mart = ensembl)
  
  xidx <- which(ensemblgenes[,3]=="X")
  yidx <- which(ensemblgenes[,3]=="Y")
  ensemblgenes[xidx, 2] <- 23
  ensemblgenes[yidx, 2] <- 24
  
  genelistgranges <- GRanges(seqnames=ensemblgenes$chromosome_name, ranges=IRanges(ensemblgenes$start_position, ensemblgenes$end_position), hgnc_symbol=ensemblgenes$hgnc_symbol
                             ,ensemble_gene_id = ensemblgenes$ensembl_gene_id)
  
  

  keeplist <- c(1:24)
  tokeep<- keeplist[which(keeplist %in% levels(factor(seqnames(genelistgranges))))]
  genelistgranges<- keepSeqlevels(genelistgranges,tokeep, pruning.mode="coarse")
  
}
 



GenelistGrange <- GeneRegions()
  
# function that finds all proteins in the potential_loss regions for each indivdual
#####  Subject = Potenital_loss
#####  Query   = GeneList
#hits returned indices over which samples have deleted genes.
hits <- findOverlaps(GenelistGrange,potential_loss, type="within")  
 
 
potential_loss <- potential_loss[subjectHits(hits)]
potential_loss@elementMetadata$DeletedGene <- GenelistGrange[queryHits(hits)]@elementMetadata$ensemble_gene_id
#GenelistGrange[queryHits(hits)]$ensemble_gene_id  Gets all the gene id:s

#in total 12692 unique genes over all samples


#---------------------------------m-RNA section ------------------------------------------
 
#Way to large file for my computer to handle
#tab5rows <- read.table("mRNA.geneEXP.tsv", header = TRUE, nrows = 2)
#classes <- sapply(tab5rows, class)
#tabAll <- read.table("mRNA.geneEXP.tsv", header = TRUE, colClasses = classes)





#
#mRNA_table <- read.table("mRNA.geneEXP.tsv", header=TRUE,sep = "\t", nrows=2)
#library(readr)
#tab100 <- read_tsv('mRNA.geneEXP.tsv',n_max=100)



