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
#CNV_table <- read.table('Xena-GDC-PANCAN.masked_cnv.tsv',header = TRUE)

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


# Sort the potential_loss object

potential_loss <- sortSeqlevels(potential_loss)
potential_loss <- sort(potential_loss, by = ~  Sample + seqnames)

#in total 12692 unique genes over all samples

#-----------------Search if any of the potential_loss candidate genes are found in list of known mutational signatures genes------------


#geneINmutlist takes a GRange obejct with genes and match it agains list of potential mutsign genes and returns matching genes as GRange object
geneINmutlist <- function(GrangeObject){
  
  
  
  mutsign_genes <- read.table("mutsign_genes.csv",header=TRUE,sep=";")

  indices <- potential_loss@elementMetadata$DeletedGene %in% mutsign_genes$Gene.ID
  found_genes <- potential_loss[indices]
  
  return (found_genes)
  
}

found_genes <- geneINmutlist(potential_loss)


#---------------------------------------MC3-MAF sectio -----------------------------------

#setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")
#MC3_DF <- read.table('mc3.v0.2.8.PUBLIC.maf.gz',header=TRUE,sep='\t',nrow=200)












#---------------------------------m-RNA section ------------------------------------------
 
#Way to large file for my computer to handle
#tab5rows <- read.table("mRNA.geneEXP.tsv", header = TRUE, nrows = 2)
#classes <- sapply(tab5rows, class)
#tabAll <- read.table("mRNA.geneEXP.tsv", header = TRUE, colClasses = classes)



# Contains 11070 columns i.e one for each patient 
#format TCGA-BL-A13J-11A-13R-A10U-07

#
#mRNA_table <- read.table("mRNA.geneEXP.tsv", header=TRUE,sep = "\t", nrows=2)
library(readr)
tab100_hits <- read_tsv('mRNA.geneEXP.tsv')
name_vector <- names(tab100_hits)
name_vector[2:length(name_vector)] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*$","",name_vector[2:length(name_vector)])
names(tab100_hits) <- name_vector



# Subset where only genes found in found_genes are selected
#tab100_hits <- tab100[tab100$gene_id %in% M[,1],]
#name_vector <- names(tab100_hits)
#name_vector[2:length(name_vector)] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*$","",name_vector[2:length(name_vector)])
#names(tab100_hits) <- name_vector

# Get sample expression and avarage expression for candidates aswell as expression/avarage as P: Returns(exp,avg_exp,P)

get_gene_expression <- function(Gene_id,Sample_id){
  
  col_name = Sample_id
  row_name = Gene_id
  T <- tab100_hits$gene_id == row_name
  A <- (tab100_hits[T ,])
  
  
  expression <- as.numeric(A[, col_name])
  

  
  avg_expression <-  as.numeric(rowMeans(A[, -(1)],na.rm=TRUE))
  
  P <- expression/avg_expression
  return_vector <- c(expression,avg_expression,P)
  #colnames(return_vector) <- c('exp','avg_exp')
  return(return_vector)
  
}

###TEST SECTION------------------------------------------------------
#Should give 66.9 and 2,287


X <- c("TCGA-OR-A5J6","TCGA-OR-A5J8")
Y <- c("?|10357","?|100134869")
M <- cbind(X,Y)

get_gene_expression("BRCA2|675","TCGA-KD-A5QS")
get_gene_expression("PMS2|5395","TCGA-13-15OO")
# Do get_gene_expression to each 
mapply(get_gene_expression,Y,X)
###--------------------------------------------------------------------------
##### Problem -- gene_id expressed as ?|1001~


Sample_id <- found_genes$Sample
gene_id <- found_genes$DeletedGene

mapply(get_gene_expression,gene_id,Sample_id)




# Want to find the expression level given a sample and ENSGnumber



