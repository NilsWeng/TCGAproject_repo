#Test of handling pancan TCGA data
rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")

library(nilsTCGA)


df2Grange <- function(dfName){
  suppressPackageStartupMessages(library("GenomicRanges"))
  
  GRange_object <- GRanges(seqnames =  dfName$Chromosome,
                           ranges=IRanges(dfName$Start,dfName$End),
                           strand='*',
                           Segment_Mean=dfName$Segment_Mean,
                           Sample=dfName$Sample)
  
  return (GRange_object)
  
  
}



CNV_table <- read.table('Xena-GDC-PANCAN.masked_cnv.tsv',header = TRUE)
#xidx <- which(CNV_table$Chromosome=="X")
#CNV_table[xidx , 2] <- 23




# Convert CN_tumor into GRange object

CNV_Grange<-df2Grange(CNV_table)

# Find all segments where the Segment_Mean < -2
potential_loss <- CNV_Grange[CNV_Grange@elementMetadata$Segment_Mean < -1,]


###--------------------------------Plotting CNVs section ---------------------------------------

library(gaia)
library(TCGAbiolinks)

#gaiaCNVplot(CNV_table[, 2:5],threshold = 0.5)



CNV_table <- CNV_table[CNV_table$Sample == 'TCGA-23-1023-01R' ,]





plot_CN <- function(CN_table){
  
  Chromo <- CNV_table$Chromosome
  CN <- CNV_table$Segment_Mean
  start <- CNV_table$Start
  end <- CNV_table$End
  
  
  
  plot (CN,
        ylim = c(-3, max(abs(CNV_table$Segment_Mean)+2)),
        type = "h",
        col = "red",
        xlab = "Chromosome",
        ylab = "Segment_mean",
        xaxt = "n")
  
  
  
  uni.chr <- unique(Chromo)
  temp <- rep(0, length(uni.chr))
  for (i in 1:length(uni.chr)) {
    temp[i] <- max(which(uni.chr[i] == Chromo))
  }
  for (i in 1:length(temp)) {
    abline(v = temp[i], col = "black", lty = "dashed")
  }
  nChroms <- length(uni.chr)
  begin <- c()
  for (d in 1:nChroms) {
    chrom <- sum(Chromo == uni.chr[d])
    begin <- append(begin, chrom)
  }
  temp2 <- rep(0, nChroms)
  for (i in 1:nChroms) {
    if (i == 1) {
      temp2[1] <- (begin[1] * 0.5)
    }
    else if (i > 1) {
      temp2[i] <- temp[i - 1] + (begin[i] * 0.5)
    }
  }
  
  uni.chr[uni.chr==23] <- "X"
  uni.chr[uni.chr==24] <- "Y"
  for (i in 1:length(temp)) {
    axis(1, at = temp2[i], labels = uni.chr[i], cex.axis = 1)
  }
  legend(x=1,y=max(Calls[,"score"]+2), y.intersp=0.8, c("Amp"), pch=15, col=c("red"), text.font=3)
  legend(x=1,y=-max(Calls[,"score"]+0.5), y.intersp=0.8, c("Del"), pch=15, col=c("blue"), text.font=3)
  
}


#plot_CN(CNV_table)



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
potential_loss@elementMetadata$DeletedGene <- GenelistGrange[queryHits(hits)]@elementMetadata$ensemble_gene_id  # Change between HGNC or ensemble
potential_loss@elementMetadata$DeletedGeneHGNC <- GenelistGrange[queryHits(hits)]@elementMetadata$hgnc_symbol

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

# Contains 11070 columns i.e one for each patient 
#format TCGA-BL-A13J-11A-13R-A10U-07


library(readr)
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")

EXP_table <- read_tsv('mRNA.geneEXP.tsv')




get_hgcn <- function(name){
  
  name <-  strsplit(name, "\\|")[[1]][1]
  return (name)
}

hgcn_Name <- sapply(EXP_table$gene_id,get_hgcn)
unname(hgcn_Name)
EXP_table$gene_id <- hgcn_Name


#Filter for gene-expression that was found in found_genes
EXP_hits <- EXP_table[EXP_table$gene_id %in% found_genes$DeletedGeneHGNC ,]




name_vector <- names(EXP_hits)
name_vector[2:length(name_vector)] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*$","",name_vector[2:length(name_vector)])
names(EXP_hits) <- name_vector




# Want to find the expression level given a sample and ENSGnumber

get_gene_expression <- function(Gene_id,Sample_id){
  
  
  #col_name ="TCGA-OR-A5J3-01A"
  #row_name = "PMS2"
  
  
  col_name = Sample_id
  row_name = Gene_id
  T <- EXP_hits$gene_id == row_name
  A <- (EXP_hits[T ,])
  
  
  expression <- as.numeric(A[, col_name])
  
  
  
  avg_expression <-  as.numeric(rowMeans(A[, -(1)],na.rm=TRUE))
  
  P <- expression/avg_expression
  return_vector <- c(expression,avg_expression,P)
  names(return_vector) <- c('exp','avg_exp',"P")
  return(return_vector)
  
}


sample_id <- as.vector(found_genes$Sample)

gene_id <- found_genes$DeletedGeneHGNC

DF <- data.frame(sample_id,gene_id)

#Crashes because not all samples are present in both.
DF <- DF[DF$sample_id %in% names(EXP_hits),]

##################### Matrix is not correct, some samples are not present in names(EXP_hits) that are found in found_genes$sample
######checked for  %in% and selected those , however i think the got shifted and its all wrong



Matrix <- as.data.frame(mapply(get_gene_expression,DF$gene_id,as.vector(DF$sample_id)))
colnames(Matrix) <- paste(as.vector(DF$sample_id),DF$gene_id, sep= "<-->")
Matrix

