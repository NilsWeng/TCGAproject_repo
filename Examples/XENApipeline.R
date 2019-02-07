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

#test

CNV_table <- read.table('Xena-GDC-PANCAN.masked_cnv.tsv',header = TRUE)
#xidx <- which(CNV_table$Chromosome=="X")
#CNV_table[xidx , 2] <- 23


# Convert CN_tumor into GRange object

CNV_Grange<-df2Grange(CNV_table)

# Find all segments where the Segment_Mean < -2
potential_loss <- CNV_Grange[CNV_Grange@elementMetadata$Segment_Mean < -2,]


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


plot_CN(CNV_table)



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
tab100 <- read_tsv('mRNA.geneEXP.tsv',n_max=4)
looking_for_gene <- c("?|100134869","?|10357") # test-vector for script
looking_for_patient <- c("TCGA-A13J-11A","TCGA-A5JA-01A")
M <- matrix(c(looking_for_gene,looking_for_patient),ncol=2)
colnames(M) <- (c("DeletedGene","Sample"))
tab100_hits <- tab100[tab100$gene_id %in% M[,1],]
name_vector <- names(tab100_hits)
#grep("TCGA-[0-9A-Z]*-[0-9A-Z]*",name_vector)
name_vector[2:length(name_vector)] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*","",name_vector[2:length(name_vector)])

#rowMeans(tab100[, -(1)],na.rm=TRUE)
# Problem sample ID is expressed in long format in the mRNA file <- Fix this !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


row_mean <- rowMeans(tab100_hits[, -(1)],na.rm=TRUE)
#sample_expression <- 
intervall <- 1:length(tab100_hits$gene_id)
for (i in intervall) {
  
  A <- M[i,'Sample']
  expression <- tab100_hits[i,A]
  p_val <- expression / row_mean[i]
  
}

# Want to find the expression level given a sample and ENSGnumber



