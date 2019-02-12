rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")

#library(BiocParallel)
library(BSgenome)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(MutationalPatterns)

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


# Read mc3-file

MC3_DF <- read.table('mc3.v0.2.8.PUBLIC.maf.gz',header=TRUE,sep='\t',nrow=200)