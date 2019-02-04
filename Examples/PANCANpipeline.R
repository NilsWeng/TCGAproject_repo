#Test of handling pancan TCGA data
rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")

library(nilsTCGA)


# try to compare the different data tables with my original
#CNV_table[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-01$", "", CNV_table[,1])
#B<-CNV_table[CNV_table$Sample == "TCGA-OR-A5KQ-01A",]

#readCNVs('PANCAN_SNP.seg')
CNV_table <- read.table('PANCAN_CNV.seg',header = TRUE)

# Extract what sample types are in the SNP_table
types <- gsub("TCGA-[A-Z0-9]*-[A-Z0-9]*-", "", (gsub("-[0-9A-Z]*-[0-9A-Z]*-01$", "", CNV_table[,1])))
unique(types)


# What types match to
normal <- c("10A", "10B", "11A", "11B")
tumor <- c("01A", "01B", "06A")

CN_tumor <- CNV_table[which(types %in% tumor),]
CN_normal <- CNV_table[which(types %in% normal),]


# Change ID to match with other data-types
CN_tumor[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN_tumor[,1])
CN_normal[,1] <- gsub("-[0-9A-Z]*-[0-9A-Z]*-[0-9A-Z]*-01$", "", CN_normal[,1])


#
mRNA_table <- read.table("mRNA.geneEXP.tsv",header=TRUE,sep = "\t", nrows=2)

seg <- read.table("TCGA_mastercalls.abs_segtabs.fixed.txt",header=TRUE,sep = "\t")

seg2 <-  read.table("TCGA_mastercalls.abs_tables_JSedit.fixed.txt",header=TRUE,sep = "\t")
