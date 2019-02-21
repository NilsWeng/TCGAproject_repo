#Double substitution mutations
#Script that reads VCF file and extract all double base substitution (DBS) mutations
rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(data.table)
library(GenomicRanges)



#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


vcf_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern="cohort.vcf")

vcf_list <- "CHOL/CHOL.cohort.vcf" # Reduce sample size for testing

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")


#Read VCF file
VCF <- fread(vcf_list)


#Find all neighbouring mutations

#Sort VCF file
VCF <- VCF[with(VCF, order(ID,as.numeric(`#CHROM`),POS)), ]


diff_vector <- as.numeric(diff(VCF$POS))
DBS_pos <- which(diff_vector == 1)
DBS_pos <- c(DBS_pos,DBS_pos+1)
DBS_pos <- sort.default(DBS_pos)
#Find all 
DBS_mut <- VCF[DBS_pos ,]



#Just to get it to G-range object......
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/DBS")
write.table(DBS_mut, file="test.vcf", row.names=FALSE, sep="\t", quote=FALSE)
DBS_mut_g <- read_vcfs_as_granges("test.vcf","test",ref_genome)


#find 3' and 5' surrounding

ranges = resize(DBS_mut_g, 3, fix = "center")





#mut_surrounding <- function()


mut_surrounding <- function(vcf,ref_genome){
  
  #Only correct for first base
  vcf <- vcf[[1]]
  
  ranges = resize(vcf, 3, fix = "center")
  vcf_context = as.character(getSeq(get(ref_genome),
                                    seqnames(vcf),
                                    start(vcf) - 1,
                                    end(vcf) + 2))
  
  return(vcf_context)
  
}




mut_surround <- mut_surrounding(DBS_mut_g,ref_genome)
mut_surround <- mut_surround[seq(1,length(mut_surround),2)]
mut_surround <- rep(as.vector(mut_surround),each=2)


DBS_mut$DBS_type <- mut_surround

#Test <- mutations_from_vcf(VCF)

#GRANGE <- read_vcfs_as_granges(vcf_list,"test",ref_genome)



