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






# From GRange object ------------------------------------------------------
rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(data.table)
library(GenomicRanges)
library(plyr)


#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


vcf_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern="cohort.vcf")

vcf_list <- tail(vcf_list,2)
DBS_freq_table <- data.frame()



for (cancer in vcf_list){
  
  print(paste("Starting on ",cancer,sep=""))
 
  cancer_name <- gsub("/..*","",cancer)
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
  Grange_vcf <- read_vcfs_as_granges(cancer,cancer_name,ref_genome)
  Grange_vcf <- Grange_vcf[[1]]
  
  
  
  
  
  
  
  DBS_mutation <- function(Grange_vcf,ref_genome){
    
    Grange_vcf <- sortSeqlevels(Grange_vcf)
    
    
    diff_vector <- as.numeric(diff(Grange_vcf@ranges@start))
    DBS_pos <- which(diff_vector == 1)
    DBS_pos <- c(DBS_pos,DBS_pos+1)
    DBS_pos <- sort.default(DBS_pos)
    #Find all 
    DBS_mut_g <- Grange_vcf[DBS_pos ,]
    
    #find 3' and 5' surrounding
    
    ranges = resize(DBS_mut_g, 3, fix = "center")
    
    
    
    
    
    #mut_surrounding <- function()
    
    
    mut_surrounding <- function(vcf_g,ref_genome){
      
      #Only correct for first base
      
      
      ranges = resize(vcf_g, 3, fix = "center")
      vcf_context = as.character(getSeq(get(ref_genome),
                                        seqnames(vcf_g),
                                        start(vcf_g) - 1,
                                        end(vcf_g) + 2))
      
      return(vcf_context)
      
    }
    
    
    
    
    mut_surround <- mut_surrounding(DBS_mut_g,ref_genome)
    mut_surround <- mut_surround[seq(1,length(mut_surround),2)]
    mut_surround <- rep(as.vector(mut_surround),each=2)
    
    from <- mut_surround
    to   <- data.frame(DBS_mut_g$ALT)
    to <- as.vector(to$value)
    to <- paste(to[seq(1,length(to),2)],to[seq(2,length(to),2)],sep="")
    to <- rep(as.vector(to),each=2)
    to <- paste(gsub("[A-Z][A-Z][A-Z]$","",from),to,sep="")
    to <- paste(to,gsub("^[A-Z][A-Z][A-Z]","",from),sep="")
    return_matrix <- data.frame(from,to,DBS_mut_g@ranges@start,DBS_mut_g@seqnames,DBS_mut_g@ranges@NAMES)
    return(unique(return_matrix))
    
    
    
    
  }
  
  
  DBS_DF <- DBS_mutation(Grange_vcf,ref_genome)
  
  
  #Make a list of what types of DBS mutations are most common
  
  from_mut <- substr(as.character(DBS_DF$from),2,3)
  to_mut <-  substr(as.character(DBS_DF$to),2,3)
  mut_type <- paste(from_mut,to_mut,sep=">")
  mut_type <- mut_type[seq(1,length(mut_type),2)]
  context <- paste(substr(as.character(DBS_DF$from),1,1),substr(as.character(DBS_DF$from),4,4),sep="..")
  context <- context[seq(1,length(context),2)]
  
  DF <- data.frame(mut_type,context)
  
  DF <- count(DF)
  DF$Cancer <- cancer_name
  DBS_freq_table <- rbind(DBS_freq_table,DF)
  #break()
  
  
  
  #Save data
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/DBS/DBS_data")
  #Save Grange-vcf file
  save(Grange_vcf,file=paste(cancer_name,"_grange_vcf.rda",sep=""))
  
  write.table(DF,file=paste(cancer_name,"_DBS_DF.txt",sep=""))
  
  
  
  
  
  
}



library(dplyr)
tot_DBS <- ALL_DBS %>% group_by(mut_type) %>% summarise(freq = sum(freq))
tot_DBS %>% arrange(desc(freq))




































Grange_vcf <- read_vcfs_as_granges(vcf_list,"CHOL",ref_genome)
Grange_vcf <- Grange_vcf[[1]]







DBS_mutation <- function(Grange_vcf,ref_genome){
  
  Grange_vcf <- sortSeqlevels(Grange_vcf)
  
  
  diff_vector <- as.numeric(diff(Grange_vcf@ranges@start))
  DBS_pos <- which(diff_vector == 1)
  DBS_pos <- c(DBS_pos,DBS_pos+1)
  DBS_pos <- sort.default(DBS_pos)
  #Find all 
  DBS_mut_g <- Grange_vcf[DBS_pos ,]
  
  #find 3' and 5' surrounding
  
  ranges = resize(DBS_mut_g, 3, fix = "center")
  
  
  
  
  
  #mut_surrounding <- function()
  
  
  mut_surrounding <- function(vcf_g,ref_genome){
    
    #Only correct for first base
    
    
    ranges = resize(vcf_g, 3, fix = "center")
    vcf_context = as.character(getSeq(get(ref_genome),
                                      seqnames(vcf_g),
                                      start(vcf_g) - 1,
                                      end(vcf_g) + 2))
    
    return(vcf_context)
    
  }
  
  
  
  
  mut_surround <- mut_surrounding(DBS_mut_g,ref_genome)
  mut_surround <- mut_surround[seq(1,length(mut_surround),2)]
  mut_surround <- rep(as.vector(mut_surround),each=2)
  
  from <- mut_surround
  to   <- data.frame(DBS_mut_g$ALT)
  to <- as.vector(to$value)
  to <- paste(to[seq(1,length(to),2)],to[seq(2,length(to),2)],sep="")
  to <- rep(as.vector(to),each=2)
  to <- paste(gsub("[A-Z][A-Z][A-Z]$","",from),to,sep="")
  to <- paste(to,gsub("^[A-Z][A-Z][A-Z]","",from),sep="")
  return_matrix <- data.frame(from,to,DBS_mut_g@ranges@start,DBS_mut_g@seqnames,DBS_mut_g@ranges@NAMES)
  return(unique(return_matrix))
  
  
  
  
}


DBS_DF <- DBS_mutation(Grange_vcf,ref_genome)


#Make a list of what types of DBS mutations are most common

from_mut <- substr(as.character(DBS_DF$from),2,3)
to_mut <-  substr(as.character(DBS_DF$to),2,3)
mut_type <- paste(from_mut,to_mut,sep=">")
mut_type <- mut_type[seq(1,length(mut_type),2)]





DBS_freq_table <- rbind(DBS_freq_table,count(mut_type))








#Create a mutational matrix for DBS
Context<-c("A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T")

Type <- c("AA","AC","AG","AT","CA","AC","AG","AT","AA","AC","AG","AT")
















