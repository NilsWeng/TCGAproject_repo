rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")

TSS2Study_DF <- read.table("TSS2Studyabb.txt",header=TRUE)




#library(BiocParallel)
library(BSgenome)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(MutationalPatterns)
library(plyr)

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


# Read mc3-file (From PanCan)
MC3_DF <- read.table('mc3.v0.2.8.PUBLIC.maf.gz',header=TRUE,sep='\t',nrow=1000);
MC3_DF$TSS.Code <-gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",gsub("TCGA-","",MC3_DF$Tumor_Sample_Barcode))
MC3_DF <- join(MC3_DF,TSS2Study_DF, by="TSS.Code")

#MC3_DF <- read.table('mc3.v0.2.8.PUBLIC.code.maf',header=TRUE,sep='\t',nrow=200);


#Save mc3-file in separate folders as VCF-files

MC3_DF_CancerType <- split (MC3_DF, f=MC3_DF$Study.Abbreviation,drop=TRUE) # mayby i need drop = TRUE ?




vcffile=vector()
subtype=vector()
sample=vector()



#Start of loops --------------------------------------------

for (cancerType in MC3_DF_CancerType){ # Loop over all cancer types
  
  #Add a catch that prevents closing loop
  if (dim(cancerType)[1]  == 0) {
    
    next()
  
  }
    

  
  
  directory_name <- paste("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files","/",unique(cancerType$Study.Abbreviation),sep="") 
  
  #crate folder if it doesnt exist
  ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
  #setwd(directory_name) 
  print (directory_name)
  setwd(directory_name)
  
  
  #Split cancer type into individual samples
  MC3_DF_Sample <- split (cancerType, f= cancerType$Tumor_Sample_Barcode,drop = TRUE) # not sure if this works
  
  
  
  
  

    for (Sample in MC3_DF_Sample){ #Loop over each sample in respective cancer type
      
      
      if (dim(Sample)[1]  == 0) { ########### Why does CancerType get split into empty samples??????????!!!!!!!!
                                  #Seem to split MC3_DF_Sample but still save old levels  # try drop = TRUE in split command
        
        next()
        
      } 
      
      
    
    sample_id <- unique(Sample$Tumor_Sample_Barcode);
    vcfdata <- matrix(".", nrow=nrow(Sample), ncol = 10);
    columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",sample_id )
    colnames(vcfdata) <- columns
    
    vcfdata[,1]=as.character(Sample$Chromosome)
    vcfdata[,2]=as.character(Sample$Start_Position)
    vcfdata[,4]=as.character(Sample$Reference_Allele)
    vcfdata[,5]=as.character(Sample$Tumor_Seq_Allele2)
    vcfdata[,9]="GT"
    vcfdata[,10]="1/0"
    vcfdata[,3]=as.character(sample_id)
    outfile=paste(sample_id,".vcf", sep="")
    vcffile <- c(vcffile, outfile)
    #outfile = paste(directory_name,outfile, sep="/")
    write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
    
    
    
    
  }
  
  
}


#Start working on subtype and sample!

#END OF LOOPS----------------------------------------------







#Split table at sample id
MC3_separated <- split( MC3_DF , f = MC3_DF$Tumor_Sample_Barcode )



vcffile=vector()
subtype=vector()
sample=vector()





































#Approach 1 ----------------------------------------------------

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3output")

for (sample_maf in MC3_separated) {
  
  
  
  sample_id <- unique(sample_maf$Tumor_Sample_Barcode);
  vcfdata <- matrix(".", nrow=nrow(sample_maf), ncol = 10);
  columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",sample_id )
  colnames(vcfdata) <- columns
  
  vcfdata[,1]=as.character(sample_maf$Chromosome)
  vcfdata[,2]=as.character(sample_maf$Start_Position)
  vcfdata[,4]=as.character(sample_maf$Reference_Allele)
  vcfdata[,5]=as.character(sample_maf$Tumor_Seq_Allele2)
  vcfdata[,9]="GT"
  vcfdata[,10]="1/0"
  vcfdata[,3]=as.character(sample_id)
  outfile=paste(sample_id,".vcf", sep="")
  write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
  vcffile <- c(vcffile, outfile)
  
  
  
}


auto = extractSeqlevelsByGroup(species="Homo_sapiens", 
                               style="UCSC",
                               group="auto")


read_vcfs_as_granges(vcffile,names,genome)




#Approach 2 - acess MC3 over and over through loop

#number_of_samples <- levels(factor(MC3_DF$Tumor_Sample_Barcode))

#for (sample_id in number_of_samples) { # Loop over number of samples
  
  
  #tumormaf=MC3_DF[which(MC3_DF$Tumor_Sample_Barcode %in% sample_id),]
  
  
  
  
  
  
  
#}



#Extract all unqiue samples and store as vcf-file





