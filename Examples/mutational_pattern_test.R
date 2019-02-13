rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")

TSS2Study_DF <- read.table("TSS2Studyabb.txt",header=TRUE)




#library(BiocParallel)
library(BSgenome)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(MutationalPatterns)
library(plyr)
library(data.table)
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


# Read mc3-file (From PanCan)

#MC3_DF <- read.table('mc3.v0.2.8.PUBLIC.maf.gz',header=TRUE,sep='\t',nrow=4000); #Has a total of 3600964 rows Around 24 min to read whole file





start_time <- Sys.time()


MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz',nrow=100000)

end_time <- Sys.time()
print("Read time")
end_time - start_time

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
      outfile <- paste(sample_id,".vcf", sep="")
      #outfile = paste(directory_name,outfile, sep="/")
      write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
      
      subtype <- c(subtype,as.character(unique(cancerType$Study.Abbreviation)))
      sample <- c(sample,as.character(sample_id))
      vcffile <- c(vcffile, paste (unique(cancerType$Study.Abbreviation),outfile,sep = "/") )
    
  }
  
  
}




#END OF LOOPS----------------------------------------------






#Save cohort_level vcf files.
for (cancerType in MC3_DF_CancerType) {
  
  #Add a catch that prevents closing loop
  if (dim(cancerType)[1]  == 0) {
    
    next()
    
  }
  
  directory_name <- paste("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files","/",unique(cancerType$Study.Abbreviation),sep="") 
  
  
  #crate folder if it doesnt exist
  ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
  print (directory_name)
  setwd(directory_name)
  
  #Save data to vcf file
  
  
  cancer_abb <- unique(cancerType$Study.Abbreviation);
  vcfdata <- matrix(".", nrow=nrow(cancerType), ncol = 10);
  columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "WHAT??" )
  colnames(vcfdata) <- columns
  
  vcfdata[,1]=as.character(cancerType$Chromosome)
  vcfdata[,2]=as.character(cancerType$Start_Position)
  vcfdata[,4]=as.character(cancerType$Reference_Allele)
  vcfdata[,5]=as.character(cancerType$Tumor_Seq_Allele2)
  vcfdata[,9]="GT"
  vcfdata[,10]="1/0"
  vcfdata[,3]=as.character(cancerType$Tumor_Sample_Barcode)
  
  outfile <- paste(cancer_abb,".cohort.vcf", sep="")
  #outfile = paste(directory_name,outfile, sep="/")
  write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
  

  
}



cohort <- data.frame("vcf"=vcffile, "subtype"=subtype, "sample"=sample)

#Count number of mutations and add to cohort table
library(R.utils)
nmut_total=sapply(cohort$vcf,countLines)-1
cohort$nmut <- nmut_total

#sum(cohort[cohort$subtype=="OV",]$nmut)




##############   Start  Cohort section-------------------------------------------------------------------------------
##############-------------------------------------------------------------------------------------------------------

#Read all cohort-level files and transform to G-range objects

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")

vcf_cohort_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern="cohort.vcf")
name_cohort <- gsub("/.*","",vcf_cohort_list)

#Read each cohort.vcf file as Grange oject
cohort_vcfs <- read_vcfs_as_granges(vcf_cohort_list,name_cohort,ref_genome)

# Create mutational matrix for each cancer type
cohort_mut_matrix = mut_matrix(cohort_vcfs,ref_genome)

#Extract mutational signatures
extract_3 = extract_signatures(cohort_mut_matrix, rank = 4) # Cant extract more signatures than number of cancer types



# Plot extracted signatures and their contribution

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Pictures/Cohort")

pdf("extracted_3.pdf")
plot_96_profile(extract_3$signatures)
dev.off()
pdf("extracted_contribution.pdf")
plot_contribution(extract_3$contribution, extract_3$signature, mode = "absolute", coord_flip = T)
dev.off()

####### End cohort section-------------------------------------------------------------
#######--------------------------------------------------------------------------------




#Open VCF-files and transform into G-ranges object-----------------------------------------------

library(MutationalPatterns)

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")


auto = extractSeqlevelsByGroup(species="Homo_sapiens", 
                               style="UCSC",
                               group="auto")



vcf_file_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern=".vcf")
read_vcfs_as_granges(cohort_table$vcf,cohort_table$sample,ref_genome)




########## ---- Mutational signatures for each sample for one cancer type  #######
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")



file_list <- as.character(cohort$vcf[which(cohort$subtype %in% "OV")])
name_list <- as.character(cohort$sample[which(cohort$subtype %in% "OV")])

GRange_VCF <- read_vcfs_as_granges(file_list,name_list,ref_genome)
GRange_VCF_auto = lapply(GRange_VCF, function(x) keepSeqlevels(x, auto, pruning.mode='coarse'))


mut_matrix_OV <- mut_matrix(GRange_VCF_auto,ref_genome)
mut_signatures_OV <- extract_signatures(mut_matrix_OV,rank=20)


#Rank optimisation
if (FALSE){
  
  library(NMF)
  #Add pseudocount to mut matrix ?
  #mut_mat <- mut_mat + 0.0001
  estim.r <- nmf(mut_matrix_OV, rank=4:20,method="brunet", nrun=10, seed=123456)
  plot(estim.r)

}



#Plot results of signature extraction

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Pictures/OV")

pdf("mut_signatures_OV1.pdf")
plot_96_profile(mut_signatures_OV$signatures,condensed = TRUE)
dev.off()
pdf("extracted_contribution.pdf") # Mayby plot in chunks of 100? Hard to see anything
plot_contribution(mut_signatures_OV$contribution, mut_signatures_OV$signature, mode = "absolute", coord_flip = T)
dev.off()


#Compare similarity of extracted results to Cosmic signatures


#Download Cosmic-signatures
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
                 "signatures_probabilities.txt", sep = "")
cosmic_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_matrix_OV), cosmic_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cosmic_signatures = cosmic_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cosmic_signatures) = cosmic_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cosmic_signatures = as.matrix(cosmic_signatures[,4:33])



#Similarity between cosmic_signatures and mut_matrix_OV
cos_sim_samples_signatures <- cos_sim_matrix(mut_matrix_OV, cosmic_signatures)


hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]

plot_cosine_heatmap(cos_sim_samples_signatures,
                     col_order = cosmic_order,
                     cluster_rows = TRUE)



#Reconstruct mut_matrix_OV by linear combination of cosmic_signatures

fit_res <- fit_to_signatures(mut_matrix_OV, cosmic_signatures)
select <- which(rowSums(fit_res$contribution) > 10)
#Plot contribution barplot
plot_contribution(fit_res$contribution[select,],
                       cosmic_signatures[,select],
                       coord_flip = FALSE,
                       mode = "absolute")



##########---------------------------------------------------------------------


  









  

