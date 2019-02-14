#Clear workspace
rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(plyr)
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

#Load cosmic_signatures from file in folder
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cosmic_signatures <- as.matrix(read.table("cosmic_signatures.txt",header=TRUE))


#Set wd
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")









# Cohort section ----------------------------------------------------------

#Read all vcf_cohort files
vcf_cohort_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern="cohort.vcf")

#Optional: Select only certain cancertypes
vcf_cohort_list <- vcf_cohort_list[gsub("/.*","",vcf_cohort_list) %in% c('CHOL','ACC','USC','UVM','DLBC','MESO')] 


name_cohort <- gsub("/.*","",vcf_cohort_list)


#Read each cohort.vcf file as Grange oject
cohort_vcfs <- read_vcfs_as_granges(vcf_cohort_list,name_cohort,ref_genome)

#Create mutational matrix for each cancer type
cohort_mut_matrix = mut_matrix(cohort_vcfs,ref_genome)


#optional, extract optimal rank for NMF
if (FALSE){
  
  library(NMF)
  #Add pseudocount to mut matrix ?
  #mut_mat <- mut_mat + 0.0001
  estim.r <- nmf(cohort_mut_matrix, rank=2:8,method="brunet", nrun=10, seed=123456)
  plot(estim.r)
  
}


#Extract mutational signatures (Cant extract more signatures than given samples)
number_of_signatures <- 5

cohort_extracted_signatures = extract_signatures(cohort_mut_matrix, number_of_signatures) 
colnames(cohort_extracted_signatures$signatures) <- c(1:number_of_signatures)

#Plot extracted signatures
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Pictures/Cohort")

pdf("Extracted_signatures.pdf")
plot_96_profile(cohort_extracted_signatures$signatures, condensed = TRUE)
dev.off()




#Compare extracted signatures with cosmic signatures
cosine_similarity = cos_sim_matrix(cohort_extracted_signatures$signatures, cosmic_signatures)
# Plot heatmap with specified signature order
hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]



pdf("Extracted_vs_Cosmic_cosine.pdf")
plot_cosine_heatmap(cosine_similarity, col_order = cosmic_order, cluster_rows = TRUE)
dev.off()


#Plot contribution of extracted signature vs cancer
pdf("Cancer_vs_Contribution_extracted_signatures.pdf")
plot_contribution(cohort_extracted_signatures$contribution, cohort_extracted_signatures$signature, mode = "absolute", coord_flip = TRUE)
dev.off()

#Plot contribution of cosmic vs cancer

fit_res <- fit_to_signatures(cohort_mut_matrix, cosmic_signatures)
# Select signatures with some contribution
select <- which(rowSums(fit_res$contribution) > 10)
# Plot contribution barplot

pdf("Cancer_vs_contribution_Cosmic_signatures.pdf")
plot_contribution(fit_res$contribution[select,],
                       cosmic_signatures[,select],
                       coord_flip = TRUE,
                       mode = "absolute")

dev.off()




# Section One cancer type at a time ---------------------------------------

#load cohort_table.txt
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
Cohort <- read.table("cohort.txt",header = TRUE)#maybe unness could just aswell list all files in folder (except cohort.vcf files)
#list.dirs(path = "C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files", full.names = TRUE, recursive = TRUE)



#loop over each cancer type


cancer_list <- as.vector(unique(Cohort$subtype))


for (cancer_type in cancer_list){
  
  
   
  print(paste("Starting on ",cancer_type,sep=""))
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
  
  VCF_files_to_read <- as.vector(Cohort$vcf[Cohort$subtype %in% cancer_type])
  VCF_files_ID      <- as.vector(Cohort$sample[Cohort$subtype %in% cancer_type])
  
  #Read VCF as Grange object
  print("reading file")
  cancer_vcfs <- read_vcfs_as_granges(VCF_files_to_read,VCF_files_ID,ref_genome) #Gets warning messages about alternative alleles
  
  
  #Create mutational matrix 
  print("creating mutational matrix")
  mutational_matrix = mut_matrix(cancer_vcfs,ref_genome)
  #Write only TCGA Patient ID
  colnames(mutational_matrix) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",(gsub("TCGA-[A-Z0-9]*-","",colnames(mutational_matrix))))
  
  #Optinal find best rank
  if (FALSE){
    
    library(NMF)
    #Add pseudocount to mut matrix ?
    #mut_mat <- mut_mat + 0.0001
    estim.r <- nmf(mutational_matrix, rank=1:8,method="brunet", nrun=20, seed=123456)
    plot(estim.r)
    
  }
  
  
  
  
  #Extract mutational signatures (Cant extract more signatures than given samples)
  #Rank set to 6 , see rankoptimisation ACC
  number_of_signatures <- 6 
  
  print("Extracting signatures")
  
  extracted_sign = extract_signatures(mutational_matrix, number_of_signatures) 
  colnames(extracted_sign$signatures) <- c(1:number_of_signatures)
  rownames(extracted_sign$contribution) <- c(1:number_of_signatures)
  
  
  #Create folder for saving plots and set wd
  directory_name <- paste("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Pictures","/",cancer_type,sep="") 
  
  #crate folder if it doesnt exist
  ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
  setwd(directory_name)
  
  
  pdf_name <- paste(cancer_type,"_all_plots.pdf",sep = "")
  pdf(pdf_name)
  
  #Plot extracted signatures
  print(plot_96_profile(extracted_sign$signatures, condensed = TRUE))
  
  
  

  #Compare extracted signatures with cosmic signatures
  cosine_similarity = cos_sim_matrix(extracted_sign$signatures, cosmic_signatures)
  # Plot heatmap with specified signature order
  hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
  # store signatures in new order
  cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
  
  print(plot_cosine_heatmap(cosine_similarity, col_order = cosmic_order, cluster_rows = TRUE))
 
  
  
  #Plot contribution of extracted signature vs cancer
  print(plot_contribution(extracted_sign$contribution, extracted_sign$signature, mode = "absolute", coord_flip = TRUE))
 
  
  #Plot contribution of cosmic vs cancer
  
  fit_res <- fit_to_signatures(mutational_matrix, cosmic_signatures)
  # Select signatures with some contribution
  select <- which(rowSums(fit_res$contribution) > 10)
 
  
  # Plot contribution barplot
  

  print(plot_contribution(fit_res$contribution[select,],
                    cosmic_signatures[,select],
                    coord_flip = TRUE,
                    mode = "absolute"))
  

  
  
  
  #Plot cosine mutational_matrix vs cosmic_signature
  
  #Similarity between cosmic_signatures and mut_matrix
  cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
 
  
  
  print(plot_cosine_heatmap(cos_sim_samples_cosmic,
                      col_order = cosmic_order,
                      cluster_rows = TRUE))
  
  
  
  #Plot cosine mutational_matrix vs extracted_signatures
  

  print(plot_contribution_heatmap(extracted_sign$contribution, cluster_samples=TRUE))
  
  
  
  dev.off()
  
  
  #Calculate number of mutations for each sample.
  #Does not get the same as the cohort$mutations -> why? 
  
  Number_of_mutation <- cancer_vcfs@unlistData@ranges@NAMES
  Number_of_mutation <- count(Number_of_mutation)
  colnames(Number_of_mutation) <- c("Sample_id","#mutations")
  table_name <- paste(cancer_type,"_number_of_mutations.txt")
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/#mutations")
  write.table(Number_of_mutation,table_name,row.names = FALSE)
  
  
  
  break() #turn on if you want to try only one sample
  
}



