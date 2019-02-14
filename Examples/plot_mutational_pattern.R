#Clear workspace
rm(list=ls())

#Import needed packages
library(MutationalPatterns)

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

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



#Load cosmic_signatures from file in folder
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cosmic_signatures <- as.matrix(read.table("cosmic_signatures.txt",header=TRUE))


#Compare extracted signatures with cosmic signatures
cosine_similarity = cos_sim_matrix(cohort_extracted_signatures$signatures, cosmic_signatures)
# Plot heatmap with specified signature order
hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]





plot_cosine_heatmap(cosine_similarity, col_order = cosmic_order, cluster_rows = TRUE)













