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

cancer_list <- cancer_list[cancer_list %in% "ACC" ]

#Alter what cancer types you want to loop over
#not_in = head(cancer_list,29)
#cancer_list <- cancer_list[!(cancer_list %in% not_in)]



for (cancer_type in cancer_list){
  
  
   
  print(paste("Starting on ",cancer_type,sep=""))
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
  
  VCF_files_to_read <- as.vector(Cohort$vcf[Cohort$subtype %in% cancer_type])
  VCF_files_ID      <- as.vector(Cohort$sample[Cohort$subtype %in% cancer_type])
  
  #Read VCF as Grange object
  print("reading file")
  cancer_vcfs <- read_vcfs_as_granges(VCF_files_to_read,VCF_files_ID,ref_genome) #Gets warning messages about alternative alleles
  
  break()
  #Create mutational matrix 
  print("creating mutational matrix")
  mutational_matrix = mut_matrix(cancer_vcfs,ref_genome)
  #Write only TCGA Patient ID
  colnames(mutational_matrix) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",(gsub("TCGA-[A-Z0-9]*-","",colnames(mutational_matrix))))
  colnames(mutational_matrix)
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
  print("Made it past extracting")
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
 
  

  #Ugly plot for studies with many patients
  k_vector <- c(1:33)
  k_vector <- c(1:length(VCF_files_to_read))
  
  while (length(k_vector) > 0) {
    
    interval_to_plot <- head(k_vector,100)
    print(plot_contribution(extracted_sign$contribution, extracted_sign$signature,
                      mode = "absolute", coord_flip = TRUE,index = interval_to_plot))
    
    
    k_vector <- k_vector[! (k_vector %in% interval_to_plot) ]
    
  }
  
  
  
  
  #Plot contribution of extracted signature vs cancer
  #print(plot_contribution(extracted_sign$contribution, extracted_sign$signature, mode = "absolute", coord_flip = TRUE))
 
  
  #Plot contribution of cosmic vs cancer
  
  fit_res <- fit_to_signatures(mutational_matrix, cosmic_signatures)
  # Select signatures with some contribution
  select <- which(rowSums(fit_res$contribution) > 10)
 
  
  # Plot contribution barplot
  
  #Ugly plot for studies with many patients
  
 
  k_vector <- c(1:length(VCF_files_to_read))
  
  while (length(k_vector) > 0) {
    
    interval_to_plot <- head(k_vector,100)
    print(plot_contribution(fit_res$contribution[select,],
                            cosmic_signatures[,select],
                            coord_flip = TRUE,
                            mode = "absolute",
                            index = interval_to_plot))
    
    
    k_vector <- k_vector[! (k_vector %in% interval_to_plot) ]
    
  }
  


  
  
  
  #Plot cosine mutational_matrix vs cosmic_signature
  
  #Similarity between cosmic_signatures and mut_matrix
  cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
 
  #error NA/NaN/Inf in foreign function call #2835-03B
  #cos_sim_s amples_cosmic[which(cos_sim_samples_cosmic == 0.0000000000)] <- 0.0000000001
  
  
  
  #######CRAAAAAAASSSSSSSSSSSH here
  print(plot_cosine_heatmap(cos_sim_samples_cosmic,
                      col_order = cosmic_order,
                      cluster_rows = TRUE))
  
  
  
  #Plot cosine mutational_matrix vs extracted_signatures
  

  print(plot_contribution_heatmap(extracted_sign$contribution, cluster_samples=TRUE))
  
  
  
  
  
  
  
  #Save-Data as R object
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/generated_data")
  
  #Create folder for saving data and set wd
  directory_name <- paste("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/generated_data","/",cancer_type,sep="") 
  
  #crate folder if it doesnt exist
  ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
  setwd(directory_name)
  
  save(mutational_matrix,file="mut_mat.rda")
  save(extracted_sign, file=paste("extracted_sign_",number_of_signatures,".rda",sep = ""))
  
  #Calculate number of mutations for each sample.
  #Does not get the same as the cohort$mutations -> why? 
  Number_of_mutation <- cancer_vcfs@unlistData@ranges@NAMES
  Number_of_mutation <- count(Number_of_mutation)
  colnames(Number_of_mutation) <- c("Sample_id","#mutations")
  table_name <- paste(cancer_type,"_number_of_mutations.txt")
  write.table(Number_of_mutation,table_name,row.names = FALSE)
  
  #Plot mutation distribution
  
  
  d <- density(Number_of_mutation$`#mutations`)
  
  print(plot(d, type="n", main="Distribution of mutations")
        ,polygon(d, col="lightgray", border="gray")
        ,rug(Number_of_mutation$`#mutations`, col="red"))
  
  
  

  
  
  
  dev.off()
  
  
  #break() #turn on if you want to try only one sample
  
}



# Total/Extract signatures from all samples at the same time --------------
rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(plyr)
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")

vcf_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern=".vcf")

#Unselect cohort.vcf files
vcf_list <- vcf_list[grep("cohort",vcf_list,invert = TRUE)]
#vcf_list <- vcf_list[grep("cohort",vcf_list)]

vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))



#############-------------------------------------------------------------------------------
#Select subgroup with high mutation-number.
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/generated_data")
mut_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/generated_data",recursive=TRUE,pattern="number_of_mutations.txt")

#Open each mutation list file and select samplenames with certain number of mutations.

number_mut_DF <- data.frame()

for (mut_file in mut_list) {
  
  
  file_df <- read.table(mut_file,header=TRUE)
  file_df$study_abb <- gsub("/..*","",mut_file)
  
  number_mut_DF <-rbind(number_mut_DF,file_df)
  
  
}


number_mut_DF <- number_mut_DF[number_mut_DF$X.mutations >= 2000 ,]


#Trim vcf_list for only high mutation samples

vcf_list <- vcf_list[vcf_list_names %in% number_mut_DF$Sample_id]
vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")

#################################--------------------------------------------------------------------



#Read vcf-files and save as R object.
start_time <- Sys.time()

print("Reading all vcfs")

all_vcfs <- read_vcfs_as_granges(vcf_list,vcf_list_names,ref_genome)

print("Done reading all vcfs!!!!!!!")


end_time <- Sys.time()
print("Read time")
end_time - start_time


#Save the all_vcfs object
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign")
save(all_vcfs,file="all_vcfs.rda")


#Extract mutational matrix
print("creating mutational matrix")
mutational_matrix = mut_matrix(all_vcfs,ref_genome)
save(mutational_matrix,file="all_mut_matrix.rda")

colnames(mutational_matrix) <- number_mut_DF$study_abb ####################For having cancer abb as names
for (i in 1:length(colnames(mutational_matrix))) {
  
  colnames(mutational_matrix)[i]<- paste(colnames(mutational_matrix)[i],i,sep="-")
  
}


#Find optimal rank
library(NMF)
#Add pseudocount to mut matrix ?
#mut_mat <- mut_mat + 0.0001
estim.r <- nmf(mutational_matrix, rank=1:10,method="brunet", nrun=10, seed=123456)
plot(estim.r)



#Clear workspace 
#rm(all_vcfs,vcf_list,vcf_list_names)
#Extract signatures
number_of_signatures <- 4
print("Extracting signatures")

extracted_sign = extract_signatures(mutational_matrix, number_of_signatures) 

colnames(extracted_sign$signatures) <- c(1:number_of_signatures)
rownames(extracted_sign$contribution) <- c(1:number_of_signatures)

save(extracted_sign, file="extracted_sign_high_mut.rda")
print("Made it past extracting")


colnames(extracted_sign$contribution) <- colnames(mutational_matrix) ####################For having cancer abb as name


#Needed for comparing to cosmic
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cosmic_signatures <- as.matrix(read.table("cosmic_signatures.txt",header=TRUE))
#Print results

#Set options here
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations")
pdf_name <- paste("High_mutation_samples.pdf")


print_mutsign <- function() {
  
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
  
  
  
  #Ugly plot for studies with many patients
  k_vector <- c(1:length(vcf_list))
  #k_vector <- c(1:length(VCF_files_to_read))
  
  while (length(k_vector) > 0) {
    
    interval_to_plot <- head(k_vector,100)
    print(plot_contribution(extracted_sign$contribution, extracted_sign$signature,
                            mode = "absolute", coord_flip = TRUE,index = interval_to_plot))
    
    
    k_vector <- k_vector[! (k_vector %in% interval_to_plot) ]
    
  }
  
  
  
  
  #Plot contribution of extracted signature vs cancer
  #print(plot_contribution(extracted_sign$contribution, extracted_sign$signature, mode = "absolute", coord_flip = TRUE))
  
  
  #Plot contribution of cosmic vs cancer
  
  fit_res <- fit_to_signatures(mutational_matrix, cosmic_signatures)
  # Select signatures with some contribution
  select <- which(rowSums(fit_res$contribution) > 10)
  
  
  # Plot contribution barplot
  
  #Ugly plot for studies with many patients
  
  
  k_vector <- c(1:length(vcf_list))
  
  while (length(k_vector) > 0) {
    
    interval_to_plot <- head(k_vector,100)
    print(plot_contribution(fit_res$contribution[select,],
                            cosmic_signatures[,select],
                            coord_flip = TRUE,
                            mode = "absolute",
                            index = interval_to_plot))
    
    
    k_vector <- k_vector[! (k_vector %in% interval_to_plot) ]
    
  }
  
  
  
  #Plot cosine mutational_matrix vs cosmic_signature
  
  #Similarity between cosmic_signatures and mut_matrix
  cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
  
  #error NA/NaN/Inf in foreign function call #2835-03B
  #cos_sim_s amples_cosmic[which(cos_sim_samples_cosmic == 0.0000000000)] <- 0.0000000001
  
  
  
  #######CRAAAAAAASSSSSSSSSSSH here
  print(plot_cosine_heatmap(cos_sim_samples_cosmic,
                            col_order = cosmic_order,
                            cluster_rows = TRUE))
  
  
  #Plot cosine mutational_matrix vs extracted_signatures
  print(plot_contribution_heatmap(extracted_sign$contribution, cluster_samples=TRUE))
  
  
  
  dev.off()
}

print_mutsign()

























#Since reading the vcf file takes long time 
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


MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz') #<- run for entire file

end_time <- Sys.time()
print("Read time")
end_time - start_time


#Select only necessary columns to save memory
library(dplyr)
MC3_DF <- select(MC3_DF, Chromosome,Tumor_Seq_Allele2,Tumor_Sample_Barcode,Start_Position,Reference_Allele)




#Add Study.abbrevation to the MC3_DF matrix.
MC3_DF$TSS.Code <-gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",gsub("TCGA-","",MC3_DF$Tumor_Sample_Barcode))
MC3_DF <- join(MC3_DF,TSS2Study_DF, by="TSS.Code")


#Split into samples

MC3_DF_samples <- split (MC3_DF, f=MC3_DF$Tumor_Sample_Barcode) # mayby i need drop = TRUE ?






# Construct into VCF format
#sample_id <- unique(Sample$Tumor_Sample_Barcode);

maf2vcf <- function(maf_item){
  
  vcfdata <- matrix(".", nrow=nrow(maf_item), ncol = 10);
  columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT","Sample_ID" )
  colnames(vcfdata) <- columns
  
  vcfdata[,1]=as.character(maf_item$Chromosome)
  vcfdata[,2]=as.character(maf_item$Start_Position)
  vcfdata[,4]=as.character(maf_item$Reference_Allele)
  vcfdata[,5]=as.character(maf_item$Tumor_Seq_Allele2)
  vcfdata[,9]="GT"
  #vcfdata[,10]="1/0"
  vcfdata[,10] = as.character(maf_item$Tumor_Sample_Barcode)
  
  return(vcfdata)
  
}


MC3_vcf <- lapply(MC3_DF_samples, maf2vcf)


# Convert to GRange

vcf2Grange <- function(VCF_item) {
  
  ref_genome <- base::get(genome)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  
  
  
  
  #GRanges(seqnames = VCF_item$Chromosome,REF )
  
  
  
  
  
}


