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

cancer_list <- cancer_list[cancer_list %in% c("UCS","UVM") ]

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







# Total/Extract signatures from highly mutated samples at the same time --------------
rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(plyr)
library(dplyr)
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")

vcf_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern=".vcf")

#Remove cohort.vcf files
vcf_list <- vcf_list[grep("cohort",vcf_list,invert = TRUE)]


vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))


#Get list with highly mutated samples
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cohort <- read.table("cohort.txt",stringsAsFactors = FALSE)

treshold_nmut <- function(cancer_name){
  
  samples <- cohort %>% filter(subtype == cancer_name)
  
  median <- median(samples$nmut)
  S <- sd(samples$nmut)
  treshold <- median + 3*S
  
  high_mut_sample <- samples %>% filter(nmut > treshold) %>% select(sample)
  
  return(high_mut_sample)
  
  
}

#Extract vector of names containing highly mutated samples
High_mut_samples <- sapply(unique(cohort$subtype), treshold_nmut)
High_mut_samples <- unlist(High_mut_samples)
High_mut_samples <- unname(High_mut_samples)






#Trim vcf_list for only high mutation samples

vcf_list <- vcf_list[vcf_list_names %in% High_mut_samples]
vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))


#--------------------------------------------------------------------


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
#Read vcf-files and save as R object.
start_time <- Sys.time()

print("Reading all vcfs")

all_vcfs <- read_vcfs_as_granges(vcf_list,vcf_list_names,ref_genome)

print("Done reading all vcfs!!!!!!!")


end_time <- Sys.time()
print("Read time")
end_time - start_time # Took about 10 min for 65 samples


#Save the all_vcfs object
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations(3s)")
save(all_vcfs,file="high_mut_vcfs.rda")


#Extract mutational matrix
print("creating mutational matrix")
mutational_matrix = mut_matrix(all_vcfs,ref_genome)
save(mutational_matrix,file="high_mut_matrix.rda")

#Add cancer abb as name
colnames(mutational_matrix) <- gsub("/.*$","",vcf_list) ####################For having cancer abb as names
#[cohort$sample %in% vcf_list_names]

#just to give unique name to each (so one can plot)
for (i in 1:length(colnames(mutational_matrix))) {
  
  colnames(mutational_matrix)[i]<- paste(colnames(mutational_matrix)[i],i,sep="-")
  
}


#Find optimal rank
library(NMF)
#Add pseudocount to mut matrix ?
#mut_mat <- mut_mat + 0.0001
estim.r <- nmf(mutational_matrix, rank=1:8,method="brunet", nrun=12, seed=123456)
plot(estim.r)



#Clear workspace 
#rm(all_vcfs,vcf_list,vcf_list_names)
#Extract signatures
number_of_signatures <- 2
print("Extracting signatures")

extracted_sign = extract_signatures(mutational_matrix, number_of_signatures) 

colnames(extracted_sign$signatures) <- c(1:number_of_signatures)
rownames(extracted_sign$contribution) <- c(1:number_of_signatures)

save(extracted_sign, file="extracted_sign_high_mut.rda")
print("Made it past extracting")


#colnames(extracted_sign$contribution) <- colnames(mutational_matrix) ####################For having cancer abb as name


#Needed for comparing to cosmic
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
#cosmic_signatures <- as.matrix(read.table("cosmic_signatures.txt",header=TRUE))
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))


#Print results

#Set options here
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations")
pdf_name <- paste("HighMut_new_cosmic.pdf")


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
  
  
  
  #plot reconstructed mutational profiles from extracted
  print(plot_compare_profiles(mutational_matrix[,1],
                          extracted_sign$reconstructed[,1],profile_names = c("Original", "Reconstructed"),condensed = TRUE))
  
  
  
  #plot reconstructed mutational profiles from cosmic
  print(plot_compare_profiles(mutational_matrix[,1], fit_res$reconstructed[,1],
                        profile_names = c("Original", "Reconstructed(cosmic)"),condensed = TRUE))
  
  
  
  dev.off()
}

print_mutsign()




# Cluster cosine-sim matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
plot(sample_cluster)
N<-5
method <- "complete"
cluster_groups <- cutree(sample_cluster,k=N)

#Plot with found clusters
plot(sample_cluster, cex = 0.6)
rect.hclust(sample_cluster, k = 5, border = 2:5)


#Determine optimal number of clusters (3 different methods)
library(factoextra)
#elbow method
fviz_nbclust(cos_sim_samples_cosmic, FUN = hcut, method = "wss")
#silhouette method
fviz_nbclust(cos_sim_samples_cosmic, FUN = hcut, method = "silhouette")
#Gap statistic method
set.seed(123)
gap_stat <- clusGap(cos_sim_samples_cosmic, FUN = kmeans, nstart = 25,K.max = 15, B = 50)
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

#tidyverse seems like a nice package !!!!

#Extract what samples are in each cluster

cluster_groups <- data.frame("Sample_id"=vcf_list_names,"Cluster"=as.vector(cluster_groups))
filename <- paste("Sample_In_Cluster",method,N,sep="_")
write.table(cluster_groups,paste(filename,".txt",sep=""),row.names = FALSE)
