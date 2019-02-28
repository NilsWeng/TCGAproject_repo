#Pipeline for plotting hyperMutaded samples and deleted genes as matrix

rm(list=ls())

#library
library(data.table)
library(dplyr)




# Get list of higly mutated samples ---------------------------------------

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cohort <- read.table("cohort.txt",stringsAsFactors = FALSE)

treshold_nmut <- function(cancer_name){
  
  samples <- cohort %>% filter(subtype == cancer_name)
  
  median <- median(samples$nmut)
  S <- sd(samples$nmut)
  treshold <- median + 5*S
  
  high_mut_sample <- samples %>% filter(nmut > treshold) %>% select(sample)
  
  return(high_mut_sample)
  
  
}
  
#Extract vector of names containing highly mutated samples
High_mut_samples <- sapply(unique(cohort$subtype), treshold_nmut)
High_mut_samples <- unlist(High_mut_samples)
High_mut_samples <- unname(High_mut_samples)
High_mut_samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",High_mut_samples)


# Get list of DNA-repair genes --------------------------------------------
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")
DNA_repair <- read.csv("DNA_repair.csv",sep=";")




# Get CopyNumber data for samples/genes -----------------------------------

read_GISTIC <- function(files){
  
  
  
  Table <- fread(files,header=TRUE,stringsAsFactors = FALSE)
  #Filter only for genes involved in DNA-repair
  Table <- Table %>% filter(Table$`Gene Symbol` %in% DNA_repair$hgnc_symbol)
  

  return(Table)
  
  
}

setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads/GISTIC")
file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)



GISTIC <- do.call("cbind",lapply(file_list ,read_GISTIC))
GISTIC <- GISTIC[, !duplicated(colnames(GISTIC))]

#select on samples present in high_mut_samples
keep <- c(High_mut_samples,c("Gene Symbol","Locus ID","Cytoband"))
  
colnames(GISTIC) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(GISTIC))


GISTIC <- GISTIC[, colnames(GISTIC) %in% keep]



# Get MC3 data for samples/genes ------------------------------------------
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")

MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz')

MC3_DF <- MC3_DF %>% select("Hugo_Symbol","Chromosome","Start_Position"
                            ,"End_Position","Tumor_Sample_Barcode","Variant_Classification")

#Filter on genes
MC3_DF <- MC3_DF[MC3_DF$Hugo_Symbol %in% DNA_repair$hgnc_symbol ,]

#Modify sample-id
MC3_DF$Tumor_Sample_Barcode <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","", MC3_DF$Tumor_Sample_Barcode)

#Filter on samples
MC3_DF <- MC3_DF[MC3_DF$Tumor_Sample_Barcode %in% High_mut_samples ,]




# Compare MC3 to GISIC keep only common genes/samples ---------------------
#intersect(MC3_DF)
common_samples <- intersect(MC3_DF$Tumor_Sample_Barcode,colnames(GISTIC))
MC3_DF1 <- MC3_DF[MC3_DF$Tumor_Sample_Barcode %in% common_samples ,]



#Should be some samples removed doing this ????????---------------------------------------------
GISTIC1 <- GISTIC %>% select("Gene Symbol","Locus ID","Cytoband",common_samples)



#setdiff(MC3_DF$Tumor_Sample_Barcode,colnames(GISTIC))

# Create a matrix ---------------------------------------------------------





# Visualise/plot matrix ---------------------------------------------------


