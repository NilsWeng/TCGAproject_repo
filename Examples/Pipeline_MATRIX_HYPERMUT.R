#Pipeline for plotting hyperMutaded samples and deleted genes as matrix

rm(list=ls())
gc()
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
#file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)

file_list <- list.files(pattern="all_thresholded.by_genes",recursive = TRUE)





GISTIC <- do.call("cbind",lapply(file_list ,read_GISTIC))
GISTIC <- GISTIC[, !duplicated(colnames(GISTIC))]

#select on samples present in high_mut_samples
keep <- c(High_mut_samples,c("Gene Symbol","Locus ID","Cytoband"))
  
colnames(GISTIC) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(GISTIC))


GISTIC <- GISTIC[, colnames(GISTIC) %in% keep]



# Get MC3 data for samples/genes ------------------------------------------
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")

#MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz')
load("MC3_DF.rda")

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



#Common samples
common_samples <- intersect(MC3_DF$Tumor_Sample_Barcode,colnames(GISTIC))

MC3_DF <- MC3_DF[MC3_DF$Tumor_Sample_Barcode %in% common_samples ,]
GISTIC <- GISTIC %>% select("Gene Symbol","Locus ID","Cytoband",common_samples)

#Common genes
common_genes <- intersect(unique(MC3_DF$Hugo_Symbol),GISTIC$`Gene Symbol`)

MC3_DF <- MC3_DF[MC3_DF$Hugo_Symbol %in% common_genes ,]
GISTIC <- GISTIC[GISTIC$`Gene Symbol` %in% common_genes ,]

MC3_DF <- MC3_DF[order(MC3_DF$Hugo_Symbol) ,]
GISTIC <- GISTIC[order(GISTIC$`Gene Symbol`) ,]
  


# Create a matrix ---------------------------------------------------------

hypermut_matrix <-  matrix(nrow=length(common_genes),ncol = length(common_samples))
colnames(hypermut_matrix) <- common_samples[order(common_samples)]
rownames(hypermut_matrix) <- common_genes[order(common_genes)]
hypermut_matrix <- as.data.frame(hypermut_matrix)

  
for (sample in common_samples){
  
  
  
  fill_vector_CN <- GISTIC %>% select(`Gene Symbol`,sample)
  colnames(fill_vector_CN) <- c("gene","CN")
  
  #Several mutations in same gene
  fill_vector_type <- MC3_DF %>% filter(MC3_DF$Tumor_Sample_Barcode == sample) %>% select(Hugo_Symbol,Variant_Classification)
  
  
  # Just write several for those with many mutations in gene
  
  #Get genes with several mutations of one kind
  fill_vector_type <- unique(fill_vector_type)
  
  
  dup <-duplicated(fill_vector_type$Hugo_Symbol)
  dup <- fill_vector_type[dup ,1]
  fill_vector_type[fill_vector_type$Hugo_Symbol %in% dup,2] <- "Several"
  fill_vector_type <- unique(fill_vector_type)
  
  
  #To handle multiple mutation in gene
  #fill_vector_type <- fill_vector_type %>%
    #group_by(Hugo_Symbol) %>%
    #summarise(Variant_Classification = paste(Variant_Classification, collapse = ","))
  
  fill_vector_type <- as.data.frame(fill_vector_type)
  colnames(fill_vector_type) <- c("gene","Type")
  

 
  
  
  fill_vector <- left_join(fill_vector_CN,fill_vector_type,by="gene")
  fill_vector$fill<- paste(fill_vector$CN,fill_vector$Type,sep=";")
  hypermut_matrix[colnames(hypermut_matrix) == sample] <- fill_vector$fill
  
  
}


# Create 2 matrix: CNVs and SNVs ------------------------------------------
#For running Malins script

#CNVs
cnvs <- hypermut_matrix
cnvs[] <- lapply(cnvs, function(x) (gsub(";..*$", "", x)))
cnvs <- as.matrix(cnvs)
class(cnvs) <- "numeric"
cnvs   <- cnvs + 2


#SNVS
snvs <- hypermut_matrix
snvs[] <- lapply(snvs, function(x) (gsub("^.*;", "", x)))
snvs <- as.matrix(snvs)

#Deal with several mutations in gene
test <- strsplit(snvs,",")
test <- 




cnvs <- t(cnvs)

snvs <- t(snvs)


# Visualise/plot matrix ---------------------------------------------------


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/Pictures")
imagefile <- "TEST.pdf"
subtypes <- 




cnvcol=c("#D55E00", "#E69F00","#ffffff","#56B4E9","#0072B2", "#08306b")
maincol=c("#737373","#bdbdbd")


#pdf(imagefile, width=7, height=7)

x = 1:(nrow(cnvs))
y = 1:ncol(cnvs)
centers <- expand.grid(y,x)
layout(matrix(c(1,2), nrow=2, ncol=1),  widths=1, heights=c(9,1))
#CNV matrix 
par(mar = c(0,4.2,4.2,1))
image(x,y,cnvs,
      col = cnvcol,
      #breaks = c(0,1,2,3,4,5),
      xaxt = 'n', 
      yaxt = 'n', 
      xlab = '', 
      ylab = '',
      ylim = c(max(y) + 0.5, min(y) - 0.5)
)

#add black lines
abline(h=y + 0.5, col=maincol[2])
abline(v=x + 0.5, col=maincol[2])
box(col=maincol[1])

 
#Add SNVs

#All variants
#levels(factor(MC3_DF$Variant_Classification))

for (l in as.integer(levels(factor(snvs)))){
  if(l>0){
    snvs_marks=which(snvs == l, arr.ind=T)
    points(snvs_marks, pch = snvsymbols[l], bg=snvcol[l], col=snvcol[l], cex=0.7, lwd=0.7)
  }
}



 
dev.off()
  
  
 