#Read and combine gene expression to large file
rm(list=ls())

library(data.table)
library(dplyr)

setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads/mRNASeq_mod")
file_list <- list.files(pattern="_normalized__data.data.txt",recursive = TRUE)


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")
DNA_repair <- read.csv("DNA_repair.csv",sep=";")




setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads/mRNASeq_mod")

read_mRNASeq <- function(files){
  
  
  
  Table <- fread(files,header=TRUE,stringsAsFactors = FALSE)
  #Filter only for genes involved in DNA-repair
  Table$`Hybridization REF` <- gsub("\\|.*$","",Table$`Hybridization REF`)
  Table <- Table %>% filter(Table$`Hybridization REF` %in% DNA_repair$hgnc_symbol)
  
  return(Table)
  
  
}

file_list <- file_list[1:3]

dataset <- do.call("cbind",lapply(file_list ,read_mRNASeq))
dataset <- dataset[, !duplicated(colnames(dataset))]

setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads")
write.table(dataset,"PANCAN_mRNA_expression.txt")


