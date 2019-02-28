#Script to read and combine gistic files to one large gistic file

rm(list=ls())
setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads/GISTIC")

library(data.table)
library(dplyr)

file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)

#file_list <- file_list[1:2]

#Read DNA_repair genes list
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data")
DNA_repair <- read.csv("DNA_repair.csv",sep=";")

setwd("C:/Users/Nils_/Downloads/Firebrowse/Downloads/GISTIC")

read_GISTIC <- function(files){
  
  
 
  Table <- fread(files,header=TRUE,stringsAsFactors = FALSE)
  #Filter only for genes involved in DNA-repair
  Table <- Table %>% filter(Table$`Gene Symbol` %in% DNA_repair$hgnc_symbol)
  return(Table)
  
  
}



dataset <- do.call("cbind",lapply(file_list ,read_GISTIC))
dataset <- dataset[, !duplicated(colnames(dataset))]



