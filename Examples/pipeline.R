#CMD
rm(list=ls())


library(nilsTCGA)



#Load required CNV-data from folder.Masked Copy Number Segment
setwd("C:/Users/Nils_/OneDrive/Skrivbord/R-stuff/Data")
cnv_files <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/R-stuff/Data",recursive=TRUE,pattern="grch38.seg.v2.txt")

my_data <- lapply(cnv_files,readCNVs)
for (file in cnv_files){
  
  cnv_files(file)
  
}

# find all positions in all samples where segment_mean < -2
potential_loss <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("GDC_Aliquot","Chromosome","Start","End","Num_Probes","Segment_Mean"))


# IDEA about naming the GRanges objects after their id. However what is it