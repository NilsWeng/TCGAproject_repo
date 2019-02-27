#Script for discovering hypermutaded samples
rm(list=ls())

#Import needed packages
#library(plyr)
library(dplyr)
#Set reference genome


#Load cosmic_signatures from file in folder
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")

cohort <- read.table("cohort.txt",stringsAsFactors = FALSE)





treshold_nmut <- function(cancer_name){
  
  samples <- cohort %>% filter(subtype == cancer_name)
  
  median <- median(samples$nmut)
  S <- sd(samples$nmut)
  treshold <- median + 10*S
  
  high_mut_sample <- samples %>% filter(nmut > treshold) %>% select(sample)
  
  return(high_mut_sample)
  
  
}
  
  HyperMut_vector <- sapply(unique(cohort$subtype), treshold_nmut)
  
  HyperMut_vector <- unlist(HyperMut_vector)
  HyperMut_vector <- unname(HyperMut_vector)


  cancer_type <- cohort %>% filter(sample %in% HyperMut_vector) %>% select(subtype)
  cancer_type <- table(cancer_type)
  cancer_type <- as.data.frame(cancer_type)
  colnames(cancer_type) <- c("name","freq")
  
 
  
  library(ggplot2)
  
  library(randomcoloR)
  n <- nrow(cancer_type)
  palette <- distinctColorPalette(n)
  
  
  bp<- ggplot(cancer_type, aes(x="", y=freq, fill=name))+geom_bar(width = 1, stat = "identity")
  pie <- bp + coord_polar("y", start=0) + scale_fill_manual(values=palette) + scale_y_continuous(labels=freq)
  pie

  
  
 pie(cancer_type$freq, labels = cancer_type$freq, col = palette,fill = palette)
  
  #ggplot(cancer_type, aes(x = freq, fill = name)) + geom_bar(width = 1)  + coord_polar(theta = "y") + theme_void()
