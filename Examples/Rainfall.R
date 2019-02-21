#Rainfall plott

rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(plyr)
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

vcf_list <- list.files(path="C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files",recursive=TRUE,pattern="cohort.vcf")
#vcf_list_names <- gsub(".cohort.vcf","",gsub("[A-Z0-9]*/","",vcf_list))

#DF_vcf_list <- data.frame(vcf_list,vcf_list_names)
#DF_vcf_list <- DF_vcf_list[2:nrow(DF_vcf_list) ,]


vcf_list <- vcf_list[2:length(vcf_list)]
#vcf_list_names <- vcf_list_names[2:length(vcf_list_names)]
#vcf_list <- vcf_list[vcf_list_names %in% c("ACC","CHOL","U" )]
#vcf_list_names <- gsub(".cohort.vcf","",gsub("[A-Z0-9]*/","",vcf_list))


vcf_list <- "CHOL/CHOL.cohort.vcf"

for (cancer in vcf_list){
  
  cancer_name <- gsub(".cohort.vcf","",gsub("[A-Z0-9]*/","",cancer))
  
  
  
  print(paste("Starting on",cancer_name,sep = " "))
  

  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
  
  print("Reading VCFS")
  
  cancer_vcfs <- read_vcfs_as_granges(cancer,cancer_name,ref_genome)
  
  print("Rainfall plotting")
  

  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Rainfall")
  pdf_name  <- paste(cancer_name,"_rainfall.pdf",sep="")
  pdf(pdf_name)
  chromosomes <- seqnames(get(ref_genome))[1:22]
  
  print(plot_rainfall(cancer_vcfs[[1]], title = names(cancer_vcfs[1]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+09))
  
  
  dev.off()
  
}







# One sample at a time with, hold on --------------------------------------


# Test plotting one sample at a time
#Rainfall plott

rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(plyr)
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")

cohort <- read.table("cohort.txt",header = TRUE)

cancer_types <- unique(cohort$subtype)

library(ggplot2)
rainfall_plot<- function(vcf,title,cex,cex_text,ylim,chromosomes) {
  
  
  # get chromosome lengths of reference genome
  chr_length = seqlengths(vcf)
  table <- read.table("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/CHR_length_hg19.csv",sep = ";",header=TRUE)
  chr_length <- table$Lenght
  names(chr_length) <- table$Chromosome
  
  # subset
  chr_length = chr_length[names(chr_length) %in% chromosomes]
  
  # cumulative sum of chromosome lengths
  chr_cum = c(0, cumsum(as.numeric(chr_length)))
  
  # Plot chromosome labels without "chr"
  names(chr_cum) = names(chr_length)
  labels = gsub("chr", "", names(chr_length))
  
  
  # position of chromosome labels
  m=c()
  for(i in 2:length(chr_cum))
    m = c(m,(chr_cum[i-1] + chr_cum[i]) / 2)
  
  
  # mutation characteristics
  type = loc = dist = chrom = c()
  
  # for each chromosome
  for(i in 1:length(chromosomes))
  {
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    n = length(chr_subset)
    if(n<=1){next}
    type = c(type, mut_type(chr_subset)[-1])
    loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    chrom = c(chrom, rep(chromosomes[i],n-1))
  }
  
  data = data.frame(type = type,
                    location = loc,
                    distance = dist,
                    chromosome = chrom)
  
  # Removes colors based on missing mutation types.  This prevents colors from
  # shifting when comparing samples with low mutation counts.
  #typesin = SUBSTITUTIONS %in% levels(data$type)
  #colors = colors[typesin]
  
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  location = NULL
  
  # make rainfall plot
  plot = ggplot(data, aes(x=location, y=distance)) +
    geom_point(aes(colour=factor(type)), cex=cex) + 
    geom_vline(xintercept = as.vector(chr_cum), linetype="dotted") +
    annotate("text", x = m, y = ylim, label = labels, cex=cex_text) +
    xlab("Genomic Location") +
    ylab("Genomic Distance") +
    scale_y_log10() +
    scale_x_continuous(expand = c(0,0), limits=c(0, max(chr_cum))) +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()) + 
    guides(colour = guide_legend(nrow = 1))
  
  
  
  return(plot)
  
  
}



for (cancer in cancer_types){
  

  #setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF")
  
  samples_in_cancer <- as.vector(cohort[cohort$subtype %in% cancer ,]$vcf)
  
  
  chromosomes <- seqnames(get(ref_genome))[1:22]
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Rainfall/Test")
  
  pdf_name <- gsub(".vcf","",gsub("[A-Z0-9]*/","",cancer))
  pdf_name  <- paste(cancer,"_rainfall.pdf",sep="")
  
  pdf(pdf_name)
  
  
  
  
  for (sample in samples_in_cancer){
    
    
  
    
    
    sample_name <- gsub("/..*","",sample)
    setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
    
    
 
    
    cancer_vcfs <- read_vcfs_as_granges(sample,sample_name,ref_genome)
    

    
    #Plot rainfall
    
    setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Rainfall/Test")
    
    #pdf_name <- gsub(".vcf","",gsub("[A-Z0-9]*/","",sample))
    #pdf_name  <- paste(pdf_name,"_rainfall.pdf",sep="")
    
    #pdf(pdf_name)
    
    
    #print(rainfall_plot(cancer_vcfs[[1]],pdf_name,cex=1.5,cex_text = 3,ylim = 1e+08,chromosomes = chromosomes)) 
    
    # rainfall doesnt seem to get seqlenghts from cancer_vcfs[[1]] try to add it and see if it works
    #lengths <- seqlengths(cancer_vcfs[[1]])
    table <- read.table("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/CHR_length_hg19.csv",sep = ";",header=TRUE)
    seqlengths(cancer_vcfs[[1]]) <- table[table$Chromosome %in% names(seqlengths(cancer_vcfs[[1]])) ,][, 2]
    
    
    print(plot_rainfall(cancer_vcfs[[1]], title = names(cancer_vcfs[1]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+09))
    
   
   
    
    
  }
  
  dev.off()
  dev.off()
  break()
  
  


  
}




#input









# get chromosome lengths of reference genome
chr_length = seqlengths(vcf)
table <- read.table("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/CHR_length_hg19.csv",sep = ";",header=TRUE)
chr_length <- table$Lenght
names(chr_length) <- table$Chromosome

# subset
chr_length = chr_length[names(chr_length) %in% chromosomes]

# cumulative sum of chromosome lengths
chr_cum = c(0, cumsum(as.numeric(chr_length)))

# Plot chromosome labels without "chr"
names(chr_cum) = names(chr_length)
labels = gsub("chr", "", names(chr_length))


# position of chromosome labels
m=c()
for(i in 2:length(chr_cum))
  m = c(m,(chr_cum[i-1] + chr_cum[i]) / 2)


# mutation characteristics
type = loc = dist = chrom = c()

# for each chromosome
for(i in 1:length(chromosomes))
{
  chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
  n = length(chr_subset)
  if(n<=1){next}
  type = c(type, mut_type(chr_subset)[-1])
  loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
  dist = c(dist, diff(start(chr_subset)))
  chrom = c(chrom, rep(chromosomes[i],n-1))
}

data = data.frame(type = type,
                  location = loc,
                  distance = dist,
                  chromosome = chrom)

# Removes colors based on missing mutation types.  This prevents colors from
# shifting when comparing samples with low mutation counts.
typesin = "SUBSTITUTIONS" %in% levels(data$type)
colors = colors[typesin]

# These variables will be available at run-time, but not at compile-time.
# To avoid compiling trouble, we initialize them to NULL.
location = NULL

# make rainfall plot
plot = ggplot(data, aes(x=location, y=distance)) +
  geom_point(aes(colour=factor(type)), cex=cex) + 
  geom_vline(xintercept = as.vector(chr_cum), linetype="dotted") +
  annotate("text", x = m, y = ylim, label = labels, cex=cex_text) +
  xlab("Genomic Location") +
  ylab("Genomic Distance") +
  scale_y_log10() +
  scale_x_continuous(expand = c(0,0), limits=c(0, max(chr_cum))) +
  ggtitle(title) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) + 
  guides(colour = guide_legend(nrow = 1))





















