#Kataegis


rm(list=ls())

#Import needed packages
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(plyr)

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
#Chromosomes length
chromosomes <- seqnames(get(ref_genome))[1:22]
table <- read.table("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/CHR_length_hg19.csv",sep = ";",header=TRUE)
chr_length <- table$Lenght
names(chr_length) <- table$Chromosome
chr_length = chr_length[names(chr_length) %in% chromosomes]

#culmulative length of chromosomes
chr_cum = c(0, cumsum(as.numeric(chr_length)))
names(chr_cum) = names(chr_length)
labels = gsub("chr", "", names(chr_length))

#Position for chromosome labels
m=c()
for(i in 2:length(chr_cum)){
  m = c(m,(chr_cum[i-1] + chr_cum[i]) / 2)
  
}







#Open cohort file
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cohort <- read.table("cohort.txt",header = TRUE)
cancer_types <- unique(cohort$subtype)




#Loop over each cancer
Kataegis_data <- data.frame()
for (cancer in cancer_types){
  
  
  samples_in_cancer <- as.vector(cohort[cohort$subtype %in% cancer ,]$vcf)
  
  for (sample in samples_in_cancer){
    
    sample_name <- gsub("/..*","",sample)
    setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/VCF_files")
    
    vcf <- read_vcfs_as_granges(sample,sample_name,ref_genome)
    vcf <- vcf[[1]]
    

    
    
    #Create a dataframe containing all relevant information for Kataegis
    type = loc = dist = chrom = id = c()
    
    for (i in 1:length(chromosomes)){
      
      chr_subset = vcf[seqnames(vcf) == chromosomes[i]]  
      n = length(chr_subset)
      sample_id <- unique(names(vcf))
      if(n<=1){next}
      type = c(type, mut_type(chr_subset)[-1])
      loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
      #Distance to previous mutation on same chromosome
      dist = c(dist, diff(start(chr_subset)))
      chrom = c(chrom, rep(chromosomes[i],n-1))
      id = c(id,rep(sample_id,n-1))
      
    
      
    }   
      
     
    sample_id <- gsub("[A-Z0-9]*/","",gsub(".vcf","",sample)) 
    
    
    
    data = data.frame(type = type,
                      location = loc,
                      distance = dist,
                      chromosome = chrom,
                      sample_id = id
                      )
    

    Kataegis_data <- rbind(Kataegis_data,data)
      
  
    
    
  }
  
  
  #Plot data
  
  #typesin = "SUBSTITUTIONS" %in% levels(data$type)
  #colors = colors(typesin)
  title <- "All at the same time - ACC"
  cex <- 1.5
  cex_text <- 3
  ylim <- 1e+09
  
  pdfname <- paste(cancer,"_Kataegis.pdf",sep = "")
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Rainfall/Test")
  pdf(pdfname)
 
  
  rainfall_plot <- function(Kataegis_data,title,cex,cex_text,ylim){
    
    
    ggplot(Kataegis_data, aes(x=location, y=distance)) +
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
    
  }
  
  #Plot all at the same time
  print(rainfall_plot(Kataegis_data,title,cex,cex_text,ylim))
  
  #Plot each sample
  for (item in split(Kataegis_data,Kataegis_data$sample_id,drop = TRUE) ){
    
    
    title <- as.vector(unique(item$sample_id))
    print(rainfall_plot(item,title,cex,cex_text,ylim))
  
  }
    

  
  
  
  dev.off()
  
    
  #Save kataegis object
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Rainfall/Test/generated_data")
  filename <- paste(cancer,"_Kataegis_table.rda",sep="")
  save(Kataegis_data,file=filename) 
  
  
  
  break()             #Remove break if you want to loop over entire dataset.
}
