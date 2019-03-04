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


 hypermut_matrix <- hypermut_matrix[1:20 , 1:20]

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



# Catagorize snvs ---------------------------------------------------------

# If there are to many types bundle together unessesary into other





# Visualise/plot matrix ---------------------------------------------------


image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, axis.pos=1, add.axis=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
  if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
  box()
  if(add.axis) {axis(axis.pos)}
}





library(randomcoloR)

#Set colour and labels 

#SNV
snv_labels <- unique(levels(factor(snvs)))
snvsymbols <- c(1:length(snv_labels))
snvcol <- distinctColorPalette(length(snv_labels))

#CNV
cnv_labels <- unique(levels(factor(cnvs)))
#cnvcol <- distinctColorPalette(length(cnv__labels))
cnvcol=c("#D55E00", "#E69F00","#ffffff","#56B4E9","#0072B2")

#Other
labelcolors <- c("#D55E00", "#0072B2")
maincol=c("#737373","#bdbdbd")

#Cancertype


get_cancertype <- function(name){
  
  cancertype <- cohort[grep(name,cohort$sample), 2]
  
  return(cancertype)
}


cancertype <- rownames(cnvs)
cancertype <- as.vector(sapply(cancertype,get_cancertype))
cancertype_labels <- levels(factor(cancertype))
cancertype_pos <- vector()
cancertype_col <- distinctColorPalette(length(cancertype_labels))

#position




#Labels
samplenames <- rownames(cnvs)
samplenames <- gsub("TCGA-","",samplenames)
samplenames <- gsub("-[A-Z0-9]*$","",samplenames)
samplecol <- character()
genenames <- colnames(cnvs)






setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/Pictures")
imagefile <- "TEST.pdf"





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
for (i in 1:length(snv_labels)){
  
  if (snv_labels[i] == "NA"){
    
    next()
    
  }
  
  snvs_marks = which(snvs == snv_labels[i],arr.ind=T)
  points(snvs_marks, pch = snvsymbols[i], bg=snvcol[i], col=snvcol[i], cex=0.7, lwd=0.7)
  
}






#Gene names
par(mar = c(0,4,4,1))
mtext(genenames, at=(1:ncol(cnvs)), side=2, adj=1, las=1, cex=0.6, font=2, col=maincol[1])

#Sample names
par(mar = c(0,4.2,3.3,1))
mtext(samplenames, at=1:nrow(cnvs), side=3, cex=0.5, las=3, adj=0, font=2)
 
#
par(mar = c(0,4.2,3.8,1))
mtext(cancertype, at=1:nrow(cnvs), side=3, cex=0.6, las=1, padj=1, adj=0.5, col=labelcolors[1], font=2)


#



par(mar=c(2,4.2,1,10))
image.scale(cnvs, col=cnvcol,  breaks =c(1:(length(cnv_labels)+1)), add.axis=FALSE)
abline(v=c(0,1,2,3,4,5))
mtext(c(0:(length(cnv_labels)-1)), at=((2:(length(cnv_labels)+1))-0.5), padj = 0.5, adj=0.5, side=1, cex=0.7, col=maincol[1])
#mtext(c(4,0,1,2,3), at=((1:(length(cnv_labels)+1))-0.5), padj = 0.5, adj=0.5, side=1, cex=0.7, col=maincol[1]) # so stupid
box(col=maincol[1])


#SNV explanation
####FRÅGA MALIN 



 