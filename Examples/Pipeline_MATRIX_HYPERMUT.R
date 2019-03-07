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
  treshold <- median + 3*S
  
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




 #hypermut_matrix <- hypermut_matrix[1:60 ,]
 
if(FALSE){
  #Optional order data after cosine-mutsign clustering
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations(3s)")
  clustered_samples <- read.table("Sample_In_Cluster_complete_7.txt",header=TRUE,stringsAsFactors = FALSE)
  clustered_samples <- clustered_samples %>% arrange(Cluster) %>% select(Sample_id)
  
  hypermut_matrix1 <- left_join(clustered_samples,hypermut_matrix,by="sample")
  
  #sample_order = rownames(cnvs)[clustered_samples]
  
  
  
  
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
  
  
  get_cancertype <- function(name){
    
    cancertype <- cohort[grep(name,cohort$sample), 2]
    
    return(cancertype)
  }
  
  cancertype <- colnames(hypermut_matrix)
  cancertype <- as.vector(sapply(cancertype,get_cancertype))
  
  colnames(hypermut_matrix) <- gsub("TCGA-","",gsub("-[A-Z0-9]*$","",colnames(hypermut_matrix)))
  colnames(hypermut_matrix) <- paste(colnames(hypermut_matrix),cancertype,sep="-")
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





cnvs <- t(cnvs)

snvs <- t(snvs)



# Catagorize snvs ---------------------------------------------------------

# If there are to many types bundle together unessesary into other
non_AA_mod <- c("3'Flank","3'UTR","5'Flank","5'UTR","Intron","Silent")
frame_shift <- c("Frame_Shift_Del","Frame_Shift_Ins")
expression <- c("Splice_Site","Translation_Start_Site")
inframe    <-c("In_Frame_Del","In_Frame_Ins")
AA_mod     <- c("Frame_shift","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation")


snvs[snvs %in% non_AA_mod]  <- "Other"
snvs[snvs %in% frame_shift] <- "Frame_shift"
snvs[snvs %in% expression]  <- "Expression regulation"
snvs[snvs %in% inframe]     <- "In frame del/ins"
snvs[snvs %in% AA_mod]      <- "AA modifying mut"

unique(levels(factor(snvs)))

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



# GGplot section ----------------------------------------------------------


#Prova GG plot istället:
library(ggplot2)
library(reshape2)
library(randomcoloR)

#Functions

get_cancertype <- function(name){
  
  cancertype <- cohort[grep(name,cohort$sample), 2]
  
  return(cancertype)
}


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

#Modify data for plotting


#Optional Cluster data before plotting
#hc.sample = hclust(dist(cnvs), method = "complete")
#sample_order = rownames(cnvs)[hc.sample$order]


#Optional order data after cosine-mutsign clusterinng
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations(3s)")
clustered_samples <- read.table("Sample_In_Cluster_complete_7.txt",header=TRUE,stringsAsFactors = FALSE)
clustered_samples <- clustered_samples %>% arrange(Cluster) 
clustered_samples$Sample_id <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",clustered_samples$Sample_id)
clustered_samples <- clustered_samples %>% filter(Sample_id %in% rownames(cnvs))
sample_order <- clustered_samples$Sample_id





#Create dataframe in format for GGplot
longData <- melt(cnvs)
longData1 <- melt(snvs)
longData$cnvs <- as.factor(longData$value)
longData$snvs <- as.factor(longData1$value)
longData$cnvs <- as.factor(longData$value)
longData <- longData %>% select(Var1,Var2,cnvs,snvs)
colnames(longData) <- c("sample","gene","cnvs","snvs")


#Get number of mutations for each sample

get_mut <- function(sample_id){
  
  hit <- grep(sample_id,cohort$sample)
  
  return(cohort[hit,4])
}

#number_mut vector 

get_mut <- function(sample_id){
  
  hit <- grep(sample_id,cohort$sample)
  
  return(cohort[hit,4])
}


number_mut <- data.frame("mut"=sapply(sample_order,get_mut),"y_pos"=c(1:length(sample_order)),
                         "x_pos"=rep((ncol(cnvs)-4),nrow(cnvs)))





#This steps order y-axis in accordance to sample_order
longData$sample = factor(longData$sample, levels = sample_order)

#Try custom X-axis
y_names <- sample_order
y_name_colour <- data.frame(cluster = clustered_samples$Cluster)
y_colour <- data.frame(cluster=c(1:length(unique(y_name_colour$cluster))),
                      colour=distinctColorPalette(length(unique(y_name_colour$cluster))))

#or choose your own
y_colour$colour <- c("#000000","#800000","#000075",
                     "#3cb44b","#9A6324","#469990","#a9a9a9")

#y_colour$colour <- c("#000000","#800000","#000075",
                     #"#3cb44b","#9A6324","#469990","#a9a9a9")

y_colour$colour <- c("#000000","#800000","#000000","#800000","#000000","#800000","#000000")

y_name_colour <- left_join(y_name_colour,y_colour,by = "cluster")
y_name_colour <- as.vector(y_name_colour$colour)
    
                      
#get cancertype for these names
cancertype <- as.vector(sapply(y_names,get_cancertype))
y_names <- gsub("TCGA-","",gsub("-[A-Z0-9]*$","",y_names))
y_names <- paste(y_names,cancertype,sep="-")
#add number_of mutations
y_names <- paste(y_names,sep="-")


#Geom_hline position
hline_pos <- c()
for (cluster in unique(clustered_samples$Cluster)){
  
  pos <- tail(grep(cluster,clustered_samples$Cluster),n=1)
  pos <- pos + 0.5
  
  hline_pos <- c(hline_pos,pos)
  
}




#snv_shapes <- c(65:(65+length(snv_labels)))  #For alphabetic char
textcol <- "grey40"
#snv_shapes <- c(21:(21+length(snv_labels)))
snv_shapes <-  c(0,1,2,3,6)



#modified ggplot
ggplot(longData,aes(x=gene,y=sample,fill=cnvs))+
  geom_tile()+
  #redrawing tiles to remove cross lines from legend
  geom_tile(colour="white",size=0.25)+
  #remove axis labels, add title
  labs(x="",y="",title="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0),labels = y_names)+  
  #scale_y_discrete(expand=c(0,0))+
  #custom breaks on x-axis
  scale_x_discrete(expand=c(0,0))+
  #custom colours for cut levels and na values
  scale_fill_manual(values=c("#D55E00", "#E69F00","grey70","#56B4E9","#0072B2"))+
 
  #equal aspect ratio x and y axis
  #coord_fixed()+
 
  #Add snv data (remove NA(ie no mutation))
  geom_point(data=subset(longData,!snvs=="NA"), aes(shape=snvs),size=1)+
  scale_shape_manual(values=snv_shapes)+
  #Add vertical lines where clusters are
  geom_hline(yintercept=hline_pos)+
  
  #Add #mut text
  #geom_text(x=number_mut$x_pos,y=number_mut$y_pos,label=number_mut$number_mut)+
  #geom_text(data = number_mut, aes(label = mut,y = y_pos, x = x_pos))
  #set base size for all font elements
  #theme_grey(base_size=10)+
  #theme options
  theme(
    #remove legend title
    legend.title=element_blank(),
    #remove legend margin
    legend.spacing  = grid::unit(0,"cm"),
    #change legend text properties
    legend.text=element_text(colour=textcol,size=7,face="bold"),
    #change legend key height
    #legend.key.height=grid::unit(0.8,"cm"),
    #set a slim legend
    legend.key.width=grid::unit(0.4,"cm"),
    #set x axis text size and colour
    axis.text.x=element_text(colour="black",angle=90, hjust=1,vjust = 0.5,size=4),
    #set y axis text colour and adjust vertical justification
    axis.text.y=element_text(vjust = 0.2,colour=y_name_colour,size=5),
    #axis.text.y=element_text(vjust = 0.2,colour="black",size=5),
    #change axis ticks thickness
    axis.ticks=element_line(size=0.4),
    #change title font, size, colour and justification
    plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
    #remove plot background
    #plot.background=element_blank(),
    #remove plot border
    panel.border = element_rect(fill = NA))
    #panel.background= element_rect(colour = "black", fill=NA, size=1))
   




    library(grid)
    dev.off()
    
    p 
    
    
    
    
    grid.text("Distinct", x = unit(0.86, "npc"), y = unit(0.80, "npc"),
              gp = gpar(fontsize = 4))
 
    
    
    
    grid.text(number_mut$mut,x=as.numeric(number_mut$x_pos),
              y=as.numeric(number_mut$y_pos))

    


 