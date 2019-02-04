readCNVs <-
function(filename){
    
    suppressPackageStartupMessages(library("GenomicRanges"))
    cnv_DF <- read.table(filename,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    
    #convert to GRanges object
    cnv_GRanges <- GRanges(seqnames =  cnv_DF$Chromosome,
                           ranges=IRanges(cnv_DF$Start,cnv_DF$End),
                           strand='*',
                           
                           Num_Probes=cnv_DF$Num_Probes,
                           Segment_Mean=cnv_DF$Segment_Mean,
                           GDC_Aliquot=cnv_DF$GDC_Aliquot)
    
    return(cnv_GRanges)
    
}
