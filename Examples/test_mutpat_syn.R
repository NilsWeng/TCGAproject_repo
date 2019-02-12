
library(BSgenome)
library(GenomicRanges)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(MutationalPatterns)
#Defining reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

# cleaning the workspace
rm(list = ls(all = TRUE))

setwd("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274")
subtypes <- c("KIRC", "KIRP", "KICH")
#Read maf file per cancer type and split to per tumor vcf file:
maf <- c(("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274/KIRC.maf"),("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274/KIRP.maf"),("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274/KICH.maf"))

vcffile=vector()
subtype=vector()
sample=vector()

for (s in 1:length(subtypes)){
  data=read.table(maf[s], header=TRUE )
  tumors=levels(factor(data$Tumor_Sample_Barcode))
  sample=c(sample, tumors)
  for (t in tumors){
    tumormaf=data[which(data$Tumor_Sample_Barcode %in% t),]
    vcfdata <- matrix(".", nrow=nrow(tumormaf), ncol = 10)
    columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", t)
    colnames(vcfdata) <- columns
    vcfdata[,1]=as.character(tumormaf$Chromosome)
    vcfdata[,2]=as.character(tumormaf$Start_Position)
    vcfdata[,4]=as.character(tumormaf$Reference_Allele)
    vcfdata[,5]=as.character(tumormaf$Tumor_Seq_Allele2)
    vcfdata[,9]="GT"
    vcfdata[,10]="1/0"
    outfile=paste(t,".vcf", sep="")
    outfile=paste(subtypes[s], outfile, sep="/")
    write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
    vcffile <- c(vcffile, outfile)
  }
  subtype <- c(subtype, rep(subtypes[s], length(tumors)))
}

cohort=data.frame("vcf"=vcffile, "subtype"=subtype, "sample"=sample)

auto = extractSeqlevelsByGroup(species="Homo_sapiens", 
                               style="UCSC",
                               group="auto")

#Code for cohort-level VCF files
kircvcf <- ("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274/KIRC.vcf")
kirpvcf <-  ("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274/KIRP.vcf")
kichvcf <- ("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/syn7824274/KICH.vcf")
vcfer=c(kircvcf, kirpvcf , kichvcf )
sam=c("KIRC", "KIRP","KICH")
vcfs <- read_vcfs_as_granges(vcfer, sam, ref_genome)
tri_n_matrix_sam = mut_matrix(vcfs,ref_genome)
extract_3 = extract_signatures(tri_n_matrix_sam, rank = 3)
extract_4 = extract_signatures(tri_n_matrix_sam, rank = 4)
pdf("extracted_3.pdf")
plot_96_profile(extract_3$signatures)
dev.off()
pdf("extracted_contribution.pdf")
plot_contribution(extract_3$contribution, extract_3$signature, mode = "absolute", coord_flip = T)
dev.off()

#Read individual vcfs

vcfs_kirc <- cohort$vcf[which(cohort$subtype %in% "KIRC")]
vcfs_kirp <- cohort$vcf[which(cohort$subtype %in% "KIRP")]
vcfs_kich <- cohort$vcf[which(cohort$subtype %in% "KICH")]

#Number of mutations per sample:
library(R.utils)
nmut_kirc=sapply(vcfs_kirc,countLines)-1
nmut_kirp=sapply(vcfs_kirp,countLines)-1
nmut_kich=sapply(vcfs_kich,countLines)-1
nmut=c(nmut_kirc, nmut_kirp, nmut_kich)
cohort$nmut=nmut

View(cohort[order(cohort$nmut, decreasing = TRUE),])

#Here is code for per-sample-vcf:
vkirc <- read_vcfs_as_granges(as.character(vcfs_kirc), cohort$sample[which(cohort$subtype %in% "KIRC")], ref_genome)
vkirp <- read_vcfs_as_granges(as.character(vcfs_kirp), cohort$sample[which(cohort$subtype %in% "KIRP")], ref_genome)
vkich <- read_vcfs_as_granges(as.character(vcfs_kich), cohort$sample[which(cohort$subtype %in% "KICH")], ref_genome)


vkirc_auto = lapply(vkirc, function(x) keepSeqlevels(x, auto, pruning.mode='coarse'))
vkirp_auto = lapply(vkirp, function(x) keepSeqlevels(x, auto, pruning.mode='coarse'))
vkich_auto = lapply(vkich, function(x) keepSeqlevels(x, auto, pruning.mode='coarse'))


tri_n_kirc = mut_matrix(vkirc_auto,ref_genome)
tri_n_kirp = mut_matrix(vkirp_auto,ref_genome)
tri_n_kich = mut_matrix(vkich_auto,ref_genome)

fit_kirc = fit_to_signatures(tri_n_kirc, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")])
fit_kirp = fit_to_signatures(tri_n_kirp, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")])
fit_kich = fit_to_signatures(tri_n_kich, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")])

fit_kirc_all = fit_to_signatures(tri_n_kirc, cancer_signatures)
fit_kirp_all = fit_to_signatures(tri_n_kirp, cancer_signatures)
fit_kich_all = fit_to_signatures(tri_n_kich, cancer_signatures)


pdf("fit_cosmic_kirc.pdf")
plot_contribution(fit_kirc$contribution[,1:100], cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
plot_contribution(fit_kirc$contribution[,101:200], cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
plot_contribution(fit_kirc$contribution[,201:ncol(fit_kirc$contribution)], cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
dev.off()

pdf("fit_cosmic_kirp.pdf")
plot_contribution(fit_kirp$contribution[,1:100], cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
plot_contribution(fit_kirp$contribution[,101:200], cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
plot_contribution(fit_kirp$contribution[,201:ncol(fit_kirp$contribution)], cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
dev.off()

pdf("fit_cosmic_kich.pdf")
plot_contribution(fit_kich$contribution, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
dev.off()

#How much of the original profiles are reconstructed with cosmic 1+5+4? 
#and how much are reconstructed with all cosmic signatures?
cos_sim_kirc = diag(cos_sim_matrix(fit_kirc$reconstructed, tri_n_kirc))
cos_sim_kirc_all = diag(cos_sim_matrix(fit_kirc_all$reconstructed, tri_n_kirc))
cos_sim_kirp = diag(cos_sim_matrix(fit_kirp$reconstructed, tri_n_kirp))
cos_sim_kirp_all = diag(cos_sim_matrix(fit_kirp_all$reconstructed, tri_n_kirp))
cos_sim_kich = diag(cos_sim_matrix(fit_kich$reconstructed, tri_n_kich))
cos_sim_kich_all = diag(cos_sim_matrix(fit_kich_all$reconstructed, tri_n_kich))
cos_sim_fits=cbind(cos_sim_kirc, cos_sim_kirc_all, cos_sim_kirp, cos_sim_kirp_all, cos_sim_kich, cos_sim_kich_all)
names=c("KIRC 1,5,4", "KIRC all", "KIRP 1,5,4","KIRP all", "KICH 1,5,4", "KICH all")
pdf("cos_sim_rek_obs.pdf")
par(cex.axis=0.8)
boxplot(cos_sim_kirc, cos_sim_kirc_all, cos_sim_kirp, cos_sim_kirp_all, cos_sim_kich, cos_sim_kich_all, names=names, main="cosine similarity observed vs reconstructed tri_n_spectrum", las=2)
dev.off()





tri_n_syn = tri_n_matrix + 0.0001
estimate_syn = nmf(tri_n_syn, rank=2:5, method="brunet", nrun=100, seed=123456)
pdf("rank_estimate.pdf")
plot(estimate_syn)
dev.off()


sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
cancer_signatures = read.table(sp_url, sep = "\t", header = T)
# reorder (to make the order of the trinucleotide changes the same)
cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
# only signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])
sig_hclust_cs = cluster_signatures(cancer_signatures)
col_order_cs = colnames(cancer_signatures)[sig_hclust_cs$order]



similarity_cosmic=cos_sim_matrix(extract_3$signature,cancer_signatures)
rownames(similarity_cosmic)=samples
pdf("similarity_to_cosmic_signatures.pdf")
plot_cosine_heatmap(similarity_cosmic, col_order_cs, cluster_rows = FALSE)
dev.off()



#Fit mutation matrix to cancer signatures.
fit_res = fit_to_signatures(tri_n_matrix, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")])

pdf("contribution_cosmic_KIRC.pdf")
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = T, mode = "absolute")
plot_contribution(fit_res$contribution[select,], cancer_signatures[,select], coord_flip = T, mode = "relative")
dev.off()

signatures=data.frame("sample" = colnames(fit_res$contribution), "Signature.1"=fit_res$contribution[1,], "Signature.5"=fit_res$contribution[2,], "Signature.4"=fit_res$contribution[3,])
colnames(signatures)[1] <- "sample"
contribution <- merge(signatures, cohort, by = "sample")
nmut=contribution$Signature.1+contribution$Signature.5+contribution$Signature.4
contribution$nmut <- contribution$Signature.1+contribution$Signature.5+contribution$Signature.4
sorted <- contribution[order(contribution$nmut, decreasing=TRUE),]

contribution_subtype=fit_res$contribution
colnames(contribution_subtype)=subtype
hypermut_subtypename=contribution_subtype[,which(contribution$nmut>100)]
hypermut_samplename=fit_res$contribution[,which(contribution$nmut>100)]
pdf("contribuition_cosmic_hypermut.pdf")
plot_contribution(hypermut_subtypename, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
plot_contribution(hypermut_samplename, cancer_signatures[,c("Signature.1", "Signature.5","Signature.4")], coord_flip = T, mode = "absolute")
dev.off()



#Read our clonal expansion data:
matchlist=read.table("/Users/malla74/projects/M_eriksson_1506/matchList.txt", header=TRUE)
sorted_by_tissue=matchlist[with(matchlist, order(tissue, age_years, individual)),]
individual = sorted_by_tissue$individual
age = sorted_by_tissue$age_class
samples = sorted_by_tissue$clone
tissue = sorted_by_tissue$tissue
vcf_files = paste("/Users/malla74/projects/M_eriksson_1506/mutationalPatterns/vcf_filtered_clean_20180226/ind_", sorted_by_tissue$individual, sep='')
vcf_files = paste(vcf_files, "_c_", sep='')
vcf_files = paste(vcf_files, sorted_by_tissue$clone, sep='')
vcf_files = paste(vcf_files, "_snps_filtered_clean.vcf", sep='')
vcfs_clone = read_vcfs_as_granges(vcf_files, samples, genome = ref_genome)

vcfs_clone = lapply(vcfs_clone, function(x) keepSeqlevels(x, auto, pruning.mode='coarse'))

tri_n_clone = mut_matrix(vcfs_clone,ref_genome)

fit_clone_extracted = fit_to_signatures(tri_n_clone, extract_3$signatures)
pdf("contribution_clones_extracted.pdf")
plot_contribution(fit_clone_extracted$contribution, extract_3$signatures, coord_flip = T, mode = "absolute")
dev.off()




#Cluster on tri-n-profiles
tri_n_all=cbind(tri_n_kirc, tri_n_kirp, tri_n_kich)
abs_mut=colSums(tri_n_all)
rel_tri_n_all=t(t(tri_n_all)/abs_mut)

#Color categories
library(RColorBrewer)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
collengts <- c(length(which(subtype %in% subtypes[1])), length(which(subtype %in% subtypes[2])), length(which(subtype %in% subtypes[3])))
subcols = col_vector[c(17,18,20)]

colSubtype[1:collengts[1]]=subcols[1]
colSubtype[(collengts[1]+1):(collengts[1]+collengts[2])]=subcols[2]
colSubtype[(collengts[1]+collengts[2]+1):(collengts[1]+collengts[2]+collengts[3])]=subcols[3]



color_legend=c("Cancer type:","KIRC","KIRP","KICH")
color_code=c("white",subcols)
#colLab=cbind(colInd, colTissue, colAge)[all,]
#colnames(colLab)=c("Individual","Tissue","Age")
library(heatmap3)
pdf("cluster_by_tri_n_matrix.pdf")
hm=heatmap3(rel_tri_n_all,  ColSideColors=colSubtype,Rowv=NA,cexRow=0.4, cexCol=0.4)
legend("left",legend=color_legend, col=color_code, pch=15, )
dev.off()
