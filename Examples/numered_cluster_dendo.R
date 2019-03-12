# try out plot

#library



#Main
rm(list=ls())
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/Global_mutsign/High_mutations(3s)")
#scource("C:/Users/Nils_/MutationalPatterns/R")
#source("C:/Users/Nils_/MutationalPatterns/R")
load("high_mut_matrix.rda")
colnames(mutational_matrix) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(mutational_matrix))
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3")
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

#function
plot_cosine_heatmap1 = function(cos_sim_matrix, col_order, cluster_rows = TRUE, method = "complete", plot_values = FALSE)
{
  # check explained argument
  if(class(cos_sim_matrix) != "matrix")
  {stop("cos_sim_matrix must be a matrix")}
  # matrix should have row and colnames
  if(length(colnames(cos_sim_matrix)) == 0)
  {stop("cos_sim_matrix is missing colnames")}
  if(length(rownames(cos_sim_matrix)) == 0)
  {stop("cos_sim_matrix is missing rownames")}
  # if no signature order is provided, use the order as in the input matrix
  if(missing(col_order))
  {
    col_order = colnames(cos_sim_matrix)
  }
  # check col_order argument
  if(class(col_order) != "character")
  {stop("col_order must be a character vector")}
  if(length(col_order) != ncol(cos_sim_matrix))
  {stop("col_order must have the same length as the number of signatures in the explained matrix")}
  
  # if cluster samples is TRUE, perform clustering
  if(cluster_rows == TRUE)
  {
    # cluster samples based on eucledian distance between relative contribution
    hc.sample = hclust(dist(cos_sim_matrix), method = method)
    # order samples according to clustering
    sample_order = rownames(cos_sim_matrix)[hc.sample$order]
  }
  else
  {
    sample_order = rownames(cos_sim_matrix)
  }
  
  Cosine.sim = NULL
  Signature = NULL
  Sample = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL
  print("it works")
  # melt
  cos_sim_matrix.m = melt(cos_sim_matrix)
  # assign variable names
  colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  
  # change factor levels to the correct order for plotting
  cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature, levels = col_order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, levels = sample_order)
  # plot heatmap
  heatmap = ggplot(cos_sim_matrix.m, aes(x=Signature, y=Sample, fill=Cosine.sim, order=Sample)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1)) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=5)) +
    labs(x=NULL, y=NULL)
  # if plot_values is TRUE, add values to heatmap
  if (plot_values)
  {
    heatmap = heatmap + geom_text(aes(label = round(Cosine.sim, 2)), size = 3)
  }
  
  # if cluster samples is TRUE, make dendrogram
  if(cluster_rows == TRUE)
  {
    # get dendrogram
    dhc = as.dendrogram(hc.sample)
    # rectangular lines
    ddata = dendro_data(dhc, type = "rectangle")
    #Add colour 
    
    clust    <- cutree(hc.sample,k=7)                   
    clust.df <- data.frame(label=names(clust), cluster=factor(clust))
    ddata[["labels"]] <- merge(ddata[["labels"]],clust.df, by="label")
    
    # plot dendrogram of hierachical clustering
    dendrogram = ggplot(segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster),size=3)+
      coord_flip() + 
      scale_y_reverse(expand = c(0.2, 0)) + 
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.line.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank(),
            legend.position = "none")
# combine plots
plot_final = plot_grid(dendrogram, heatmap, align='h', rel_widths=c(0.3,1),rel_heights = c(2,1))
  }
  else
  {
    plot_final = heatmap +
      # reverse order of the samples such that first is up
      ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  
  return(plot_final)
}


print(plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE))


#Library 
library(reshape2)
library(ggplot2)
library(ggdendro)
library(cowplot)




hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]

print(plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE))








#plotting
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
plot(sample_cluster)
N<-7
method <- "complete"
cluster_groups <- cutree(sample_cluster,k=N)

#Plot with found clusters
plot(sample_cluster, cex = 0.6)
rect.hclust(sample_cluster, k = N, border = 2:5)

#Test 

hc = hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
dendr <- dendro_data(hc, type="rectangle")
clust    <- cutree(hc,k=7)                    # find 2 clusters
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")









ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
            size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())












