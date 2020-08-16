rm(list=ls())

library(tidyverse)
#require(ggforce)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gplots)
library(gtable)
library(grid)


directory <- "Vitispathways_Analysis/TSV_Files"
ABSLFCTable <- read_tsv("Cluster_Analysis_ClustersbyABSMean.tsv")

ABSLFCTable$MeanL2FC <- round(ABSLFCTable$MeanL2FC,3)

clusterOrder <- ABSLFCTable %>% 
  arrange(desc(MeanL2FC)) %>% 
  select(cluster) %>%
  pull()

#clusterOrder <- clusterOrder[!clusterOrder %in% c("1:0.147", "10:0.069")]

pathwaysFiles <- list.files(directory, full.names = TRUE, pattern = ".tsv$")

pathFileList <- c()

for (file in pathwaysFiles){
  
  fileRoot <- file %>% 
    str_replace(".*/", "") %>%
    str_replace("\\.tsv", "")
  
  pathFileList <- c(fileRoot, pathFileList)
  
  assign(fileRoot, read_tsv(file))
  assign(fileRoot, cbind("FileName"=rep(fileRoot, nrow(get(fileRoot))), get(fileRoot)))

  
}

AllPathways <- bind_rows(mget(pathFileList)) %>% 
  filter(`permuted p-value` < 0.05)

AllPathways$Pathway <- str_replace(AllPathways$Pathway, "vv\\d{1,5}","")
head(AllPathways)
AllPathways$cluster <- as.numeric(str_extract(AllPathways$FileName, "\\d+$"))


joinAllPathways <- AllPathways %>% 
  left_join(ABSLFCTable, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC))


pathwayHeatMapset <- joinAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value") %>%
  select("Pathway", clusterOrder + 1)

pathwayHeatMapsetMatrix <- as.matrix(pathwayHeatMapset[2:length(colnames(pathwayHeatMapset))])
rownames(pathwayHeatMapsetMatrix) <- pathwayHeatMapset$Pathway
head(pathwayHeatMapsetMatrix)


# # Upregulated clusters of interest ----------------------------------------
# 
# ## Upregulated clusters of interest
# 
# clusterOfInterest <- c(3)
# 
# coiABSLFCTable <- ABSLFCTable %>% 
#   filter(cluster %in% clusterOfInterest)
# 
# upregClusters <- coiABSLFCTable %>%
#   arrange(desc(MeanL2FC)) %>%
#   tibble::rownames_to_column("OrderId")%>% 
#   arrange(desc(MeanL2FC)) %>%
#   filter(MeanL2FC > 0) %>%
#   pull(OrderId) %>%
#   as.numeric()
# 
# 
# joinCOIAllPathways <- AllPathways %>% 
#   filter(cluster %in% clusterOfInterest) %>%
#   left_join(ABSLFCTable, by = c("cluster" = "cluster")) %>% 
#   arrange(desc(MeanL2FC)) %>% 
#   unite(cluster, MeanL2FC, col = cluster, sep = ":") 
# 
# 
# uppathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
#   spread("cluster", "permuted p-value")
# 
# uppathwayHeatMapsetMatrix <- as.matrix(uppathwayHeatMapset[2:length(colnames(uppathwayHeatMapset))])
# rownames(uppathwayHeatMapsetMatrix) <- uppathwayHeatMapset$Pathway
# 
# upregMatrix <- uppathwayHeatMapsetMatrix
# 
# rownames(upregMatrix) %>%
#   write.csv(file.path(directory, "VitisNet_COI_UpregPathways.csv"))
# 
# heatmapOutfile <- file.path(directory, paste("VitisPathways_Sigresults_UpReg_", length(clusterOfInterest), "Clusters_of_Interest_heatmap_by_Cluster.pdf", sep = ""))
# breaksList = seq(0, 15, by = 1)
# 
# color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
# pheatmap(upregMatrix, 
#          color = color, 
#          na_col = 'white',
#          fontsize = 5,
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE,
#          main = paste(length(clusterOfInterest), " UpReg Clusters of Interest Vitis Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
#          legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
#          filename = heatmapOutfile)
# 
# # Downregulated clusters of interest ----------------------------------------
# 
# ## Downregulated clusters of interest
# 
# clusterOfInterest <- c(2)
# 
# coiABSLFCTable <- ABSLFCTable %>% 
#   filter(cluster %in% clusterOfInterest)
# 
# downregClusters <- coiABSLFCTable %>%
#   arrange(desc(MeanL2FC)) %>%
#   tibble::rownames_to_column("OrderId")%>% 
#   arrange(desc(MeanL2FC)) %>%
#   filter(MeanL2FC < 0) %>%
#   pull(OrderId) %>%
#   as.numeric()
# 
# 
# joinCOIAllPathways <- AllPathways %>% 
#   filter(cluster %in% clusterOfInterest) %>%
#   left_join(ABSLFCTable, by = c("cluster" = "cluster")) %>% 
#   arrange(desc(MeanL2FC)) %>% 
#   unite(cluster, MeanL2FC, col = cluster, sep = ":") 
# 
# 
# downpathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
#   spread("cluster", "permuted p-value")
# 
# downpathwayHeatMapsetMatrix <- as.matrix(downpathwayHeatMapset[2:length(colnames(downpathwayHeatMapset))])
# rownames(downpathwayHeatMapsetMatrix) <- downpathwayHeatMapset$Pathway
# 
# downregMatrix <- downpathwayHeatMapsetMatrix
# 
# rownames(downregMatrix) %>%
#   write.csv(file.path(directory, "VitisNet_COI_DownregPathways.csv"))
# 
# heatmapOutfile <- file.path(directory, paste("VitisPathways_Sigresults_DownReg_", length(clusterOfInterest), "Clusters_of_Interest_heatmap_by_Cluster.pdf", sep = ""))
# breaksList = seq(0, 15, by = 1)
# 
# color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
# pheatmap(downregMatrix, 
#          color = color, 
#          na_col = 'white',
#          fontsize = 5,
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE,
#          main = paste(length(clusterOfInterest), " DownReg Clusters of Interest Vitis Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
#          legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
#          filename = heatmapOutfile)
# 
# # Maxregulated clusters of interest ----------------------------------------
# 
# ## Maxregulated clusters of interest
# 
# clusterOfInterest <- c(13)
# 
# coiABSLFCTable <- ABSLFCTable %>% 
#   filter(cluster %in% clusterOfInterest)
# 
# MaxregClusters <- coiABSLFCTable %>%
#   arrange(desc(MeanL2FC)) %>%
#   tibble::rownames_to_column("OrderId")%>% 
#   arrange(desc(MeanL2FC)) %>%
#   filter(MeanL2FC < 0) %>%
#   pull(OrderId) %>%
#   as.numeric()
# 
# 
# joinCOIAllPathways <- AllPathways %>% 
#   filter(cluster %in% clusterOfInterest) %>%
#   left_join(ABSLFCTable, by = c("cluster" = "cluster")) %>% 
#   arrange(desc(MeanL2FC)) %>% 
#   unite(cluster, MeanL2FC, col = cluster, sep = ":") 
# 
# 
# MaxpathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
#   spread("cluster", "permuted p-value")
# 
# MaxpathwayHeatMapsetMatrix <- as.matrix(MaxpathwayHeatMapset[2:length(colnames(MaxpathwayHeatMapset))])
# rownames(MaxpathwayHeatMapsetMatrix) <- MaxpathwayHeatMapset$Pathway
# 
# MaxregMatrix <- MaxpathwayHeatMapsetMatrix
# 
# rownames(MaxregMatrix) %>%
#   write.csv(file.path(directory, "VitisNet_COI_MaxregPathways.csv"))
# 
# heatmapOutfile <- file.path(directory, paste("VitisPathways_Sigresults_MaxReg_", length(clusterOfInterest), "Clusters_of_Interest_heatmap_by_Cluster.pdf", sep = ""))
# breaksList = seq(0, 15, by = 1)
# 
# color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
# pheatmap(MaxregMatrix, 
#          color = color, 
#          na_col = 'white',
#          fontsize = 5,
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE,
#          main = paste(length(clusterOfInterest), " MaxReg Clusters of Interest Vitis Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
#          legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
#          filename = heatmapOutfile)

#pathwayHeatMapsetMatrix[which(pathwayHeatMapsetMatrix == 0)] <- 0.00000001

# This massive command filters pathways with only one p-value in all of the clusters
#filterPathways <- rownames(pathwayHeatMapsetMatrix[which(rowSums(as.data.frame(is.na(pathwayHeatMapsetMatrix))) < (length(colnames(pathwayHeatMapsetMatrix)) - 1)),])

##### Outputs


# Output TIFF -------------------------------------------------------------

## All Clusters in heatmap
heatmapOutfile <- file.path(directory, paste("VitisPathways_Sigresults_", length(colnames(pathwayHeatMapsetMatrix)), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.tiff", sep = ""))
tiff(heatmapOutfile, width = 1220, height = 1144, res = 140, compression = "lzw")
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(pathwayHeatMapsetMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 8,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_col = 8,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste(length(colnames(pathwayHeatMapsetMatrix)), " Clusters VitisPathways Results\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()

## Upregulated clusters 
upregClusters <- ABSLFCTable %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()

upregMatrix <- pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters]
write.csv(rownames(upregMatrix), file.path(directory, "VitisNet_24hpc_UpregPathways.csv"))


# colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("VitisPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.tiff", sep = ""))

tiff(heatmapOutfile, width = 1220, height = 1144, res = 140, compression = "lzw")
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(upregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 8,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_col = 8,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(upregMatrix, na.rm = TRUE)),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("VitisPathways Results - ", length(upregClusters), " Upregulated Clusters\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=10))
dev.off()


## Downregulated clusters 
downregClusters <- ABSLFCTable %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()

downregMatrix <- pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,downregClusters], 1, function(x) sum(is.na(x))) < length(downregClusters)),downregClusters]
write.csv(rownames(downregMatrix), file.path(directory, "VitisNet_24hpc_DownregPathways.csv"))

#colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("VitisPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.tiff", sep = ""))

tiff(heatmapOutfile, width = 1220, height = 1144, res = 140, compression = "lzw")
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(downregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 8,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_col = 8,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(downregMatrix, na.rm = TRUE)),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("VitisPathways Results - ", length(downregClusters), " Downregulated Clusters\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=10))
dev.off()



# This code is used to edit the x-axis labels -----------------------------

# trace(pheatmap:::draw_colnames, edit=TRUE)
# 
# function (coln, gaps, ...)
# {
#   coord = find_coordinates(length(coln), gaps)
#   x = coord$coord - 0.5 * coord$size
#   res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,
#                                                         "bigpts"), vjust = 0.5, hjust = 0, rot = 300, gp = gpar(...))
#   return(res)
# }

