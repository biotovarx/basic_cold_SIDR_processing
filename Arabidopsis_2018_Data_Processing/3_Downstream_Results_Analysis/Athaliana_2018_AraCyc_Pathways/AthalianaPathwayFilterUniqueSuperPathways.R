

# Initiate Script ---------------------------------------------------------

rm(list=ls())

library(tidyverse)
#require(ggforce)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gplots)


directory <- "AT_Pathways_Files_TSV_Files"
ABSLFCTable <- read_tsv("Cluster_Analysis_ClustersbyABSMean.tsv")

ABSLFCTable <- arrange(ABSLFCTable, desc(MeanL2FC))
ABSLFCTable$MeanL2FC <- round(ABSLFCTable$MeanL2FC,3)


clusterOrder <- ABSLFCTable %>% 
  arrange(desc(MeanL2FC)) %>% 
  unite(cluster, MeanL2FC, col = cluster, sep = ":") %>%
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

AllPathways$Pathway <- str_extract(AllPathways$Pathway, "\\((.*?)\\)$") %>% str_sub(2, -2) %>% str_to_title()
tail(AllPathways$Pathway)


# Clean Up Text Formatting ------------------------------------------------
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "<I>", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "</I>", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Iv", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Iii", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Ii", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Udp", "UDP")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Gdp", "UDP")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, " Of ", " of ")



# AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, " IV|III|II|I", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "^Superpathway of", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "I |I$", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Vi|V", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "\\(.*\\)", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "\\s+$", "")





AllPathways <- AllPathways %>% 
  left_join(ABSLFCTable, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC))


# Filter unique pathways
AllPathways %>% 
  filter(MeanL2FC > 0) %>%
  select(Pathway) %>%
  unique()

AllPathways %>% 
  filter(MeanL2FC < 0) %>%
  select(Pathway) %>%
  unique()



View(unique(AllPathways$Pathway))



AllPathways$cluster <- as.numeric(str_extract(AllPathways$FileName, "\\d+$"))









# Pathway Clustering ------------------------------------------------------





# pathwayHeatMapset <- AllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
#   spread("cluster", "permuted p-value") %>%
#   select(c("Pathway", clusterOrder))

# pathwayHeatMapsetMatrix <- as.matrix(pathwayHeatMapset[2:length(colnames(pathwayHeatMapset))])
# rownames(pathwayHeatMapsetMatrix) <- pathwayHeatMapset$Pathway
# head(pathwayHeatMapsetMatrix)

#pathwayHeatMapsetMatrix[which(pathwayHeatMapsetMatrix == 0)] <- 0.00000001

# This massive command filters pathways with only one p-value in all of the clusters
#filterPathways <- rownames(pathwayHeatMapsetMatrix[which(rowSums(as.data.frame(is.na(pathwayHeatMapsetMatrix))) < (length(colnames(pathwayHeatMapsetMatrix)) - 1)),])

##### Outputs

# PDF

# ## All Clusters in heatmap
# heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_", length(colnames(pathwayHeatMapsetMatrix)), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.pdf", sep = ""))
# breaksList = seq(0, 15, by = 1)
# color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
# pheatmap(pathwayHeatMapsetMatrix, 
#          color = color, 
#          na_col = 'white',
#          fontsize = 4,
#          legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix)),
#          cluster_rows = FALSE, 
#          cluster_cols = FALSE,
#          main = paste(length(colnames(pathwayHeatMapsetMatrix)), " Clusters Athaliana Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
#          legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value\n"),
#          filename = heatmapOutfile)


## Upregulated clusters 
upregClusters <- ABSLFCTable %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()







# colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.pdf", sep = ""))
breaksList = seq(0, 15, by = 1)

color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
pheatmap(pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters], 
         color = color, 
         na_col = 'white',
         fontsize = 5,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = paste(length(upregClusters), " UpReg Clusters Athaliana Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
         filename = heatmapOutfile)

## Downregulated clusters 
downregClusters <- ABSLFCTable %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()

#colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.pdf", sep = ""))
breaksList = seq(0, 15, by = 1)

color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
pheatmap(pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,downregClusters], 1, function(x) sum(is.na(x))) < length(downregClusters)),downregClusters], 
         color = color, 
         na_col = 'white',
         fontsize = 5,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = paste(length(downregClusters), " DownReg Clusters Athaliana Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
         filename = heatmapOutfile)


# PNG

## All Clusters in heatmap
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_", length(colnames(pathwayHeatMapsetMatrix)), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.png", sep = ""))
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
pheatmap(pathwayHeatMapsetMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 4,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = paste(length(colnames(pathwayHeatMapsetMatrix)), " Clusters Athaliana Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
         filename = heatmapOutfile)


## Upregulated clusters 
upregClusters <- ABSLFCTable %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()

# colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.png", sep = ""))
breaksList = seq(0, 15, by = 1)

color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
pheatmap(pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters], 
         color = color, 
         na_col = 'white',
         fontsize = 5,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = paste(length(upregClusters), " UpReg Clusters Athaliana Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
         filename = heatmapOutfile)

## Downregulated clusters 
downregClusters <- ABSLFCTable %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()

#colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_heatmap_by_Cluster_rank_MeanL2FC.png", sep = ""))
breaksList = seq(0, 15, by = 1)

color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
pheatmap(pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,downregClusters], 1, function(x) sum(is.na(x))) < length(downregClusters)),downregClusters], 
         color = color, 
         na_col = 'white',
         fontsize = 5,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = paste(length(downregClusters), " DownReg Clusters Athaliana Pathways Results\nClusters Ranked by MeanL2FC", sep = ""),
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\np-value"),
         filename = heatmapOutfile)




## This code is used to edit the x-axis labels

# trace(pheatmap:::draw_colnames, edit=TRUE)
# 
# 
# function (coln, gaps, ...)
# {
#   coord = find_coordinates(length(coln), gaps)
#   x = coord$coord - 0.5 * coord$size
#   res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,
#                                                         "bigpts"), vjust = 0.5, hjust = 0, rot = 300, gp = gpar(...))
#   return(res)
# }

# trace(pheatmap:::draw_rownames, edit=TRUE)
# 
# 
# function (rown, gaps, ...) 
# {
#   coord = find_coordinates(length(rown), gaps)
#   y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
#   res = textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, 
#                                           hjust = 0, rot = 45, gp = gpar(...))
#   return(res)
# }

