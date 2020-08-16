rm(list=ls())

library(tidyverse)
#require(ggforce)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gplots)
library(gtable)
library(grid)


directory <- "AT_Pathways_Files_TSV_Files"
ABSLFCTable12hpc <- read_tsv("Cluster_Analysis_ClustersbyABSMean.tsv")

# ABSLFCTable <- arrange(ABSLFCTable, desc(MeanL2FC))
# ABSLFCTable$MeanL2FC <- round(ABSLFCTable$MeanL2FC,3)

# ABSLFCTable12hpc <- filter(ABSLFCTable12hpc, cluster %in% c(2,3,4,6,9,10,11)) %>%
#    arrange(desc(MeanL2FC))

ABSLFCTable12hpc <- arrange(ABSLFCTable12hpc, desc(MeanL2FC))

ABSLFCTable12hpc$MeanL2FC <- round(ABSLFCTable12hpc$MeanL2FC,3) 

clusterOrder <- ABSLFCTable12hpc %>% 
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

AllPathways$Pathway <- str_extract(AllPathways$Pathway, "\\((.*?)\\)$") %>% str_sub(2, -2) %>% str_to_title()
tail(AllPathways$Pathway)


# Clean up text formatting
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "<I>", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "</I>", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Iv", "IV")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Vi", "VI")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Iii", "III")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Ii", "II")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Udp", "UDP")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Gdp", "UDP")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, " Of ", " of ")


# AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, " IV|III|II|I", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "^Superpathway of ", "")


#View(AllPathways)



AllPathways$cluster <- as.numeric(str_extract(AllPathways$FileName, "\\d+$"))


joinAllPathways <- AllPathways %>% 
  left_join(ABSLFCTable12hpc, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC))


pathwayHeatMapset <- joinAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value") %>%
  select(c("Pathway", clusterOrder + 1))

pathwayHeatMapsetMatrix <- as.matrix(pathwayHeatMapset[2:length(colnames(pathwayHeatMapset))])
rownames(pathwayHeatMapsetMatrix) <- pathwayHeatMapset$Pathway
head(pathwayHeatMapsetMatrix)

#pathwayHeatMapsetMatrix[which(pathwayHeatMapsetMatrix == 0)] <- 0.00000001

# This massive command filters pathways with only one p-value in all of the clusters
#filterPathways <- rownames(pathwayHeatMapsetMatrix[which(rowSums(as.data.frame(is.na(pathwayHeatMapsetMatrix))) < (length(colnames(pathwayHeatMapsetMatrix)) - 1)),])

##### Outputs

# PNG

## All Clusters in heatmap
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_", length(colnames(pathwayHeatMapsetMatrix)), "Clusters_heatmap_by_Cluster_12hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(pathwayHeatMapsetMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(pathwayHeatMapsetMatrix)), " Clusters\n(Clusters Ranked by MeanL2FC at 12hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()



# ABSLFCTable12hpc <- filter(ABSLFCTable12hpc, cluster %in% c(2,3,4,6,9,10,11)) %>%
#    arrange(desc(MeanL2FC))

## Upregulated clusters 
upregClusters <- ABSLFCTable12hpc %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()

upregMatrix <- pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters]

rownames(upregMatrix) %>%
str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_12hpc_All_UpregPathways.csv"))

# colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_heatmap_by_Cluster_12hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(upregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ",length(colnames(upregMatrix)), " Up-regulated Clusters\n(Clusters Ranked by MeanL2FC at 12hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()

# Upregulated clusters of interest ----------------------------------------

## Upregulated clusters of interest

clusterOfInterest <- c(2,4,6)

coiABSLFCTable12hpc <- ABSLFCTable12hpc %>% 
  filter(cluster %in% clusterOfInterest)

upregClusters <- coiABSLFCTable12hpc %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()


joinCOIAllPathways <- AllPathways %>% 
  filter(cluster %in% clusterOfInterest) %>%
  left_join(ABSLFCTable12hpc, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC)) %>% 
  unite(cluster, MeanL2FC, col = cluster, sep = ":") 


uppathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value")

uppathwayHeatMapsetMatrix <- as.matrix(uppathwayHeatMapset[2:length(colnames(uppathwayHeatMapset))])
rownames(uppathwayHeatMapsetMatrix) <- uppathwayHeatMapset$Pathway

upregMatrix <- uppathwayHeatMapsetMatrix[which(apply(uppathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters]

rownames(upregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_12hpc_COI_UpregPathways.csv"))

heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_of_Interest_heatmap_by_Cluster_12hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(upregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(upregMatrix)), " Clusters of Interest\n(Clusters Ranked by MeanL2FC at 12hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()


## Downregulated clusters 
downregClusters <- ABSLFCTable12hpc %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()


downregMatrix <- pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,downregClusters], 1, function(x) sum(is.na(x))) < length(downregClusters)),downregClusters]

rownames(downregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_12hpc_All_DownregPathways.csv"))


#colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_heatmap_by_Cluster_12hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(downregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(downregMatrix)), " Down-regulated Clusters\n(Clusters Ranked by MeanL2FC at 12hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()

# Downregulated clusters of interest ----------------------------------------

## Downregulated clusters of interest

clusterOfInterest <- c(9)

coiABSLFCTable12hpc <- ABSLFCTable12hpc %>% 
  filter(cluster %in% clusterOfInterest)

downregClusters <- coiABSLFCTable12hpc %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()


joinCOIAllPathways <- AllPathways %>% 
  filter(cluster %in% clusterOfInterest) %>%
  left_join(ABSLFCTable12hpc, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC)) %>% 
  unite(cluster, MeanL2FC, col = cluster, sep = ":") 


downpathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value")

downpathwayHeatMapsetMatrix <- as.matrix(downpathwayHeatMapset[2:length(colnames(downpathwayHeatMapset))])
rownames(downpathwayHeatMapsetMatrix) <- downpathwayHeatMapset$Pathway

downregMatrix <-downpathwayHeatMapsetMatrix

rownames(downregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_12hpc_COI_DownregPathways.csv"))

heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_of_Interest_heatmap_by_Cluster_12hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(downregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(downregMatrix)), " Down-regulated Clusters of Interest\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()




# Heatmap rank by mean LFC 24hpc ------------------------------------------


ABSLFCTable24hpc <- read_tsv("Cluster_Analysis_ClustersbyABSMean.tsv")

# ABSLFCTable24hpc <- filter(ABSLFCTable24hpc, cluster %in% c(2,3,4,6,9,10,11)) %>%
#    arrange(desc(MeanL2FC))

ABSLFCTable24hpc <- arrange(ABSLFCTable24hpc, desc(MeanL2FC))


ABSLFCTable24hpc$MeanL2FC <- round(ABSLFCTable24hpc$MeanL2FC,3) 

# ABSLFCTable <- arrange(ABSLFCTable, desc(MeanL2FC))
# ABSLFCTable$MeanL2FC <- round(ABSLFCTable$MeanL2FC,3)

# ABSLFCTable12hpc <- filter(ABSLFCTable12hpc, cluster %in% c(2,3,4,6,9,10,11)) %>%
#    arrange(desc(MeanL2FC))

ABSLFCTable24hpc <- arrange(ABSLFCTable24hpc, desc(MeanL2FC))

ABSLFCTable24hpc$MeanL2FC <- round(ABSLFCTable24hpc$MeanL2FC,3) 

clusterOrder <- ABSLFCTable24hpc %>% 
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

AllPathways$Pathway <- str_extract(AllPathways$Pathway, "\\((.*?)\\)$") %>% str_sub(2, -2) %>% str_to_title()
tail(AllPathways$Pathway)


# Clean up text formatting
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "<I>", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "</I>", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Iv", "IV")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Iii", "III")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Ii", "II")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Udp", "UDP")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "Gdp", "UDP")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, " Of ", " of ")


# AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, " IV|III|II|I", "")
AllPathways$Pathway <- str_replace_all(AllPathways$Pathway, "^Superpathway of ", "")


#View(AllPathways)



AllPathways$cluster <- as.numeric(str_extract(AllPathways$FileName, "\\d+$"))


joinAllPathways <- AllPathways %>% 
  left_join(ABSLFCTable24hpc, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC))


pathwayHeatMapset <- joinAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value") %>%
  select(c("Pathway", clusterOrder + 1))

pathwayHeatMapsetMatrix <- as.matrix(pathwayHeatMapset[2:length(colnames(pathwayHeatMapset))])
rownames(pathwayHeatMapsetMatrix) <- pathwayHeatMapset$Pathway
head(pathwayHeatMapsetMatrix)

#pathwayHeatMapsetMatrix[which(pathwayHeatMapsetMatrix == 0)] <- 0.00000001

# This massive command filters pathways with only one p-value in all of the clusters
#filterPathways <- rownames(pathwayHeatMapsetMatrix[which(rowSums(as.data.frame(is.na(pathwayHeatMapsetMatrix))) < (length(colnames(pathwayHeatMapsetMatrix)) - 1)),])

## All Clusters in heatmap
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_", length(colnames(pathwayHeatMapsetMatrix)), "Clusters_heatmap_by_Cluster_24_hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(pathwayHeatMapsetMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results", length(colnames(pathwayHeatMapsetMatrix)), " Clusters\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()


## Upregulated clusters 
upregClusters <- ABSLFCTable24hpc %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()

upregMatrix <- pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters]


rownames(upregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_24hpc_All_UpregPathways.csv"))


# colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_heatmap_by_Cluster_24hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(upregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(upregMatrix)), " Up-regulated Clusters\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()

# Upregulated clusters of interest ----------------------------------------

## Upregulated clusters of interest

clusterOfInterest <- c(4,9,10)

coiABSLFCTable24hpc <- ABSLFCTable24hpc %>% 
  filter(cluster %in% clusterOfInterest)

upregClusters <- coiABSLFCTable24hpc %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC > 0) %>%
  pull(OrderId) %>%
  as.numeric()


joinCOIAllPathways <- AllPathways %>% 
  filter(cluster %in% clusterOfInterest) %>%
  left_join(ABSLFCTable24hpc, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC)) %>% 
  unite(cluster, MeanL2FC, col = cluster, sep = ":") 


uppathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value")

uppathwayHeatMapsetMatrix <- as.matrix(uppathwayHeatMapset[2:length(colnames(uppathwayHeatMapset))])
rownames(uppathwayHeatMapsetMatrix) <- uppathwayHeatMapset$Pathway

upregMatrix <- uppathwayHeatMapsetMatrix[which(apply(uppathwayHeatMapsetMatrix[,upregClusters], 1, function(x) sum(is.na(x))) < length(upregClusters)),upregClusters]

rownames(upregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_24hpc_COI_UpregPathways.csv"))

heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_UpReg_", length(upregClusters), "Clusters_of_Interest_heatmap_by_Cluster_24hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(upregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(upregMatrix)), " Up-regulated Clusters of Interest\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()



## Downregulated clusters 
downregClusters <- ABSLFCTable24hpc %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()

downregMatrix <- pathwayHeatMapsetMatrix[which(apply(pathwayHeatMapsetMatrix[,downregClusters], 1, function(x) sum(is.na(x))) < length(downregClusters)),downregClusters]

rownames(downregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_24hpc_All_DownregPathways.csv"))

#colRange <- 6:length(colnames(pathwayHeatMapsetMatrix))
heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_heatmap_by_Cluster_24hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(downregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(downregMatrix)), " Down-regulated Clusters\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()

# Downregulated clusters of interest ----------------------------------------

## Downregulated clusters of interest

clusterOfInterest <- c(6, 3,11)

coiABSLFCTable24hpc <- ABSLFCTable24hpc %>% 
  filter(cluster %in% clusterOfInterest)

downregClusters <- coiABSLFCTable24hpc %>%
  arrange(desc(MeanL2FC)) %>%
  tibble::rownames_to_column("OrderId")%>% 
  arrange(desc(MeanL2FC)) %>%
  filter(MeanL2FC < 0) %>%
  pull(OrderId) %>%
  as.numeric()


joinCOIAllPathways <- AllPathways %>% 
  filter(cluster %in% clusterOfInterest) %>%
  left_join(ABSLFCTable24hpc, by = c("cluster" = "cluster")) %>% 
  arrange(desc(MeanL2FC)) %>% 
  unite(cluster, MeanL2FC, col = cluster, sep = ":") 


downpathwayHeatMapset <- joinCOIAllPathways[,c("Pathway", "permuted p-value", "cluster")] %>%
  spread("cluster", "permuted p-value")

downpathwayHeatMapsetMatrix <- as.matrix(downpathwayHeatMapset[2:length(colnames(downpathwayHeatMapset))])
rownames(downpathwayHeatMapsetMatrix) <- downpathwayHeatMapset$Pathway

downregMatrix <-downpathwayHeatMapsetMatrix

rownames(downregMatrix) %>%
  str_replace_all(" I+|V$", "") %>%
  unique() %>% 
  write.csv(file.path(directory, "AraCyc_24hpc_COI_DownregPathways.csv"))

heatmapOutfile <- file.path(directory, paste("AthalianaPathways_Sigresults_DownReg_", length(downregClusters), "Clusters_of_Interest_heatmap_by_Cluster_24hpc_rank_MeanL2FC.png", sep = ""))
png(heatmapOutfile, width = 1220, height = 1144, res = 140)
breaksList = seq(0, 15, by = 1)
color <- colorRampPalette(rev(brewer.pal(n = 5, name = "Greens")))(length(breaksList))
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.98, height=0.96, name="vp", just=c("right","top"))), action="prepend")
pheatmap(downregMatrix, 
         color = color, 
         na_col = 'white',
         fontsize = 6,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         legend_breaks = c(0, 0.01, 0.02, 0.03, 0.04, max(pathwayHeatMapsetMatrix, na.rm = TRUE)),
         fontsize_col = 8,
         legend_labels = c("0", "0.01", "0.02", "0.03", "0.04", "Permuted\n P-Value"),
         main = paste("AraCyc Pathway Results for ", length(colnames(downregMatrix)), " Down-regulated Clusters of Interest\n(Clusters Ranked by MeanL2FC at 24hpc)", sep = ""))
setHook("grid.newpage", NULL, "replace")
grid.text("Clusters", x=0.35, y=-0.01, gp=gpar(fontsize=12))
dev.off()



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

