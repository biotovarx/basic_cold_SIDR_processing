
# Read in results from Arabidopsis AgriGO analysis clusters, filter,and create some nice plots


# Clear environment
rm(list=ls())

# CRAN packages
library(tidyverse)
library(ggforce)
library(gtable)  
library(grid)
library(scales)

# Directory to which output will be saved and check for and if needed create output directory:
output_dir <- paste(getwd(), "/Athalaian_AgriGO_Processing", sep="")



if(file_test(op='-d',output_dir)){
  print(noquote(paste('Output will be saved to the existing directory ',output_dir,'.',sep='')))
  print(noquote('___________________________________________________________________________'))
  print(noquote(''))
}else{
  print(noquote(paste('Creating directory',output_dir,'for output.')))
  print(noquote('___________________________________________________________________________'))
  print(noquote(''))
  dir.create(output_dir)
}


# Read and filter GO files ------------------------------------------------

# Find cluster files, read, combine single dataframe
clusterDirectory <- "Clustered_Genes_All_experiments_13_clusters"
clusterFileEnding <- "_genes.csv"
clusterFiles <- list.files(clusterDirectory, full.names = TRUE, pattern = clusterFileEnding)

# Read in cluster files
FullClusterData <- clusterFiles %>% 
  map_df(~read_csv(.))


# Find PANTHER GO files, read, combine single dataframe
goTermDirectory <- "AgriGO_Analysis"
goTermFileEnding <- "_cluster.*.tsv"
goTermFiles <- list.files(goTermDirectory, full.names = TRUE, pattern = goTermFileEnding)

# goFullTerm <- goTermFiles %>%
#   map_df(~read_tsv(., col_types = cols(.default = "c")))


agriFileList <- c()

for (file in goTermFiles){
  
  fileRoot <- file %>% 
    str_replace(".*/", "") %>%
    str_replace("\\.tsv", "")
  
  agriFileList <- c(fileRoot, agriFileList)
  # file <- "AgriGO_Analysis/Agrigo_Biological_Process_cluster10.tsv"
  tempFile <- read_tsv(file) %>%
    cbind("FileName"=rep(fileRoot, nrow(.)), .) %>%
    filter(term_type == "P") %>%
    filter(FDR < 0.05) 
  # %>%
  #   separate_rows(entries, sep = "//")
  
    tempFile$entries <- str_replace_all(tempFile$entries, " // ", ",")
    tempFile$entries <- str_replace_all(tempFile$entries, "// ", "")
    tempFile <- filter(tempFile, entries != "")
  assign(fileRoot, tempFile)

  
}

allAgriGOData <- bind_rows(mget(agriFileList))

#View(Agrigo_Biological_Process_cluster2)


allAgriGOData$cluster <- as.numeric(str_extract(allAgriGOData$FileName, "(?<=cluster)\\d+"))



#allAgriGOData[which(allAgriGOData$cluster == 11), "cluster"] <- 3




# Script functions --------------------------------------------------------

# Label functions
genLabelClusterData <- function(clust) {
  # Starts with the initial cluster df 
  # Returns # gene in cluster
  FullClusterData %>% 
    filter(cluster == clust) %>% 
    count() %>% 
    pull()
}

genLabelTermData <- function(clust, mainOntology) {
  # Starts with the gathered ontology cluster df 
  # Returns # of genes with terms
  goGatheredClusterData %>% 
    filter(Ontology == mainOntology) %>%
    filter(cluster == clust) %>% 
    distinct(GeneID) %>% 
    count() %>% 
    pull()
}

genLabelPlotGenes <- function(clust, goGenesCounts) {
  # Starts with the ontology count df and gathered ontology cluster df 
  # Returns # of genes shown in plot
  clusterTerm <- goGenesCounts %>%
    filter(cluster == clust) %>% 
    pull(Term_Description)
  
  goGatheredClusterData %>% 
    filter(cluster == clust) %>% 
    select(GeneID, Term_Description) %>% 
    filter(Term_Description %in% clusterTerm) %>%
    distinct(GeneID) %>%
    count() %>%
    pull()
  
}

genLabelTotalTermData <- function(clust, mainOntology) {
  # Starts with the gathered ontology cluster df 
  # Returns # of terms for cluster
  goGatheredClusterData %>% 
    filter(Ontology == mainOntology) %>% 
    filter(cluster == clust) %>% 
    group_by(Term_Description) %>%
    summarize(n = n()) %>%
    count() %>% 
    pull()
}

genLabelLowTermData <- function(clust, mainOntology) {
  # Starts with the gathered ontology cluster df 
  # Returns # of terms with ', lowTermCutoff, ' gene
  goGatheredClusterData %>% 
    filter(Ontology == mainOntology) %>% 
    filter(cluster == clust) %>% 
    group_by(Term_Description) %>%
    summarize(n = n()) %>%
    filter(n == 1) %>%
    count() %>% 
    pull()
}

# plot output function
plotClustersofInterest <- function(geneData, mainOntology, plotName, titleName){
  
  # p <- ggplot(geneData, 
  #             aes(x = order, y = n)) + 
  #   geom_bar(stat = "identity", position = "dodge") + 
  #   coord_flip() + 
  #   facet_wrap(~clusterLabel, scale = "free_y") +
  #   scale_x_continuous(
  #     breaks = geneData$order,
  #     labels = geneData$Term_Description,
  #     expand = c(0,0)) +
  #   scale_y_continuous(breaks = pretty_breaks()) +
    
   p <- ggplot(geneData, 
           aes(x = queryitem/querytotal, 
               y = order, col = FDR, size = queryitem)) + 
    geom_point() + 
    facet_wrap(~clusterLabel, scale = "free") +
    scale_color_gradient(low = "red", high = "blue") +
    #scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(
      breaks = geneData$order,
      # labels = function(x) str_wrap(geneData$Term, width = 15),
      labels = geneData$Term,
      expand = c(0,0.5)) +
    # facet_wrap_paginate(~clusterLabel, 
    #                     nrow = 1, 
    #                     ncol = 1,  
    #                     scales = 'free', page = i) +
    theme(
      axis.text.y = element_text(color="#993333", size=10, face="bold"),
      axis.text.x = element_text(face="bold"),
      strip.text.x = element_text(size = 9, face="bold")) +
    labs(x = "Gene Ratio (# genes related to GO term / total number of genes in cluster)",
         y = paste("Gene Ontology ", mainOntology, " Complete Terms", sep = ""),
         size = "# of Genes")
  
  z <- ggplotGrob(p)
  
  #  New strip at the top
  z <- gtable_add_rows(z, z$height[7], pos = 6)  # New row added below row 6
  
  # Check the layout
  gtable_show_layout(z)   # New strip goes into row 7 
  # New strip spans columns 5 to 9
  
  z <- gtable_add_grob(z, 
                       list(rectGrob(gp = gpar(col = NA, 
                                               fill = "gray85", 
                                               size = 1)),
                            textGrob(paste("A. thaliana ", 
                                           titleName
                                           , " Gene Entrichment Analysis Via Gene Ontology ", 
                                           mainOntology, 
                                           " Terms\n(FDR cutoff < 0.05, Terms with > 5 genes)", sep = ""), 
                                     gp = gpar(cex = .75, 
                                               fontface = 'bold', 
                                               col = "black"))), 
                       # Adjust the size of the plot title based on number of facets available
                       t=7, l=5, b=7, r=length(levels(geneData$clusterLabel)) * 5, name = c("a", "b"))
  
  # Add small gap between strips - below row 6
  z <- gtable_add_rows(z, unit(2/10, "line"), 7)
  
  for(i in which(grepl("strip-r", z$layout$name))){
    z$grobs[[i]]$layout$clip <- "off"
  }
  
  
  # # Draw it
  # grid.newpage()
  # grid.draw(z)

  filelabel <- paste("Athaliana_Gene_Ontology_",
                     paste(mainOntology, collapse = "_"),
                     "_",
                     plotName,
                     "_Plot.png", 
                     sep = "")
  outfile <- file.path(output_dir, filelabel)

  png(outfile, width = 1920, height = 1144, res = 140)
  # Draw it
  grid.newpage()
  grid.draw(z)
  dev.off()
  
}



# Order GO terms in each cluster. Necessary to display terms on y-axis decending by gene counts
allAgriGOData <- allAgriGOData %>%
  filter(queryitem > 5) %>%
  arrange(cluster, queryitem) %>%
  mutate(order = row_number())

# 12hpc_Upreg_Genes -------------------------------------------------------

mainOntology = "Biological Process"
Upreg12hpcCluster <- filter(allAgriGOData, cluster %in% c(2,4,6))

labels <- c(paste('Cluster #2\n# gene in cluster:', genLabelClusterData(2)),
            paste('Cluster #4\n# gene in cluster:', genLabelClusterData(4)),
            paste('Cluster #6\n# gene in cluster:', genLabelClusterData(6)))

Upreg12hpcCluster$clusterLabel <- factor(Upreg12hpcCluster$cluster, levels = c(2,4,6), labels = labels)

plotClustersofInterest(Upreg12hpcCluster, mainOntology, "12hpc_Upreg", "12hpc Up-regulated")


# 12hpc_Downreg_Genes -------------------------------------------------------

Downreg12hpcCluster <- filter(allAgriGOData, cluster %in% c(9))

labels <- c(paste('Cluster #9\n# gene in cluster:', genLabelClusterData(9)))

Downreg12hpcCluster$clusterLabel <- factor(Downreg12hpcCluster$cluster, levels = c(9), labels = labels)

plotClustersofInterest(Downreg12hpcCluster, mainOntology, "12hpc_Downreg", "12hpc Down-regulated")

# 24hpc_Upnreg_Genes -------------------------------------------------------

Upreg24hpcCluster <- filter(allAgriGOData, cluster %in% c(4,9,10))

labels <- c(paste('Cluster #4\n# gene in cluster:', genLabelClusterData(4)),
            paste('Cluster #9\n# gene in cluster:', genLabelClusterData(9)),
            paste('Cluster #10\n# gene in cluster:', genLabelClusterData(10)))

Upreg24hpcCluster$clusterLabel <- factor(Upreg24hpcCluster$cluster, levels = c(4,9,10), labels = labels)

plotClustersofInterest(Upreg24hpcCluster, mainOntology, "24hpc_Upreg", "24hpc Up-regulated")

# 24hpc_Downreg_Genes -------------------------------------------------------

Upreg24hpcCluster <- filter(allAgriGOData, cluster %in% c(3, 6, 11))

labels <- c(paste('Cluster #3\n# gene in cluster:', genLabelClusterData(3)),
            paste('Cluster #6\n# gene in cluster:', genLabelClusterData(6)),
            paste('Cluster #11\n# gene in cluster:', genLabelClusterData(11)))


Upreg24hpcCluster$clusterLabel <- factor(Upreg24hpcCluster$cluster, levels = c(3,6,11), labels = labels)


plotClustersofInterest(Upreg24hpcCluster, mainOntology, "24hpc_Downreg", "24hpc Down-regulated")



