# Read in genes vitis analysis clusters and create some nice plot


# Clear environment
rm(list=ls())

# CRAN packages
library(tidyverse)
library(ggforce)
library(gtable)  
library(grid)
library(scales)

# Directory to which output will be saved and check for and if needed create output directory:
output_dir <- paste(getwd(), "/Vitis_Analysis_Clusters_Processing", sep="")



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


# Read and process GO files ------------------------------------------------

# Find cluster files, read, combine single dataframe
clusterDirectory <- "Clustered_Genes_No_E3_and_No_Outliers_13_clusters/"
clusterFileEnding <- "_genes.csv"
clusterFiles <- list.files(clusterDirectory, full.names = TRUE, pattern = clusterFileEnding)

# Read and process cluster files
# Cluster files have functional information from VitisNet
FullClusterData <- clusterFiles %>% 
  map_df(~read_csv(., col_types = cols(.default = "c")))
# Limit to cluster data to genes with GO data, make the data tidy
slimGenes <- FullClusterData %>%
  select(realGenename, cluster, Gene_Ontology) %>%
  separate_rows(Gene_Ontology, sep = "\\|") %>%
  filter(Gene_Ontology != "0") %>%
  drop_na(Gene_Ontology) %>%
  separate(Gene_Ontology, into = c("GO_Term", "Term_Description", "Ontology"), sep = ";")


# Read and Process GO terms from PANTHER
vitisPantherGOTerms <- read_tsv("VviniferapantherGeneList2_25_19.tsv")
vitisGatheredPantherGo <- vitisPantherGOTerms %>%
  gather(GO_database_MF_Complete, GO_database_BP_Complete, GO_database_CC_Complete, key = "Ontology", value = "GO_Term") %>%
  separate_rows(GO_Term, sep = ";") %>%
  filter(!is.na(GO_Term)) %>%
  separate(GO_Term, into = c("Term_Description", "GO_Term"), sep = "\\(GO") %>%
  select(Mapped_IDs, Ontology, Term_Description, GO_Term)


# Clean up GO data ---------------------------------------------------------

# Change VitisNet ontology names
slimGenes$Ontology <- ifelse(slimGenes$Ontology == "F", 
                             "Molecular Function", 
                             ifelse(slimGenes$Ontology == "P", 
                                    "Biological Process", 
                                    ifelse(slimGenes$Ontology == "C", 
                                           "Cellular Component", 
                                           NA
                                           )
                                    )
                             )

# Change PANTHER ontology names
vitisGatheredPantherGo$Ontology <- ifelse(vitisGatheredPantherGo$Ontology == "GO_database_MF_Complete",
                             "Molecular Function",
                             ifelse(vitisGatheredPantherGo$Ontology == "GO_database_BP_Complete",
                                    "Biological Process",
                                    ifelse(vitisGatheredPantherGo$Ontology == "GO_database_CC_Complete",
                                           "Cellular Component",
                                           NA)
                                    )
                             )
vitisGatheredPantherGo <- FullClusterData %>%
  select(realGenename, cluster) %>%
  right_join(vitisGatheredPantherGo, by = c("realGenename" = "Mapped_IDs"))

# Cleanup PANTHER names
vitisGatheredPantherGo$GO_Term <- str_replace(vitisGatheredPantherGo$GO_Term, "^:", "GO:")
vitisGatheredPantherGo$GO_Term <- str_replace(vitisGatheredPantherGo$GO_Term, "\\)$", "")

# Output gene with blank VitisNet Ontologies, use as imput to PANTHER
setdiff(FullClusterData$realGenename, slimGenes$realGenename) %>%
  as.data.frame() %>%
  write_csv(file.path(output_dir, "Vvinifera_Vitisnet_No_GO_Terms_Updated.csv"))

# Combine the two tables (VitisNet + PANTHER) into single dataset
slimGenes <- bind_rows(vitisGatheredPantherGo, slimGenes)



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
  slimGenes %>% 
    filter(Ontology == mainOntology) %>%
    filter(cluster == clust) %>% 
    distinct(realGenename) %>% 
    count() %>% 
    pull()
}

genLabelPlotGenes <- function(clust, slimGenesCounts) {
# Starts with the ontology count df and gathered ontology cluster df 
# Returns # of genes shown in plot
  clusterTerm <- slimGenesCounts %>%
    filter(cluster == clust) %>% 
    pull(Term_Description)
  
  slimGenes %>% 
    filter(cluster == clust) %>% 
    select(realGenename, Term_Description) %>% 
    filter(Term_Description %in% clusterTerm) %>%
    distinct(realGenename) %>%
    count() %>%
    pull()
  
}

genLabelTotalTermData <- function(clust, mainOntology) {
# Starts with the gathered ontology cluster df 
# Returns # of terms for cluster
  slimGenes %>% 
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
  slimGenes %>% 
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
  
  p <- ggplot(geneData, 
              aes(x = order, y = n)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    coord_flip() + 
    facet_wrap(~clusterLabel, scale = "free_y") +
    scale_x_continuous(
      breaks = geneData$order,
      labels = geneData$Term_Description,
      expand = c(0,0)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    # facet_wrap_paginate(~clusterLabel, 
    #                     nrow = 1, 
    #                     ncol = 1,  
    #                     scales = 'free', page = i) +
    theme(
      axis.text.y = element_text(color="#993333", size=12, face="bold"),
      axis.text.x = element_text(face="bold"),
      strip.text.x = element_text(size = 9, face="bold")) +
    labs(x = paste("Gene Ontology ", mainOntology, " Complete Terms", sep = ""),
         y = "Number of Genes Associated with Terms")
  
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
          textGrob(paste("V. vinifera ", titleName, " Functional Annotation Via Gene Ontology ", mainOntology, " Terms\n(Limited to Terms with 2 or more genes)", sep = ""), 
              gp = gpar(cex = 1, 
              fontface = 'bold', 
              col = "black"))), 
# Adjust the size of the plot title based on number of facets available
            t=7, l=5, b=7, r=length(levels(geneData$clusterLabel)) * 5, name = c("a", "b"))
  
  # Add small gap between strips - below row 6
  z <- gtable_add_rows(z, unit(2/10, "line"), 7)
  
  for(i in which(grepl("strip-r", z$layout$name))){
    z$grobs[[i]]$layout$clip <- "off"
  }
  

  filelabel <- paste("Vvinifera_Gene_Ontology_",
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

getGOTermsClustersofInterest <- function(mainOntology){
  
  slimGenesCounts <- slimGenes %>%
  filter(Ontology == mainOntology) %>%
  #  filter(Term_Description != "biological_process") %>%
  group_by(cluster, Term_Description) %>%
  summarise(n = n()) %>%
  ungroup %>%
  filter(n >= 2) %>%
  arrange(cluster, n) %>%
  mutate(order = row_number())

  
  # 24hpc_Upreg_Genes -------------------------------------------------------
  
  intCluster <- 3
  Upreg24hpcCounts <- filter(slimGenesCounts, cluster %in% c(intCluster))
  
  
  labels <- c(paste('Cluster #3\n# gene in cluster:', genLabelClusterData(intCluster),
                    '\n# of genes with terms:', genLabelTermData(intCluster, mainOntology),
                    '\n# of genes shown in plot:', genLabelPlotGenes(intCluster, slimGenesCounts),
                    '\n# of terms for cluster:', genLabelTotalTermData(intCluster, mainOntology),
                    '\n# of terms with ', lowTermCutoff, ' gene:', genLabelLowTermData(intCluster, mainOntology)))
  
  Upreg24hpcCounts$clusterLabel <- factor(Upreg24hpcCounts$cluster, levels = c(intCluster), labels = labels)
  
  plotClustersofInterest(Upreg24hpcCounts, mainOntology, "24hpc_Upreg", "24hpc Up-regulated")
  
  # 24hpc_Downreg_Genes -------------------------------------------------------
  
  intCluster <- 2
  Downreg24hpcCounts <- filter(slimGenesCounts, cluster %in% c(intCluster))
  
  labels <- c(paste('Cluster #2\n# gene in cluster:', genLabelClusterData(intCluster),
                    '\n# of genes with terms:', genLabelTermData(intCluster, mainOntology),
                    '\n# of genes shown in plot:', genLabelPlotGenes(intCluster, slimGenesCounts),
                    '\n# of terms for cluster:', genLabelTotalTermData(intCluster, mainOntology),
                    '\n# of terms with ', lowTermCutoff, ' gene:', genLabelLowTermData(intCluster, mainOntology)))
  
  Downreg24hpcCounts$clusterLabel <- factor(Downreg24hpcCounts$cluster, levels = c(intCluster), labels = labels)
  
  plotClustersofInterest(Downreg24hpcCounts, mainOntology, "24hpc_Downreg", "24hpc Down-regulated")
  
  # Max_Upreg_Genes -------------------------------------------------------

  intCluster <- 13
  MaxUpregCounts <- filter(slimGenesCounts, cluster %in% c(intCluster))

  labels <- c(paste('Cluster #13\n# gene in cluster:', genLabelClusterData(intCluster),
                    '\n# of genes with terms:', genLabelTermData(intCluster, mainOntology),
                    '\n# of genes shown in plot:', genLabelPlotGenes(intCluster, slimGenesCounts),
                    '\n# of terms for cluster:', genLabelTotalTermData(intCluster, mainOntology),
                    '\n# of terms with ', lowTermCutoff, ' gene:', genLabelLowTermData(intCluster, mainOntology)))


  MaxUpregCounts$clusterLabel <- factor(MaxUpregCounts$cluster, levels = c(intCluster), labels = labels)


  plotClustersofInterest(MaxUpregCounts, mainOntology, "Max_Upreg", "Max L2FC Up-regulated")

}


# Run Processing ----------------------------------------------------------

# Set variables and process clusters of interest for main ontologies

lowTermCutoff <- "1"

mainOntologies <- c("Biological Process", "Molecular Function","Cellular Component")

map(mainOntologies,getGOTermsClustersofInterest)

  
  
#   
# extractUPregGenesbyOntology <- function(genes, ontology){
#   genes %>% 
#     filter(L2FC24vsUTC > 0) %>% 
#     filter(Ontology == ontology)
#   
# }
# 
# for (i in c("Biological_Process", "Cellular_Component", "Molecular_Function")){
#   
#   assign(i, extractUPregGenesbyOntology(slimGenes, i))
#   
# }
# 
# extractDOWNregGenesbyOntology <- function(genes, ontology){
#   genes %>% 
#     filter(L2FC24vsUTC < 0) %>% 
#     filter(Ontology == ontology)
#   
# }
# 
# for (i in c("Biological_Process", "Cellular_Component", "Molecular_Function")){
#   
#   assign(i, extractDOWNregGenesbyOntology(slimGenes, i))
#   
# }
  
  
  
# 
# 
# 
# 
# 
# 
# clustersGathered <- gather(FullClusterData, "L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC", key = "pairGroup", value = "Log2FoldChange", -"cluster")
# 
# clustersGathered$pairGroup <- factor(clustersGathered$pairGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))
# levels(clustersGathered$pairGroup) <- c("12hpc", "24hpc", "36hpc", "48hpc")
# 
# labels <- c(paste('Cluster #1\nGenes:', FullClusterData %>% filter(cluster == 1) %>% count() %>% pull()),
#             paste('Cluster #2\nGenes:', FullClusterData %>% filter(cluster == 2) %>% count() %>% pull()),
#             paste('Cluster #3\nGenes:', FullClusterData %>% filter(cluster == 3) %>% count() %>% pull()),
#             paste('Cluster #4\nGenes:', FullClusterData %>% filter(cluster == 4) %>% count() %>% pull()),
#             paste('Cluster #5\nGenes:', FullClusterData %>% filter(cluster == 5) %>% count() %>% pull()),
#             paste('Cluster #6\nGenes:', FullClusterData %>% filter(cluster == 6) %>% count() %>% pull()),
#             paste('Cluster #7\nGenes:', FullClusterData %>% filter(cluster == 7) %>% count() %>% pull()),
#             paste('Cluster #8\nGenes:', FullClusterData %>% filter(cluster == 8) %>% count() %>% pull()),
#             paste('Cluster #9\nGenes:', FullClusterData %>% filter(cluster == 9) %>% count() %>% pull()),
#             paste('Cluster #10\nGenes:', FullClusterData %>% filter(cluster == 10) %>% count() %>% pull()),
#             paste('Cluster #11\nGenes:', FullClusterData %>% filter(cluster == 11) %>% count() %>% pull()),
#             paste('Cluster #12\nGenes:', FullClusterData %>% filter(cluster == 12) %>% count() %>% pull()),
#             paste('Cluster #13\nGenes:', FullClusterData %>% filter(cluster == 13) %>% count() %>% pull()))
# 
# clustersGathered$cluster <- factor(clustersGathered$cluster, levels = c(1:13), labels = labels)
# 
# 
# 
# 
# 
# ggplot(clustersGathered, aes(x = pairGroup, y = Log2FoldChange, group = V3name)) +
#   geom_line() +
#   ylim(-4.2, 4.2) + 
#   theme(axis.text.x = element_text(angle = 65, hjust = 1)) +
#   facet_wrap(~cluster) +
#   labs(title = "Grape Gene Clusters") + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   xlab("Treatment Group Relative to UTC") +
#   ylab("Log2 Fold Change")
# 
# 
# 
# 
# # 
# # +
# #   labs(title = paste("Expression Profiles Cluster #", cluster, sep = ""),
# #        subtitle = plotSubtitle,
# #        x= "Pair-Wise Comparison") +
# #   theme(plot.subtitle=element_text(size=9, hjust=0.5, face="italic", color="black")
# 
# 
# 

