

# Read in results from Arabidopsis AgriGO analysis clusters, filter, and map to gene names. Outputs results as TSVs.


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

# Read GOC terms for gene names

GOC_names <- read_tsv("AT_Gene_Ontology_Terms/All_Clusters_PantherGenes.tsv", col_types = cols(.default = "c")) %>%
  select(`Mapped IDs`, `Gene Name`)


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
 
  assign(fileRoot, tempFile)
  
  tempFile$entries <- str_replace_all(tempFile$entries, " // ", ",")
  tempFile$entries <- str_replace_all(tempFile$entries, "// ", "")
  
  tempFile <- tempFile %>% separate_rows(entries, sep = ",")
  tempFile <- filter(tempFile, entries != "")
  tempFile <- tempFile %>% 
    left_join(GOC_names, by = c("entries" = "Mapped IDs"))
  
  FullClusterData %>%
    right_join(tempFile, by = c("GeneID" = "entries")) %>%
    write_tsv(file.path(output_dir, paste(fileRoot, "_GeneNames.tsv")))

  
}

allAgriGOData <- bind_rows(mget(agriFileList))

#View(Agrigo_Biological_Process_cluster2)


allAgriGOData$cluster <- as.numeric(str_extract(allAgriGOData$FileName, "(?<=cluster)\\d+"))
