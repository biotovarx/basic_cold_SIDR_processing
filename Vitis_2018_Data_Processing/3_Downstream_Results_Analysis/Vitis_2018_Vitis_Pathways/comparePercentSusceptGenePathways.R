

#' Investigate the pathways overlaps for vitis signifcant vitis genes found by multiple models.
#' Genes are in the context of the % susceptibility model.


rm(list=ls())

library(tidyverse)
#require(ggforce)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(gplots)
library(gtable)
library(grid)
library(VennDiagram)


# Functions ---------------------------------------------------------------

preProcessPathwayTable <- function(pathTable){
  
  # Filter and clean pathway names
  pathTable %>% 
    filter(`permuted p-value` < 0.05)
}



# Process Files -----------------------------------------------------------


treatmentUnion3Model <- read_tsv("Percent_Suscept_Pathways_Analysis_TSV_Files/VitisPathways_results_1341_Treatment_and_3Model_Union_136_Genes.tsv") %>% 
  preProcessPathwayTable()
treatmentUnion3Model$Pathway <- str_replace(treatmentUnion3Model$Pathway, "vv\\d{1,5}","")

treatmentIntersect <- read_tsv("Percent_Suscept_Pathways_Analysis_TSV_Files/VitisPathways_results_1342_Treatment_Intersect_81_Genes.tsv") %>% 
  preProcessPathwayTable()
treatmentIntersect$Pathway <- str_replace(treatmentIntersect$Pathway, "vv\\d{1,5}","")

threeModelIntersect <- read_tsv("Percent_Suscept_Pathways_Analysis_TSV_Files/VitisPathways_results_1343_3Model_Intersect_55_Genes.tsv") %>%
  preProcessPathwayTable()
threeModelIntersect$Pathway <- str_replace(threeModelIntersect$Pathway, "vv\\d{1,5}","")


outfile <- file.path("ThreeGenesListOverlap.png")
venn.diagram(list("treatmentUnion3Model" = treatmentUnion3Model$Pathway, "Treatment" = treatmentIntersect$Pathway, "3ModelIntersect" = threeModelIntersect$Pathway),
             filename = outfile,
             cex = 1.35,
             cat.cex = 1.30,
             cat.fontface = "bold",
             imagetype = "png",
             fill = c("red", "blue", "green"),
             alpha = 0.20,
             col = "transparent")

