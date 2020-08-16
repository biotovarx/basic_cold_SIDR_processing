
rm(list=ls())


library(DESeq2)
library(ggplot2)
require(graphics)
require(gplots)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
library(tidyverse)
library(gridExtra)


source("Iterative_Leave_One_Out_Approach.R")

fullCountMatrix <- read.csv("Vitis2018FullCountMatrix.csv", row.names = 1, header = TRUE)
#smallCountMatrix <- read.csv("Vitis2018FullCountMatrix.csv", row.names = 1, header = TRUE)
colnames(fullCountMatrix) <- sub("_counts.txt", "", colnames(fullCountMatrix))

samples <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")




# Run iLoo on full matix; 75 samples total
# posOutliers <- iLOO(fullCountMatrix) # Generated 169 gene outliers
# 
# outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
# 
# # Tally the total number of flagged genes per sample
# perOutlierSum <- colSums(table(outliersSamps))
# 
# # Change labels to sample names
# for (name in 1:length(names(perOutlierSum))){
#   #print(names(perOutlierSum)[name])
#   #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#   #print(as.numeric(names(perOutlierSum)[name]))
#   names(perOutlierSum)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])]
# }
# perOutlierSum <- perOutlierSum[order(names(perOutlierSum))]
# perOutlierSum <- as.data.frame(perOutlierSum) %>% tibble::rownames_to_column("SampleNames")
# colnames(perOutlierSum) <- c("SampleNames", "FlaggedGenes")
# perOutlierSum <- separate(perOutlierSum, 'SampleNames', c("Organism", "treatment", "rep", "experiment", "discNumber"), remove = FALSE) %>% unite("experiment", "discNumber", col = "UniqueDisc", remove = FALSE)
# perOutlierSum
# group_by(perOutlierSum[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)





# # Run iLoo on treatment groups; 15 samples per group
# treatments <- c("12hpc", "24hpc", "36hpc", "48hpc", "UTC")
# 
# for (treatment in treatments){
# 
#   assign(paste(treatment, "Matrix", sep = ""), as.data.frame(fullCountMatrix[grep(treatment, names(fullCountMatrix))]))
#   posOutliers <- iLOO(get(paste(treatment, "Matrix", sep = "")))
#   outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
#   treatmentOutput <- paste(treatment, "PerOutlierSum", sep = "")
#   assign(treatmentOutput, colSums(table(outliersSamps)))
#   # Change labels to sample names
#   for (name in 1:length(names(get(treatmentOutput)))){
#     #print(names(perOutlierSum)[name])
#     #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#     #print(as.numeric(names(perOutlierSum)[name]))
#     assign(names(get(treatmentOutput))[name], colnames(posOutliers)[as.numeric(names(get(treatmentOutput))[name])])
#   }
#   assign(treatmentOutput, treatmentOutput[order(names(treatmentOutput))])
#   
#   break
# }


# Run iLoo on treatment groups; 15 samples per group
treatments <- c("12hpc", "24hpc", "36hpc", "48hpc", "UTC")


# ************************************************************
# The following section is very computationally intensive. 
# - Uncomment and run once.
# - The results will be saved and can be loaded if visualization need to be redone.
# - recomment out after.


stop("Don't rerun this section, load your data!!!")

# # Run 12hpctreatment
#   matrix12hpc <- as.data.frame(fullCountMatrix[grep("12hpc", names(fullCountMatrix))])
#   posOutliers <- iLOO(matrix12hpc)
#   rm(matrix12hpc)
#   outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
#   
#   # Tally the total number of flagged genes per sample
#   perOutlierSum12 <- colSums(table(outliersSamps))
#   perOutlierSum12[samples[!samples %in% names(perOutlierSum12)]] <- 0
#   # Change labels to sample names
#   for (name in 1:length(names(perOutlierSum12))){
#     #print(names(perOutlierSum)[name])
#     #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#     #print(as.numeric(names(perOutlierSum)[name]))
#     names(perOutlierSum12)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum12)[name])]
#   }
#   perOutlierSum12 <- perOutlierSum12[order(names(perOutlierSum12))]
#   perOutlierSum12 <- as.data.frame(perOutlierSum12) %>% tibble::rownames_to_column("SampleNames")
#   colnames(perOutlierSum12) <- c("SampleNames", "FlaggedGenes")
#   perOutlierSum12 <- separate(perOutlierSum12, 'SampleNames', c("Organism", "treatment", "rep", "experiment", "discNumber"), remove = FALSE) %>% unite("experiment", "discNumber", col = "UniqueDisc", remove = FALSE)
#   perOutlierSum12
#   #group_by(perOutlierSum12[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
# 
#   
# # Run 24hpctreatment
#   matrix24hpc <- as.data.frame(fullCountMatrix[grep("24hpc", names(fullCountMatrix))])
#   posOutliers <- iLOO(matrix24hpc)
#   rm(matrix24hpc)
#   outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
#   
#   # Tally the total number of flagged genes per sample
#   perOutlierSum24 <- colSums(table(outliersSamps))
#   perOutlierSum24[samples[!samples %in% names(perOutlierSum24)]] <- 0
#   # Change labels to sample names
#   for (name in 1:length(names(perOutlierSum24))){
#     #print(names(perOutlierSum)[name])
#     #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#     #print(as.numeric(names(perOutlierSum)[name]))
#     names(perOutlierSum24)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum24)[name])]
#   }
#   perOutlierSum24 <- perOutlierSum24[order(names(perOutlierSum24))]
#   perOutlierSum24 <- as.data.frame(perOutlierSum24) %>% tibble::rownames_to_column("SampleNames")
#   colnames(perOutlierSum24) <- c("SampleNames", "FlaggedGenes")
#   perOutlierSum24 <- separate(perOutlierSum24, 'SampleNames', c("Organism", "treatment", "rep", "experiment", "discNumber"), remove = FALSE) %>% unite("experiment", "discNumber", col = "UniqueDisc", remove = FALSE)
#   perOutlierSum24
#   #group_by(perOutlierSum24[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
#   
#   
# # Run 36hpctreatment
#   matrix36hpc <- as.data.frame(fullCountMatrix[grep("36hpc", names(fullCountMatrix))])
#   posOutliers <- iLOO(matrix36hpc)
#   rm(matrix36hpc)
#   outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
#   
#   # Tally the total number of flagged genes per sample
#   perOutlierSum36 <- colSums(table(outliersSamps))
#   perOutlierSum36[samples[!samples %in% names(perOutlierSum36)]] <- 0
#   # Change labels to sample names
#   for (name in 1:length(names(perOutlierSum36))){
#     #print(names(perOutlierSum)[name])
#     #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#     #print(as.numeric(names(perOutlierSum)[name]))
#     names(perOutlierSum36)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum36)[name])]
#   }
#   perOutlierSum36 <- perOutlierSum36[order(names(perOutlierSum36))]
#   perOutlierSum36 <- as.data.frame(perOutlierSum36) %>% tibble::rownames_to_column("SampleNames")
#   colnames(perOutlierSum36) <- c("SampleNames", "FlaggedGenes")
#   perOutlierSum36 <- separate(perOutlierSum36, 'SampleNames', c("Organism", "treatment", "rep", "experiment", "discNumber"), remove = FALSE) %>% unite("experiment", "discNumber", col = "UniqueDisc", remove = FALSE)
#   perOutlierSum36
#   #group_by(perOutlierSum36[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
# 
# 
# # Run 48hpctreatment
#   matrix48hpc <- as.data.frame(fullCountMatrix[grep("48hpc", names(fullCountMatrix))])
#   posOutliers <- iLOO(matrix48hpc)
#   rm(matrix48hpc)
#   outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
#   
#   # Tally the total number of flagged genes per sample
#   perOutlierSum48 <- colSums(table(outliersSamps))
#   perOutlierSum48[samples[!samples %in% names(perOutlierSum48)]] <- 0
#   # Change labels to sample names
#   for (name in 1:length(names(perOutlierSum48))){
#     #print(names(perOutlierSum)[name])
#     #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#     #print(as.numeric(names(perOutlierSum)[name]))
#     names(perOutlierSum48)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum48)[name])]
#   }
#   perOutlierSum48 <- perOutlierSum48[order(names(perOutlierSum48))]
#   perOutlierSum48 <- as.data.frame(perOutlierSum48) %>% tibble::rownames_to_column("SampleNames")
#   colnames(perOutlierSum48) <- c("SampleNames", "FlaggedGenes")
#   perOutlierSum48 <- separate(perOutlierSum48, 'SampleNames', c("Organism", "treatment", "rep", "experiment", "discNumber"), remove = FALSE) %>% unite("experiment", "discNumber", col = "UniqueDisc", remove = FALSE)
#   perOutlierSum48
#   #group_by(perOutlierSum48[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
#   
#   # Run UTC
#   matrixUTC <- as.data.frame(fullCountMatrix[grep("UTC", names(fullCountMatrix))])
#   posOutliers <- iLOO(matrixUTC)
#   rm(matrixUTC)
#   outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
#   
#   # Tally the total number of flagged genes per sample
#   perOutlierSumUTC <- colSums(table(outliersSamps))
#   perOutlierSumUTC[samples[!samples %in% names(perOutlierSumUTC)]] <- 0
#   # Change labels to sample names
#   for (name in 1:length(names(perOutlierSumUTC))){
#     #print(names(perOutlierSum)[name])
#     #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
#     #print(as.numeric(names(perOutlierSum)[name]))
#     names(perOutlierSumUTC)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSumUTC)[name])]
#   }
#   perOutlierSumUTC <- perOutlierSumUTC[order(names(perOutlierSumUTC))]
#   perOutlierSumUTC <- as.data.frame(perOutlierSumUTC) %>% 
#     tibble::rownames_to_column("SampleNames")
#   colnames(perOutlierSumUTC) <- c("SampleNames", "FlaggedGenes")
#   perOutlierSumUTC <- separate(perOutlierSumUTC, 'SampleNames', c("Organism", "treatment", "rep", "experiment", "discNumber"), remove = FALSE) %>% unite("experiment", "discNumber", col = "UniqueDisc", remove = FALSE)
#   perOutlierSumUTC
#   #group_by(perOutlierSum36[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
# 
# totalOutlier <- bind_rows(perOutlierSum12, 
#                           perOutlierSum24, 
#                           perOutlierSum36, 
#                           perOutlierSum48, 
#                           perOutlierSumUTC)
# 
# #saveRDS(totalOutlier, "Vitis_Outliers/Vitis_Total_Outliers.rds")


# Run Section -------------------------------------------------------------


totalOutlier <- readRDS("Vitis_Total_Outliers.rds")

subset(totalOutlier, FlaggedGenes > 90) %>% 
  pull(SampleNames) %>%
  sub("(.*)", "\\1_counts.txt", .) %>%

  write.csv("Vitis_2018_Samples_Above_90_Flagged_Genes.csv", row.names = FALSE)

theme_Publication <- function(base_size=14, base_family="serif") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.x = element_text(angle = 90, vjust = 0.3, size = 7),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


# Output TIFF -------------------------------------------------------------


outfile <- "Vitis_2018_Flagged_Genes_By_Samples.tiff"

p <- ggplot(totalOutlier, aes(x = SampleNames, y = FlaggedGenes)) + 
  geom_col() + 
  geom_hline(yintercept=90) +
  labs(x="Sample", 
       y="# of flagged genes", 
       title="Vitis 2018 Flagged Genes by Sample", 
       subtitle="Flagged genes following an iterative leave-one-out approach")
p <- p + scale_colour_Publication()+ theme_Publication()
ggsave(filename = outfile, plot = p, device = "tiff", compression = "lzw")




# Output PDF --------------------------------------------------------------


outfile <- "Vitis_2018_Flagged_Genes_By_Sample_Colored_By_Experiment.pdf"
pdf(outfile)
ggplot(totalOutlier, aes(x = UniqueDisc, 
                         y = FlaggedGenes)) + 
  geom_col(aes(fill = treatment)) + 
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) + 
  labs(x="Sample - Unique Disc", 
       y="# of flagged genes",
       title=expression(paste(italic("V. vinifera"), " Flagged Genes by Leaf Sample")), 
       subtitle="Flagged genes following an iterative leave-one-out approach")
dev.off()


c("Vitis_12hpc_3_E3_D13", "Vitis_12hpc_5_E3_D15", "Vitis_24hpc_5_E3_D15", "Vitis_24hpc_6_E5_D11", "Vitis_36hpc_1_E3_D11", "Vitis_36hpc_9_E5_D14", "Vitis_48hpc_8_E5_D13")
