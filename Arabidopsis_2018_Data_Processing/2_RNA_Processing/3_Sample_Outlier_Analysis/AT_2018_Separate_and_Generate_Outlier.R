
rm(list=ls())


library(DESeq2)
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)
require(graphics)
require(gplots)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
library(tidyr)

setwd("Sample_Outlier_Analysis/")

source("Iterative_Leave_One_Out_Approach.R")


ATfullCountMatrix <- read.csv("AT2018FullCountMatrix.csv",
                              row.names = 1,
                              header = TRUE)
#smallCountMatrix <- read.csv("Vitis2018FullCountMatrix.csv", row.names = 1, header = TRUE)
colnames(ATfullCountMatrix) <- sub("_AF_counts.txt", "", 
                                   colnames(ATfullCountMatrix))
colnames(ATfullCountMatrix) <- sub("At_", "", 
                                   colnames(ATfullCountMatrix))

samples <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")

for (geno in c("Col0", "Pen1")){
  
  
  fullCountMatrix <- ATfullCountMatrix[grep(geno, colnames(ATfullCountMatrix))]


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

# Run 12hpctreatment
  matrix12hpc <- as.data.frame(fullCountMatrix[grep("12hpc", names(fullCountMatrix))])
  posOutliers <- iLOO(matrix12hpc)
  rm(matrix12hpc)
  outliersSamps <- as.data.frame(which(!is.na(posOutliers), 
                                       arr.ind = TRUE))
  
  # Tally the total number of flagged genes per sample
  perOutlierSum12 <- colSums(table(outliersSamps))
  perOutlierSum12[samples[!samples %in% names(perOutlierSum12)]] <- 0
  # Change labels to sample names
  for (name in 1:length(names(perOutlierSum12))){
    #print(names(perOutlierSum)[name])
    #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
    #print(as.numeric(names(perOutlierSum)[name]))
    names(perOutlierSum12)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum12)[name])]
  }
  perOutlierSum12 <- perOutlierSum12[order(names(perOutlierSum12))]
  perOutlierSum12 <- as.data.frame(perOutlierSum12) %>%
    tibble::rownames_to_column("SampleNames")
  colnames(perOutlierSum12) <- c("SampleNames", "FlaggedGenes")
  perOutlierSum12 <- separate(perOutlierSum12, 
                              'SampleNames', 
                              c("genotype", "treatment", "rep", "experiment", "tube"), 
                              remove = FALSE) %>% 
    unite("genotype", 
          "treatment", 
          "rep", 
          col = "UniqueSample", 
          remove = FALSE)
  perOutlierSum12
  #group_by(perOutlierSum12[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)

  
# Run 24hpctreatment
  matrix24hpc <- as.data.frame(fullCountMatrix[grep("24hpc", names(fullCountMatrix))])
  posOutliers <- iLOO(matrix24hpc)
  rm(matrix24hpc)
  outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
  
  # Tally the total number of flagged genes per sample
  perOutlierSum24 <- colSums(table(outliersSamps))
  perOutlierSum24[samples[!samples %in% names(perOutlierSum24)]] <- 0
  # Change labels to sample names
  for (name in 1:length(names(perOutlierSum24))){
    #print(names(perOutlierSum)[name])
    #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
    #print(as.numeric(names(perOutlierSum)[name]))
    names(perOutlierSum24)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum24)[name])]
  }
  perOutlierSum24 <- perOutlierSum24[order(names(perOutlierSum24))]
  perOutlierSum24 <- as.data.frame(perOutlierSum24) %>% 
    tibble::rownames_to_column("SampleNames")
  colnames(perOutlierSum24) <- c("SampleNames", "FlaggedGenes")
  perOutlierSum24 <- separate(perOutlierSum24,
                              'SampleNames',
                              c("genotype", "treatment", "rep", "experiment", "tube"),
                              remove = FALSE) %>% 
    unite("genotype", 
          "treatment", 
          "rep", 
          col = "UniqueSample", 
          remove = FALSE)
  perOutlierSum24
  #group_by(perOutlierSum24[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
  
  
# Run 36hpctreatment
  matrix36hpc <- as.data.frame(fullCountMatrix[grep("36hpc", names(fullCountMatrix))])
  posOutliers <- iLOO(matrix36hpc)
  rm(matrix36hpc)
  outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
  
  # Tally the total number of flagged genes per sample
  perOutlierSum36 <- colSums(table(outliersSamps))
  perOutlierSum36[samples[!samples %in% names(perOutlierSum36)]] <- 0
  # Change labels to sample names
  for (name in 1:length(names(perOutlierSum36))){
    #print(names(perOutlierSum)[name])
    #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
    #print(as.numeric(names(perOutlierSum)[name]))
    names(perOutlierSum36)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum36)[name])]
  }
  perOutlierSum36 <- perOutlierSum36[order(names(perOutlierSum36))]
  perOutlierSum36 <- as.data.frame(perOutlierSum36) %>% 
    tibble::rownames_to_column("SampleNames")
  colnames(perOutlierSum36) <- c("SampleNames", "FlaggedGenes")
  perOutlierSum36 <- separate(perOutlierSum36, 
                              'SampleNames', 
                              c("genotype", "treatment", "rep", "experiment", "tube"), 
                              remove = FALSE) %>% 
    unite("genotype", 
          "treatment", 
          "rep", 
          col = "UniqueSample", 
          remove = FALSE)
  perOutlierSum36
  #group_by(perOutlierSum36[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)


# Run 48hpctreatment
  matrix48hpc <- as.data.frame(fullCountMatrix[grep("48hpc", names(fullCountMatrix))])
  posOutliers <- iLOO(matrix48hpc)
  rm(matrix48hpc)
  outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
  
  # Tally the total number of flagged genes per sample
  perOutlierSum48 <- colSums(table(outliersSamps))
  perOutlierSum48[samples[!samples %in% names(perOutlierSum48)]] <- 0
  # Change labels to sample names
  for (name in 1:length(names(perOutlierSum48))){
    #print(names(perOutlierSum)[name])
    #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
    #print(as.numeric(names(perOutlierSum)[name]))
    names(perOutlierSum48)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSum48)[name])]
  }
  perOutlierSum48 <- perOutlierSum48[order(names(perOutlierSum48))]
  perOutlierSum48 <- as.data.frame(perOutlierSum48) %>% 
    tibble::rownames_to_column("SampleNames")
  colnames(perOutlierSum48) <- c("SampleNames", "FlaggedGenes")
  perOutlierSum48 <- separate(perOutlierSum48, 
                              'SampleNames', 
                              c("genotype", "treatment", "rep", "experiment", "tube"), 
                              remove = FALSE) %>% 
    unite("genotype", 
          "treatment", 
          "rep", 
          col = "UniqueSample", 
          remove = FALSE)
  perOutlierSum48
  #group_by(perOutlierSum48[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)
  
  # Run UTC
  matrixUTC <- as.data.frame(fullCountMatrix[grep("UTC", names(fullCountMatrix))])
  posOutliers <- iLOO(matrixUTC)
  rm(matrixUTC)
  outliersSamps <- as.data.frame(which(!is.na(posOutliers), arr.ind = TRUE))
  
  # Tally the total number of flagged genes per sample
  perOutlierSumUTC <- colSums(table(outliersSamps))
  perOutlierSumUTC[samples[!samples %in% names(perOutlierSumUTC)]] <- 0
  # Change labels to sample names
  for (name in 1:length(names(perOutlierSumUTC))){
    #print(names(perOutlierSum)[name])
    #print(colnames(posOutliers)[as.numeric(names(perOutlierSum)[name])])
    #print(as.numeric(names(perOutlierSum)[name]))
    names(perOutlierSumUTC)[name] <- colnames(posOutliers)[as.numeric(names(perOutlierSumUTC)[name])]
  }
  perOutlierSumUTC <- perOutlierSumUTC[order(names(perOutlierSumUTC))]
  perOutlierSumUTC <- as.data.frame(perOutlierSumUTC) %>% 
    tibble::rownames_to_column("SampleNames")
  colnames(perOutlierSumUTC) <- c("SampleNames", "FlaggedGenes")
  perOutlierSumUTC <- separate(perOutlierSumUTC, 
                               'SampleNames', 
                               c("genotype", "treatment", "rep", "experiment", "tube"), 
                               remove = FALSE) %>% 
    unite("genotype", 
          "treatment", 
          "rep", 
          col = "UniqueSample", 
          remove = FALSE)
  perOutlierSumUTC
  #group_by(perOutlierSum48[c(5,8)], UniqueDisc) %>% summarise(sum = sum(FlaggedGenes)) %>% arrange(desc(sum)) %>% separate(UniqueDisc, c("experiment","Disc")) %>% select(experiment, sum) %>% count(experiment)

  # Combine the results for all treatments, filter for samples with flagged genes >90, output to file
  totalOutlier <- bind_rows(perOutlierSum12, perOutlierSum24, perOutlierSum36, perOutlierSum48, perOutlierSumUTC)
  subset(totalOutlier, FlaggedGenes > 400) %>% 
    pull(SampleNames) %>%
    sub("(.*)", "At_\\1_counts.txt", .) %>%
    t() %>%
    write.csv(paste("AT_", geno, "_2018_Samples_Above_400_Flagged_Genes.txt", sep = ""), row.names = FALSE)
    
  assign(paste(geno, "_OutlierSamples", sep = ""), totalOutlier)
  write_tsv(totalOutlier, paste(geno, "_OutlierSamples.tsv", sep = ""))
  
  # Output plots
  outfile <- paste("AT_", geno, "_2018_Flagged_Genes_By_Sample.pdf", sep = "")
  pdf(outfile)
  print(ggplot(totalOutlier, 
             aes(x = UniqueSample, 
                 y = FlaggedGenes)) + 
        geom_col() + 
        geom_hline(yintercept=400) + 
        geom_text(aes(label=FlaggedGenes), 
                  vjust=0) +  
        theme(axis.text.x = element_text(angle = 55, 
                                         hjust = 1)) + 
        labs(x="Sample", 
             y="# of flagged genes", 
             title=paste("AT ", 
                         geno, 
                         " 2018 Flagged Genes by Sample", 
                         sep = ""), 
             subtitle="Flagged genes following an iterative leave-one-out approach"))
  dev.off()

  # Plot the number of flag genes
  outfile <- paste("AT_", geno, "_2018_Flagged_Genes_With_Sample_Colored_By_Treatment.pdf",  sep = "")
  pdf(outfile)
  print(ggplot(totalOutlier, 
             aes(x = UniqueSample, 
                 y = FlaggedGenes)) + 
        geom_col(aes(fill = treatment)) + 
        geom_hline(yintercept=400) + 
        geom_text(aes(label=FlaggedGenes), 
                  vjust=0) +  
        theme(axis.text.x = element_text(angle = 55,
                                         hjust = 1)) + 
        labs(x="Sample", 
             y="# of flagged genes", 
             title=paste("AT ",
                         geno, 
                         " 2018 Flagged Genes by Sample", 
                         sep = ""), 
             subtitle="Flagged genes following an iterative leave-one-out approach"))
  dev.off()


  
}



