
# Clear environment
rm(list=ls())


######################################################################

# This script will perform a DEA using DEseq2, utilizing the Wald test for pairwise comparison

######################################################################


# Bioconductor packages
library(DESeq2)
library(BiocParallel)

# CRAN packages
library(readr)
library(tibble)
library(dplyr)
library(ggplot2)
require(graphics)
require(gplots)
library(RColorBrewer)
library(tidyr)
library(scales)
library(pheatmap)
library(Mfuzz)


# Register parallel cores
register(MulticoreParam(workers=3))

### Initate date and time for out directory and log
outDate <- format(Sys.time(), "%m_%d_%Y_%H_%M_%S")


## Set path to gene count files
directory <- "Athaliana_Col0_and_Pen1_2018_HTseq/"
## Percent susept values determined by Bill
Athalian_percent_Susept <- read_csv("AT_ColdSIDR_PercentSucept.csv")

## COR Genes
UpCORGenes <- read_csv("Arabidopsis_Supporting_Data/Shi_et_al_2018_All_Upreg_COR_genes_identified_in_Col.csv")

DownCORGenes <- read_csv("Shi_et_al_2018_All_Downreg_COR_genes_identified_in_Col.csv")

AllCORGenes <- c(UpCORGenes$Gene_ID, DownCORGenes$Gene_ID)

pen1Related <- read_tsv("Arabidopsis_Supporting_Data/PEN1_Associated_Genes/Pen1AssociatedGenes.tsv")

# Outlier Samples determined by iLOO approach; samples had more than 400 flagged genes
#all_outliers <- c("At_Pen1_12hpc_6_E2_T41_counts.txt","At_Pen1_36hpc_7_E3_T5_counts.txt", "At_Col0_24hpc_4_E2_T7_counts.txt","At_Col0_36hpc_1_E1_T4_counts.txt","At_Col0_UTC_6_E2_T43_counts.txt")

pen1Outliers <- c("At_Pen1_UTC_1_E1_T14_counts.txt", "At_Col0_UTC_9_E3_T43_counts.txt", #UTC
                  "At_Col0_24hpc_2_E1_T22_counts.txt", #24hpc
                  "At_Pen1_48hpc_2_E1_T17_counts.txt")

## Retrieve list of sample files and filter based on criteria
sampleFiles <- grep("At",list.files(directory), value = TRUE)

## DEG shrinkage method
shrinkageType <- "apeglm"

## Filter dataset
#sampleFiles <- sampleFiles[grep("Pen1", sampleFiles, invert = TRUE)]
#sampleFiles <- sampleFiles[grep("E1", sampleFiles, invert = TRUE)]
#sampleFiles <- sampleFiles[grep("E2", sampleFiles, invert = TRUE)]
#sampleFiles <- sampleFiles[!sampleFiles %in% all_outliers]
sampleFiles <- sampleFiles[!sampleFiles %in% pen1Outliers]


# Label for output
LOGsampleFiles <- "No_Outliers_Using_PEN1"
outLabel <- paste(LOGsampleFiles, "_with_Intersects_", shrinkageType, "_Shrink", sep = "")

# Directory to which output will be saved:

output_dir <- paste(getwd(), "/Athaliana_Col0_VS_Pen1_2018_DESeq_", outLabel, "_", outDate, sep="")

# Check for and if needed create output directory:
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


### Create log file
logFile <- file.path(output_dir, paste("Analysis_Log_", outLabel, "_", outDate, ".txt", sep = ''))
  cat("At_Col0_Vs_Pen1_DE_Analysis",file=logFile,sep="\n")
  cat(outLabel,file=logFile,append=TRUE)
  cat(outDate,file=logFile,append=TRUE)
  cat("\n\n",file=logFile,append=TRUE)

### Create log for session (package version) info
sessionLog <- file.path(output_dir, paste("Analysis_R_Session_Log_", outLabel, "_", outDate, ".txt", sep = ''))
### Setup session log heading
  cat("Arabidopsis_DE_Analysis",file=sessionLog,sep="\n")
  cat(outDate,file = sessionLog,append=TRUE)
  cat("\n\n",file=sessionLog,append=TRUE)

### Divert R output to the session log file
sink(sessionLog, append = TRUE)
  sessionInfo()
sink()




# ### Update log
# cat("Drop Samples:",file=logFile,sep="\n", append = TRUE)
# cat("Outlier Samples determined by iLOO approach; samples had more than 90 flagged genes",file=logFile,append = TRUE)
# cat(all_outliers,file=logFile,append=TRUE)
# cat("\n\n",file=logFile,append=TRUE)



# Retrieve metadata from file names using Regex ---------------------------

## Metadata example: At_Pen1_48hpc_7_E3_T2_AF_counts.txt

## Organism
sampleOrganism <- sub("At_(.*?)_(.*?)_.*","\\1", sampleFiles)
## Treatment
sampleCondition <- sub("At_(.*?)_(.*?)_.*","\\2", sampleFiles)
## Treatment time
treatmentTime <- as.numeric(sub("At_(.*?)_(.*?)hpc_.*","\\2", sampleFiles))
treatmentTime[is.na(treatmentTime)] <- 0
## Experiment
experiment <-  sub("At_(.*?)_(.*?)_(.*?)_(.*?)_.*","\\4", sampleFiles)
## Tube Number
tubeNumber <-  sub("At_(.*?)_(.*?)_(.*?)_(.*?)_(.*?)_.*","\\4_\\5", sampleFiles)

realtubeNumber <-  sub("At_(.*?)_(.*?)_(.*?)_(.*?)_(.*?)_.*","\\5", sampleFiles)
## Circadian 
time_clock <- ifelse(grepl("At_.*?_(12hpc|36hpc)_.*", sampleFiles),"Night",  grepl("At_.*?_(UTC|24hpc|48hpc)_.*", sampleFiles, "Day", NA))
time_clock[which(time_clock == "TRUE")] <- "Day"


## Create dataframe with sample metadata
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, genotype = sampleOrganism, condition = sampleCondition, tubeNumber = tubeNumber, experiment = experiment, realtubeNumber = realtubeNumber, treatTime = treatmentTime,  time_clock = time_clock)
sampleTable

## Add percent Sucept information to metadata table
sampleTable <- transform(sampleTable, experiment = as.character(experiment))
sampleTable <- left_join(sampleTable, Athalian_percent_Susept, by = c("experiment" = "Experiment", "condition" = "Treatment", "genotype" = "Species"))

## Print out samples separated by genotype and condition 
table(sampleTable[,c("condition", "genotype")])


## Create DESeq2 class variable
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ experiment + genotype + condition)

## Change condition and experiment to factors and relevel
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = c("UTC", "12hpc", "24hpc", "36hpc", "48hpc"))
colData(ddsHTSeq)$experiment <- factor(colData(ddsHTSeq)$experiment, levels = c("E1", "E2", "E3"))
colData(ddsHTSeq)$genotype <- factor(colData(ddsHTSeq)$genotype, levels = c("Col0", "Pen1"))


## Prefiltering based on low gene counts; speeds up analysis, does not impact DEGs
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


### Append experiment breakdown to log file
cat(paste(c("Conditions: ", levels(ddsHTSeq$condition)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Tube Numbers: ", levels(ddsHTSeq$tubeNumber)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Experiment: ", levels(ddsHTSeq$experiment)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Organism: ", levels(ddsHTSeq$organism)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Full Design: ", design(ddsHTSeq)), collapse = " "), sep = "\n", file=logFile,append=TRUE)

cat("\n\n\n", file=logFile,append=TRUE)

### Output log2 fold changes to txt file
### Update filename above based on samples included
sink(logFile, append = TRUE)
  invisible(print(LOGsampleFiles, quote = FALSE))
  invisible(print("", quote = FALSE))
#  invisible(print(paste("Pen1 Outlier Samples Removed: ", paste(pen1Outliers, collapse = " "), sep = " "), quote = FALSE))
  invisible(print("", quote = FALSE))
  invisible(print(paste("Number of Samples: ", length(colnames(assay(ddsHTSeq))), sep = ""), quote = FALSE))
  invisible(print("", quote = FALSE))
  invisible(print(table(sampleTable[,c("condition", "genotype")])))
sink()

## output shrinkage set above to log
cat("\n\n\n", file=logFile,append=TRUE)
cat(paste0("Shrinkage Type: ", shrinkageType, "\n"), file = logFile, append = TRUE)

# Run DESeq analysis ------------------------------------------------------

dds <- DESeq(ddsHTSeq, parallel = TRUE)

## Outlier Analysis
# 
# maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
# idx <- which(rownames(dds) == "ENSMUSG00000076609")
# unname(counts(dds)[idx, ])



## Run QC visualization steps from external file
## Sample PCA, hierarchical clustering, a dendrogram, MA plot...
source("runDESeq2_AT_Sample_QC_Visualization_A.R")


## Creates the list of results names, will be used to retrieve the pairwise comparisons when using Wald Test
nameList <- resultsNames(dds)
nameList <- nameList[!nameList %in% "Intercept"]

timeCourseGeneList <- list()
maxOutlierGeneList <- c()
minOutlierGeneList <- c()


# Output Default Pairwise Analysis ----------------------------------------


## Loop through the default pairwise results and write each to a CSV file
for (i in nameList){
  print(i)
  cat(i, file=logFile, sep = "\n",append=TRUE)

  res <- results(dds, name = i)
   ## Sig genes
   testOut <- as.data.frame(res[which(res$padj < 0.1),]) %>%
     tibble::rownames_to_column("Araport11Genename")

   if(grepl("condition_", i)){
      timeCourseGeneList[[i]][["Araport11Genename"]] <- testOut[,"Araport11Genename"]
      maxOutlierGeneList <- c(testOut[which.max(testOut$log2FoldChange), "Araport11Genename"], maxOutlierGeneList)
      minOutlierGeneList <- c(testOut[which.min(testOut$log2FoldChange), "Araport11Genename"], minOutlierGeneList)
   }
  outfile <- file.path(output_dir, paste(i, "_", outDate,".csv", sep=""))
  write_csv(testOut, outfile)

   ### Output table dimensions to log file
   cat(dim(testOut)[1], sep = "\n", file=logFile,append=TRUE)
   cat("\n\n",file=logFile,append=TRUE)

}


# Create list of 12hpc + 24hpc genes --------------------------------------

sig12hpcGenes <- results(dds, name = "condition_12hpc_vs_UTC", tidy = TRUE) %>%
  filter(padj < 0.1) %>%
  select(row)

sig24hpcGenes <- results(dds, name = "condition_24hpc_vs_UTC", tidy = TRUE) %>%
  filter(padj < 0.1) %>%
  select(row)

intersectGenes <- unlist(intersect(sig12hpcGenes, sig24hpcGenes))


# Create Unshrunken Significant DEG Timecourse ---------------------------------------


## Create a table of the signicant DEGs from the time points to show their change in expression over time.

timeCompareList <- c("condition_12hpc_vs_UTC",
                     "condition_24hpc_vs_UTC",
                     "condition_36hpc_vs_UTC",
                     "condition_48hpc_vs_UTC")

# Loop through the default pairwise results and write each to a CSV file
for (i in timeCompareList){
  print(i)
  assign(i, results(dds, name = i, tidy = TRUE))
}


tConditionTest <- data.frame("GeneID"=as.character(condition_12hpc_vs_UTC$row), 
                             "L2FC12vsUTC" = condition_12hpc_vs_UTC$log2FoldChange, 
                             "L2FC24vsUTC" = condition_24hpc_vs_UTC$log2FoldChange, 
                             "L2FC36vsUTC" = condition_36hpc_vs_UTC$log2FoldChange, 
                             "L2FC48vsUTC" = condition_48hpc_vs_UTC$log2FoldChange)
tConditionTest$GeneID <- as.character(tConditionTest$GeneID)


## Filter genes
tConditionTest <- tConditionTest[which(tConditionTest$GeneID %in% as.vector(unlist(timeCourseGeneList))),]


tGather <- gather(tConditionTest,
                  "L2FC12vsUTC",
                  "L2FC24vsUTC",
                  "L2FC36vsUTC",
                  "L2FC48vsUTC", 
                  key = "sumGroup", value = "Log2FoldChange")

tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))

ggplot(tGather, 
       aes(x = sumGroup, 
           y = Log2FoldChange, 
           group = GeneID)) + 
  geom_line()


# Output unshrunken plot --------------------------------------------------

## Treatments Relative To UTC
outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Unshrunken_Expression_Profiles_TreatmentVSUTC.pdf", sep = ""))

pdf(outfile)
ggplot(tGather, aes(x = sumGroup,
                    y = Log2FoldChange,
                    group = GeneID)) + 
  geom_line() + 
  labs(title = "Unshrunken Expression Profiles for Genes", 
       subtitle = paste("Dataset: ", LOGsampleFiles, "   Sig Genes: ", dim(tConditionTest)[1], sep = ""), 
       x= "Pair-Wise Comparison")
dev.off()

ggplot(tGather, aes(x = sumGroup,
                    y = Log2FoldChange,
                    group = GeneID)) + 
  geom_line() + 
  labs(title = "Unshrunken Expression Profiles for Genes",
       subtitle = paste("Dataset: ", LOGsampleFiles, "   Sig Genes: ", dim(tConditionTest)[1], sep = ""),
       x= "Pair-Wise Comparison")


outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Unshrunken_Expression_Profiles_TreatmentVSUTC.png", sep = ""))
ggsave(outfile, device = "png")


# Output unshrunken intersect plot --------------------------------------------------

## Treatments Relative To UTC

intersecGather <- filter(tGather, GeneID %in% intersectGenes)

ggplot(intersecGather, aes(x = sumGroup,
                           y = Log2FoldChange,
                           group = GeneID)) + 
  geom_line() + 
  labs(title = "Unshrunken Profiles for Intersect Genes",
       subtitle = paste("Dataset: ", LOGsampleFiles, "   Sig Genes: ", length(intersectGenes), sep = ""),
       x= "Pair-Wise Comparison")


outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Unshrunken_Expression_Intersect_Genes_TreatmentVSUTC.png", sep = ""))
ggsave(outfile, device = "png")



rm(list = (timeCompareList))



# Plot counts for gene of interest ----------------------------------------


#Enter gene of interest to plot counts with
#pickGene <- "AT3G11820" # Pen1 loci
pickGene <- "AT2G42540"

d <- plotCounts(dds,
                gene=pickGene,
                intgroup=c("condition", "PercentSuscept", "experiment", "genotype"),
                returnData=TRUE)

dMean <- d %>%
  group_by(condition) %>% 
  summarise(count = mean(count))

dCOL0Mean <- d %>% 
  filter(genotype == "Col0") %>%
  group_by(condition) %>% 
  summarise(count = mean(count))

dPEN1Mean <- d %>% 
  filter(genotype == "Pen1") %>%
  group_by(condition) %>% 
  summarise(count = mean(count))


p <- ggplot() +
  geom_point(data = d, 
             aes(x=condition,
                 y=count,
                 col = experiment,
                 shape = genotype),
             size = 3,
             position=position_jitter(w=0.1,h=0),
             alpha = .6) +
  geom_point(data = dMean,
             aes(x = condition,
                 y = count),
             size = 4,
             shape = 4) +
  geom_point(data = dCOL0Mean,
             aes(x = condition,
                 y = count),
             size = 4,
             shape = 16) +
  geom_point(data = dPEN1Mean,
             aes(x = condition,
                 y = count),
             size = 4,
             shape = 17) +
  
  labs(title = paste("Gene:", pickGene, sep = ""), x= "Sample")

p
outfile <-  file.path(output_dir, paste(LOGsampleFiles, "_", pickGene, "_Gene", ".png", sep = ""))
ggsave(outfile, p, device = "png")




# Unshrunken Max LFC Genes ------------------------------------------------


outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Unshrunken_Max_LFC_Genes", ".pdf", sep = ""))
pdf(outfile, onefile = TRUE)

for (INTGene in maxOutlierGeneList){
  
  d <- plotCounts(dds, gene=INTGene, intgroup=c("condition", "PercentSuscept", "experiment", "genotype"), returnData=TRUE)
  dMean <- d %>% group_by(condition) %>% summarise(count = mean(count))
  
  print(ggplot() +
          geom_point(data = d, 
                     aes(x=condition, 
                         y=count, 
                         col = experiment, 
                         shape = genotype), 
                     size = 3, 
                     position=position_jitter(w=0.1,h=0),
                     alpha = .6) +
          geom_point(data = dMean, 
                     aes(x = condition, 
                         y = count), 
                     size = 4) +
          labs(title = paste("Gene:", INTGene, sep = ""), x= "Sample"))
  

}
dev.off()



# Unshrunken Min LFC Genes ------------------------------------------------


outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Unshrunken_Min_LFC_Genes", ".pdf", sep = ""))
pdf(outfile, onefile = TRUE)

for (INTGene in minOutlierGeneList){
  
  d <- plotCounts(dds, gene=INTGene, intgroup=c("condition", "PercentSuscept", "experiment", "genotype"), returnData=TRUE)
  dMean <- d %>% group_by(condition) %>% summarise(count = mean(count))
  
  print(ggplot() +
    geom_point(data = d, 
               aes(x=condition, 
                   y=count, 
                   col = experiment, 
                   shape = genotype), 
               size = 3, 
               position=position_jitter(w=0.1,h=0),
               alpha = .6) +
    geom_point(data = dMean, 
               aes(x = condition, 
                   y = count), 
               size = 4) +
    labs(title = paste("Gene:", INTGene, sep = ""), x= "Sample"))

}
dev.off()




# Create Shrunken Significant DEG Timecourse ------------------------------

## Create a table of the signicant DEGs from the time points to show their change in expression over time.

timeCompareList <- c("condition_12hpc_vs_UTC",
                     "condition_24hpc_vs_UTC",
                     "condition_36hpc_vs_UTC",
                     "condition_48hpc_vs_UTC",
                     "genotype_Pen1_vs_Col0")


# Loop through the default pairwise results and write each to a CSV file

for (i in timeCompareList){
  print(i)
  
  assign(i, as.data.frame(lfcShrink(dds, coef = i, type = shrinkageType)))
  
  get(i) %>%
    tibble::rownames_to_column("Araport11Genename") %>%
    filter(padj < 0.1) %>%
    write_csv(file.path(output_dir, paste(i, "_Shrunken_Genes_", outDate, ".csv", sep = "")))
  
}


tConditionTest <- data.frame("GeneID" = as.character(row.names(condition_12hpc_vs_UTC)),
                             "L2FC12vsUTC" = condition_12hpc_vs_UTC$log2FoldChange,
                             "L2FC24vsUTC" = condition_24hpc_vs_UTC$log2FoldChange,
                             "L2FC36vsUTC" = condition_36hpc_vs_UTC$log2FoldChange,
                             "L2FC48vsUTC" = condition_48hpc_vs_UTC$log2FoldChange)

tConditionTest$GeneID <- as.character(tConditionTest$GeneID)

tConditionTest <- tConditionTest[which(tConditionTest$GeneID %in% as.vector(unlist(timeCourseGeneList))),]
#tConditionTest <- tConditionTest[as.vector(unlist(timeCourseGeneList)),]

write_csv(tConditionTest, path = file.path(output_dir, paste(outLabel, "_All_Wald_Sig_genes_Combined_shrunken.csv", sep = "")))

inner_join(tConditionTest, pen1Related, by = c("GeneID" = "Araport11Gene")) %>%
  write_csv(path = file.path(output_dir, paste(outLabel, "_Genes_Related_to_Pen1.csv", sep = "")))

tGather <- gather(tConditionTest,
                  "L2FC12vsUTC",
                  "L2FC24vsUTC",
                  "L2FC36vsUTC",
                  "L2FC48vsUTC",
                  key = "sumGroup", value = "Log2FoldChange")
tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))
ggplot(tGather, aes(x = sumGroup, y = Log2FoldChange, group = GeneID)) + geom_line()

## Treatments Relative To UTC
outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Shrunken_Expression_Profiles_TreatmentVSUTC.png", sep = ""))
ggplot(tGather, aes(x = sumGroup, 
                    y = Log2FoldChange, 
                    group = GeneID)) + 
  geom_line() + 
  labs(title = "Shrunken Expression Profiles for Genes", 
       subtitle = paste("Dataset: ", LOGsampleFiles, "   Sig Genes: ", dim(tConditionTest)[1], sep = ""),
       x= "Pair-Wise Comparison")

ggsave(outfile, device = "png")


## Treatments Relative To UTC
outfile <- file.path(output_dir, paste(LOGsampleFiles, "_L2FC_Vs_UTC_Heatmap.png", sep = ""))
pheatmap(tConditionTest[2:5],
         cluster_cols=F, 
         cluster_rows = T, 
         cexCol=0.75, 
         cexRow=0.75, 
         show_rownames = FALSE, 
         main = "Log2 Fold Change of Treatments Relative To UTC", 
         filename= outfile)




## Treatments Relative To UTC

intersecGather <- filter(tGather, GeneID %in% intersectGenes)

write_csv(tConditionTest[which(tConditionTest$GeneID %in% intersectGenes),],
          path = file.path(output_dir, paste(outLabel, "_12hpc_24hpc_Intersect_Genes.csv", sep = "")))


ggplot(intersecGather, aes(x = sumGroup,
                           y = Log2FoldChange,
                           group = GeneID)) + 
  geom_line() + 
  labs(title = "Shrunken Profiles for Intersect Genes",
       subtitle = paste("Dataset: ", LOGsampleFiles, "   Sig Genes: ", length(intersectGenes), sep = ""),
       x= "Pair-Wise Comparison")


outfile <- file.path(output_dir, paste(LOGsampleFiles, "_Shrunken_Expression_Intersect_Genes_TreatmentVSUTC.png", sep = ""))
ggsave(outfile, device = "png")




rm(list = (timeCompareList))




# Cluster Analysis --------------------------------------------------------

set.seed(1) ## Necessary for reproducibiliy of clusters

maxCluster <- 15
minCluster <- 2

averageClusterDisp <- data.frame(cluster = c(minCluster:maxCluster), averageDisp = c(minCluster:maxCluster))
averageClusterDisp$averageDisp <- NA
clusterNameLog <- c()

for (testCluster in c(minCluster:maxCluster))
{
  
  ex.m <- as.matrix(tConditionTest[2:5])
  eset <- new('ExpressionSet', exprs=ex.m)
  kl <- kmeans2(eset, k=testCluster, iter.max=10000)
  
  test <- bind_cols(tConditionTest, cluster=kl$cluster)
  
  tGather <- gather(test,
                    "L2FC12vsUTC",
                    "L2FC24vsUTC",
                    "L2FC36vsUTC",
                    "L2FC48vsUTC",
                    key = "sumGroup", value = "Log2FoldChange")
  
  tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))
  levels(tGather$sumGroup) <- c("12hpc", "24hpc", "36hpc", "48hpc")
  
  clusterNameVar <- paste('cluster', testCluster, '_LFCtableClusters', sep = "")
  assign(clusterNameVar, kl$cluster)
  clusterNameLog <- c(clusterNameVar, clusterNameLog)
  
  averageClusterDisp[which(averageClusterDisp$cluster == testCluster), "averageDisp"] <- tGather %>% 
    group_by(cluster, sumGroup) %>%
    summarise(MeanL2FC=mean(Log2FoldChange),
              STD=(sd(Log2FoldChange))) %>%
    group_by(cluster) %>%
    summarize(MeanclusterSTD = median(STD)) %>%
    summarize(AveSTDForcluster = median(MeanclusterSTD))
  
}  

optCluster <- averageClusterDisp[which.min(averageClusterDisp$averageDisp),'cluster']

clusterKeepName <- paste('cluster', optCluster, '_LFCtableClusters', sep = "")

rm(list = (clusterNameLog[clusterNameLog != clusterKeepName]))


## Full Genes Cluster Analysis

optClusterLFC <- bind_cols(tConditionTest, "cluster" = get(clusterKeepName))


tGather <- gather(optClusterLFC, "L2FC12vsUTC",
                  "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC", key = "sumGroup", value = "Log2FoldChange", -"cluster")
tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))
levels(tGather$sumGroup) <- c("12hpc", "24hpc", "36hpc", "48hpc")

clustersByABSMean <- tGather %>% group_by(cluster, sumGroup) %>%
  summarise(MeanL2FC=mean(Log2FoldChange)) %>%
  filter(sumGroup == "12hpc") %>% 
  mutate(AbsMeanLFC = abs(MeanL2FC)) %>% 
  arrange(desc(AbsMeanLFC))



# Directory to which output will be saved and check for and if needed create output directory:
output_Cluster_dir <- paste(output_dir, "/Clustered_Genes_", LOGsampleFiles,"_", optCluster, "_clusters", sep="")

if(file_test(op='-d',output_Cluster_dir)){
  print(noquote(paste('Output will be saved to the existing directory ',output_Cluster_dir,'.',sep='')))
  print(noquote('___________________________________________________________________________'))
  print(noquote(''))
}else{
  print(noquote(paste('Creating directory',output_Cluster_dir,'for output.')))
  print(noquote('___________________________________________________________________________'))
  print(noquote(''))
  dir.create(output_Cluster_dir)
}

write_tsv(clustersByABSMean, file.path(output_Cluster_dir, paste(LOGsampleFiles, "_Wald_Analysis_Full_Cluster_Analysis_ClustersbyABSMean.tsv", sep = "")))

outfile <- file.path(output_Cluster_dir, paste(LOGsampleFiles, "_Wald_Analysis_Full_Cluster_Analysis_", optCluster, "_Clusters.pdf", sep = ""))
pdf(outfile, onefile = TRUE)
for (cluster in clustersByABSMean$cluster){
  
  
  plotSubtitle <- paste("Wald Analysis, Dataset: ", LOGsampleFiles, ", Sig Genes: ", dim(optClusterLFC[which(optClusterLFC$cluster == cluster),])[1],"/", dim(optClusterLFC)[1], ", ABS Mean L2FC @12hrs: ", round(clustersByABSMean[which(clustersByABSMean$cluster == cluster), "AbsMeanLFC"], 2), sep = "")
  print(ggplot(tGather[which(tGather$cluster == cluster),], aes(x = sumGroup, y = Log2FoldChange, group = GeneID)) +
          geom_line() +
          ylim(min(tGather$Log2FoldChange - 0.2), max(tGather$Log2FoldChange) + 0.2) +
          labs(title = paste("Expression Profiles Cluster #", cluster, sep = ""),
               subtitle = plotSubtitle,
               x= "Pair-Wise Comparison") +
          theme(plot.subtitle=element_text(size=9, hjust=0.5, face="italic", color="black")))
  
  optClusterLFC[which(optClusterLFC$cluster == cluster),] %>%
    write_csv(path = file.path(output_Cluster_dir, paste(outLabel, "_Full_Cluster_Analysis_cluster_", cluster, "_genes.csv", sep = "")))
  
  optClusterLFC[which(optClusterLFC$cluster == cluster),] %>%
    filter(GeneID %in% AllCORGenes) %>%
    write_csv(path = file.path(output_Cluster_dir, paste("COR_", outLabel, "_Full_Cluster_Analysis_cluster_", cluster, "_genes.csv", sep = "")))
  
}
dev.off()

## Output data from analysis so analysis does not need to be re-run
outDATA <- file.path(output_dir, paste("Analysis_RData_", outLabel, "_", outDate, ".RData", sep = ''))
save.image(file = outDATA)

