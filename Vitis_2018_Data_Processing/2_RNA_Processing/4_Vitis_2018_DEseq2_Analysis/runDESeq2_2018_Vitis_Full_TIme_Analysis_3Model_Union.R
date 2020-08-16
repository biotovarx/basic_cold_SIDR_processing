
# # Clear environment
# rm(list=ls())


### This script will generate the exclusive union between the three models


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
library(VennDiagram)


# Register parallel cores
register(MulticoreParam(workers=2))

### Initate date and time for out directory and log
outDate <- format(Sys.time(), "%m_%d_%Y_%H_%M_%S")


# Read in necessary external files ----------------------------------------

## Fennall spreadsheet with Vitis info
Fennall_Grape_data <- read_tsv("Fennall_Grape_Data.tsv")
Fennal_Names <- colnames(Fennall_Grape_data)[2:length(colnames(Fennall_Grape_data))]
## Spreadsheet of real gene names
realGeneName <- read_csv("Vitis_V3Genes_Corresponding_Usable_VIT_Gene_Names.csv")
## Set path to gene count files
directory <- "Vitis_Cold_Count_Data_Files/"
## Percent susept values determined by Bill
vitis_percent_Susept <- read_csv("Vitis_ColdSIDR_PercentSucept.csv")




# Dataset filtering, output directory, log setup --------------------------


# Outlier Samples determined by iLOO approach; samples had more than 90 flagged genes
all_outliers <- c("Vitis_12hpc_3_E3_D13_counts.txt", 
                  "Vitis_12hpc_5_E3_D15_counts.txt", 
                  "Vitis_24hpc_5_E3_D15_counts.txt", 
                  "Vitis_24hpc_6_E5_D11_counts.txt", 
                  "Vitis_36hpc_1_E3_D11_counts.txt", 
                  "Vitis_36hpc_9_E5_D14_counts.txt", 
                  "Vitis_48hpc_8_E5_D13_counts.txt")

#all_outliers <- c("Vitis_12hpc_3_E3_D13_counts.txt", "Vitis_12hpc_5_E3_D15_counts.txt", "Vitis_24hpc_5_E3_D15_counts.txt", "Vitis_36hpc_1_E3_D11_counts.txt", "Vitis_36hpc_9_E5_D14_counts.txt")

## Retrieve list of sample files and filter based on criteria
sampleFiles <- grep("Vitis",list.files(directory), value = TRUE)

sampleFiles <- sampleFiles[grep("48hold", sampleFiles, invert = TRUE)]
#sampleFiles <- sampleFiles[grep("E3", sampleFiles, invert = TRUE)]
#sampleFiles <- sampleFiles[grep("E6", sampleFiles, invert = TRUE)]
#sampleFiles <- sampleFiles[!sampleFiles %in% all_outliers]
#sampleFiles <- sampleFiles[grep("E3_D13|E3_D15", sampleFiles, invert = TRUE, perl = TRUE)]


# Label for output
LOGsampleFiles <- "All_Experiments"
outLabel <- paste(LOGsampleFiles, "_Union_3Models_Combination", sep = "")

# Fennal output
#outname <- c("V3name", "realGenename", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", Fennal_Names)



# Directory to which output will be saved and check for and if needed create output directory:
output_dir <- paste(getwd(), "/Vitis_2018_DESeq_", outLabel, "_", outDate, sep="")

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
logFile <- file.path(output_dir, paste("Analysis_Log_", outDate, ".txt", sep = ''))
cat("Vitis_DE_Analysis",file=logFile,sep="\n")
cat(outDate,file=logFile,append=TRUE)
cat("\n\n",file=logFile,append=TRUE)

### Create log for session (package version) info
sessionLog <- file.path(output_dir, paste("Analysis_R_Session_Log_", outDate, ".txt", sep = ''))
### Setup session log heading
cat("Vitis_DE_Analysis",file=sessionLog,sep="\n")
cat(outDate,file = sessionLog,append=TRUE)
cat("\n\n",file=sessionLog,append=TRUE)

### Divert R output to the session log file
sink(sessionLog, append = TRUE)
sessionInfo()
sink()




### Update log
cat("Drop Samples:",file=logFile,sep="\n", append = TRUE)
cat("Outlier Samples determined by iLOO approach; samples had more than 90 flagged genes",file=logFile,append = TRUE)
cat(all_outliers,file=logFile,append=TRUE)
cat("\n\n",file=logFile,append=TRUE)



# Retrieve metadata from file names using Regex ---------------------------

## Treatment
sampleCondition<- sub("Vitis_(.*?)_.*","\\1", sampleFiles)
## Treatment time
treatmentTime <- as.numeric(sub("Vitis_(.*?)hpc_.*","\\1", sampleFiles))
treatmentTime[is.na(treatmentTime)] <- 0
## Experiment
experiment <-  sub("Vitis_(.*?)_(.*?)_(.*?)_.*","\\3", sampleFiles)
## Generate just disc numbers
trueDiscNumber <-  sub("Vitis_(.*?)_(.*?)_(.*?)_(.*?)_.*","\\4", sampleFiles)
## Generate disc numbers with experiment; UNIQUE
discNumber <-  sub("Vitis_(.*?)_(.*?)_(.*?)_(.*?)_.*","\\3_\\4", sampleFiles)
## Circadian 
time_clock <- ifelse(grepl("Vitis_(12hpc|36hpc)_.*", sampleFiles),"Night",  grepl("Vitis_(UTC|24hpc|48hpc)_.*", sampleFiles, "Day", NA))
time_clock[which(time_clock == "TRUE")] <- "Day"

## Create dataframe with sample metadata
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition, 
                          treatTime = treatmentTime, 
                          treatTime2 = treatmentTime^2, 
                          discNumber = discNumber, 
                          trueDiscNumber = trueDiscNumber, 
                          experiment = experiment, 
                          time_clock = time_clock)

sampleTable

## Add percent Sucept information to metadata table
sampleTable <- transform(sampleTable, experiment = as.character(experiment))
sampleTable <- left_join(sampleTable, vitis_percent_Susept, by = c("experiment" = "Experiment", "condition" = "Treatment"))



## Create DESeq2 class variable

## Design #1: Full: ~ discNumber + condition
##           Reduced: ~ discNumber

## Design #2: Full: ~ discNumber + time_clock + percentSucept
##           Reduced: ~ discNumber + time_clock

## Design #3: Full: ~ discNumber + time1 + time2
##           Reduced: ~ discNumber

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ discNumber + condition)

## Change condition and experiment to factors and relevel
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = c("UTC", "12hpc", "24hpc", "36hpc", "48hpc"))
colData(ddsHTSeq)$experiment <- factor(colData(ddsHTSeq)$experiment, levels = c("E3", "E5", "E6"))

## Prefiltering based on low gene counts; speeds up analysis, does not impact DEGs
keep <- rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq <- ddsHTSeq[keep,]


### Append experiment breakdown to log file
cat(paste(c("Conditions: ", levels(ddsHTSeq$condition)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Disc Numbers: ", levels(ddsHTSeq$discNumber)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Experiment: ", levels(ddsHTSeq$experiment)), collapse = " "), sep = "\n", file=logFile,append=TRUE)
cat(paste(c("Full Design: ", design(ddsHTSeq)), collapse = " "), sep = "\n", file=logFile,append=TRUE)



# Run DESeq analysis for condition ----------------------------------------


cat("Reduced Design: =~ discNumber", sep = "\n", file=logFile,append=TRUE)
ddsCondition <- DESeq(ddsHTSeq, test = "LRT", reduced =~ discNumber, parallel = TRUE)



## Run QC steps from external file 
# ## Sample PCA, hierarchical clustering, a dendrogram, MA plot...
# source("runDESeq2_Vitis_Sample_QC_Visualization.R")



a <- ddsCondition@rowRanges@partitioning@NAMES
z <- results(ddsCondition)$padj
b <- rowMeans(t(t(assays(ddsCondition)$mu[,which(colData(ddsCondition)$condition == "UTC")])/sizeFactors(ddsCondition)[which(colData(ddsCondition)$condition == "UTC")]))
c <- rowMeans(t(t(assays(ddsCondition)$mu[,which(colData(ddsCondition)$condition == "12hpc")])/sizeFactors(ddsCondition)[which(colData(ddsCondition)$condition == "12hpc")]))
d <- rowMeans(t(t(assays(ddsCondition)$mu[,which(colData(ddsCondition)$condition == "24hpc")])/sizeFactors(ddsCondition)[which(colData(ddsCondition)$condition == "24hpc")]))
e <- rowMeans(t(t(assays(ddsCondition)$mu[,which(colData(ddsCondition)$condition == "36hpc")])/sizeFactors(ddsCondition)[which(colData(ddsCondition)$condition == "36hpc")]))
f <- rowMeans(t(t(assays(ddsCondition)$mu[,which(colData(ddsCondition)$condition == "48hpc")])/sizeFactors(ddsCondition)[which(colData(ddsCondition)$condition == "48hpc")]))
## Table "t" contains the expected gene expression levels for each treatment group
## Log2fold changes can be calculated by taking the log2 of one group over the other i.e. log2(qSum12hpc/qSumUTC)
tCondition <- data.frame("GeneID"=as.character(a),
                         "padjValue"=z, 
                         "qMeanUTC"=b, 
                         "qMean12hpc"=c, 
                         "qMean24hpc"=d, 
                         "qMean36hpc"=e, 
                         "qMean48hpc"=f)
#tCondition$LFC2 <- log(t$qMean24hpc/t$qMeanUTC, 2)

tUnionCondition <- tCondition

## Subset to adjusted pvalue < 0.1
tCondition <- tCondition[which(tCondition$padjValue < 0.1),]


conditionLRT <- tCondition$GeneID



# Run DESeq analysis for % suscept ----------------------------------------

design(ddsHTSeq) <- ~discNumber + time_clock + PercentSuscept
cat(paste(c("Full Design: ", design(ddsHTSeq)), collapse = " "), sep = "\n", file=logFile,append=TRUE)

cat("Reduced Design: =~ discNumber + time_clock", sep = "\n", file=logFile,append=TRUE)

ddsSucept <- DESeq(ddsHTSeq, test = "LRT", reduced = ~ discNumber + time_clock, parallel = TRUE)

a <- ddsSucept@rowRanges@partitioning@NAMES
z <- results(ddsSucept)$padj
b <- rowMeans(t(t(assays(ddsSucept)$mu[,which(colData(ddsSucept)$condition == "UTC")])/sizeFactors(ddsSucept)[which(colData(ddsSucept)$condition == "UTC")]))
c <- rowMeans(t(t(assays(ddsSucept)$mu[,which(colData(ddsSucept)$condition == "12hpc")])/sizeFactors(ddsSucept)[which(colData(ddsSucept)$condition == "12hpc")]))
d <- rowMeans(t(t(assays(ddsSucept)$mu[,which(colData(ddsSucept)$condition == "24hpc")])/sizeFactors(ddsSucept)[which(colData(ddsSucept)$condition == "24hpc")]))
e <- rowMeans(t(t(assays(ddsSucept)$mu[,which(colData(ddsSucept)$condition == "36hpc")])/sizeFactors(ddsSucept)[which(colData(ddsSucept)$condition == "36hpc")]))
f <- rowMeans(t(t(assays(ddsSucept)$mu[,which(colData(ddsSucept)$condition == "48hpc")])/sizeFactors(ddsSucept)[which(colData(ddsSucept)$condition == "48hpc")]))

tSucept <- data.frame("GeneID"=as.character(a),
                      "padjValue"=z, 
                      "qMeanUTC"=b, 
                      "qMean12hpc"=c, 
                      "qMean24hpc"=d, 
                      "qMean36hpc"=e, 
                      "qMean48hpc"=f)
# Add in column showing the log2 of max counts over min counts
#t$LFC2 <- log(t$qSum24hpc/t$qSumUTC, 2)

tUnionSucept <- tSucept

tSucept <- tSucept[which(tSucept$padjValue < 0.1),]

percentSuceptLRT <- tSucept$GeneID


# Run DESeq analysis for time^2 -------------------------------------------

#Setup for design formula = ~ discNumber + treatTime + treatTime2
design(ddsHTSeq) <- ~discNumber + treatTime + treatTime2

cat(paste(c("Full Design: ", design(ddsHTSeq)), collapse = " "), sep = "\n", file=logFile,append=TRUE)

cat("Reduced Design: =~ discNumber", sep = "\n", file=logFile,append=TRUE)

ddsTime <- DESeq(ddsHTSeq, test = "LRT", reduced = ~ discNumber, parallel = TRUE)

a <- ddsTime@rowRanges@partitioning@NAMES
z <- results(ddsTime)$padj
b <- rowMeans(t(t(assays(ddsTime)$mu[,which(colData(ddsTime)$condition == "UTC")]) / sizeFactors(ddsTime)[which(colData(ddsTime)$condition == "UTC")]))
c <- rowMeans(t(t(assays(ddsTime)$mu[,which(colData(ddsTime)$condition == "12hpc")])/sizeFactors(ddsTime)[which(colData(ddsTime)$condition == "12hpc")]))
d <- rowMeans(t(t(assays(ddsTime)$mu[,which(colData(ddsTime)$condition == "24hpc")])/sizeFactors(ddsTime)[which(colData(ddsTime)$condition == "24hpc")]))
e <- rowMeans(t(t(assays(ddsTime)$mu[,which(colData(ddsTime)$condition == "36hpc")])/sizeFactors(ddsTime)[which(colData(ddsTime)$condition == "36hpc")]))
f <- rowMeans(t(t(assays(ddsTime)$mu[,which(colData(ddsTime)$condition == "48hpc")])/sizeFactors(ddsTime)[which(colData(ddsTime)$condition == "48hpc")]))


tTime <- data.frame("GeneID"=as.character(a),
                    "padjValue"=z, 
                    "qSumUTC"=b,
                    "qSum12hpc"=c, 
                    "qSum24hpc"=d, 
                    "qSum36hpc"=e, 
                    "qSum48hpc"=f)
# Add in column showing the log2 of max counts over min counts
#t$LFC2 <- log((apply(t[,3:7], 1, max)/apply(t[,3:6], 1, min)), 2)

tUnionTime <- tTime

tTime <- tTime[which(tTime$padjValue < 0.1),]

allLRTTime2 <- tTime$GeneID



#  Examine the plot counts of a single gene for each model ----------------


TriMODgene <- as.character(tUnionCondition[which.min(tUnionCondition$padjValue), "GeneID"])

TriMODgene <- "Vitvi05g00242"

## Plot counts for ddsCondtion
d <- plotCounts(ddsCondition, 
                gene=TriMODgene, 
                intgroup="condition", 
                returnData=TRUE)
dMean <- d %>% 
  group_by(condition) %>% 
  summarise(count = mean(count))

ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  geom_point(alpha = .4) +
  geom_point(data = dMean, size = 4) 
 #+ scale_y_log10(breaks=c(25,100,400))
  


## Plot counts for Percent suscept
d <- plotCounts(ddsSucept, gene=TriMODgene, intgroup=c("condition", "PercentSuscept", "experiment"), returnData=TRUE)

ggplot(d, 
       aes(x=PercentSuscept, 
           y=count, 
           shape = experiment, 
           col = condition)) + 
  geom_point(size = 3) + 
  scale_y_log10(breaks=c(25,100,400))


## Plot counts for treatment time
d <- plotCounts(ddsSucept, 
                gene=TriMODgene, 
                intgroup=c("condition", "PercentSuscept", "experiment", "treatTime"), 
                returnData=TRUE)

ggplot(d, aes(x=treatTime, y=count, shape = experiment, col = condition)) + 
  geom_point(size = 3) + 
  scale_y_log10(breaks=c(25,100,400))

## Plot counts for treatment time ^2
d <- plotCounts(ddsSucept, 
                gene=TriMODgene, 
                intgroup=c("condition", "PercentSuscept", "experiment", "treatTime2"), 
                returnData=TRUE)

ggplot(d, 
       aes(x=treatTime2, 
           y=count, 
           shape = experiment, 
           col = condition)) + 
  geom_point(size = 3) + 
  scale_y_log10(breaks=c(25,100,400))



# Subset only the gene expression levels for each treatment group ---------


tTableUnionCondition <- tUnionCondition[, c(3:7)]



outfile <- file.path(output_dir, paste(LOGsampleFiles, "_VennDiagram_3Model_Intersect.pdf", sep = ""))
pdf(outfile)
venn(list("%Sucept" = percentSuceptLRT, 
          "Treatment" = conditionLRT, 
          "Time^2"=allLRTTime2))
dev.off()

outfile <- file.path(output_dir, paste(LOGsampleFiles, "_VennDiagram_3Model_Intersect.tiff", sep = ""))
venn.diagram(list("Susceptibility" = percentSuceptLRT, "Treatment" = conditionLRT,  "Time" = allLRTTime2),
             filename = outfile,
             cex = 1.55,
             cat.cex = 1.80,
             cat.fontface = "bold",
             imagetype = "tiff",
             fill = c("red", "blue", "green"),
             alpha = 0.20,
             euler.d = TRUE,
             margin = 0.07,
             cat.dist	= 0.075,
             cat.just = list(c(0.5,0.5), c(0.5,0.5), c(0.5, -0.5)),
             category.names = c(
               expression(bold("Susceptibility")),
               expression(bold("Treatment")),
               expression(bold(~ Time^{2}))),
             col = "transparent")




# Retrieve overlapping genes
overlapGenes3Model <- Reduce(union, list(percentSuceptLRT,conditionLRT,allLRTTime2))
#overlapGenes2Model <- Reduce(union, list(percentSuceptLRT,allLRTTime2))

tTableUnionCondition <- tTableUnionCondition[which(rownames(tTableUnionCondition) %in% overlapGenes3Model),]

tTablePersucept <- tCondition[rownames(tSucept),c(3:7)]
tTablePersucept <- tTablePersucept[complete.cases(tTablePersucept),]

## Create new table with columns showing LFC relative to UTC
LFCtableUTC <- tibble("V3name" = row.names(tTableUnionCondition), 
                      "L2FC12vsUTC" = log(tTableUnionCondition$qMean12hpc/tTableUnionCondition$qMeanUTC, 2), 
                      "L2FC24vsUTC" = log(tTableUnionCondition$qMean24hpc/tTableUnionCondition$qMeanUTC, 2), 
                      "L2FC36vsUTC" = log(tTableUnionCondition$qMean36hpc/tTableUnionCondition$qMeanUTC, 2), 
                      "L2FC48vsUTC" = log(tTableUnionCondition$qMean48hpc/tTableUnionCondition$qMeanUTC, 2))

summary(LFCtableUTC)



# Distribution of LFC for sig genes ---------------------------------------

# Get count for each LFC > 2
Above_2 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] >= 2))))
# Get count for each LFC between 1.5 and 2
Between_1.5_and_2 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] >= 1.5 & LFCtableUTC[,column] < 2))))
# Get count for each LFC between 1.5 and 2
Between_1_and_1.5 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] >= 1 & LFCtableUTC[,column] < 1.5))))
# Get count for each LFC between .5 and 1
Between_0.5_and_1 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] >= 0.5 & LFCtableUTC[,column] < 1))))
# Get count for each LFC between 0 and .5
Between_0_and_0.5 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] > 0 & LFCtableUTC[,column] < 0.5))))

Equal_0 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] == 0))))

Below_Neg2 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] <= -2))))
# Get count for each LFC between -1.5 and -2
Between_Neg2_and_Neg1.5 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] <= -1.5 & LFCtableUTC[,column] > -2))))
# Get count for each LFC between -1 and -1.5
Between_Neg1.5_and_Neg1 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] <= -1 & LFCtableUTC[,column] > -1.5))))
# Get count for each LFC between -0.5 and -1
Between_Neg1_and_Neg0.5 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] <= -0.5 & LFCtableUTC[,column] > -1))))
# Get count for each LFC between 0 and -0.5
Between_Neg0.5_and_0 <- unname(sapply(colnames(LFCtableUTC)[2:5], function(column) length(which(LFCtableUTC[,column] < 0 & LFCtableUTC[,column] > -0.5))))


LFCTally_Out <- data.frame(rbind(Below_Neg2, Between_Neg2_and_Neg1.5, Between_Neg1.5_and_Neg1, Between_Neg1_and_Neg0.5, Between_Neg0.5_and_0 , Equal_0, Between_0_and_0.5 ,Between_0.5_and_1, Between_1_and_1.5, Between_1.5_and_2, Above_2))


colnames(LFCTally_Out) <- c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC")
LFCTally_Out

### Write the table to LFC Tally Log file
LFC_Tally_LogFile <- file.path(output_dir, paste(outLabel, ".tsv", sep = ""))

cat(outLabel, file = LFC_Tally_LogFile, sep = "\n")
cat(paste("Number of Genes: ", 
          dim(LFCtableUTC)[1], 
          sep = ""), 
    file = LFC_Tally_LogFile, 
    append = TRUE, 
    sep = "\n")
LFCTally_Out %>% tibble::rownames_to_column("LFC_Range") %>%
  write_tsv(LFC_Tally_LogFile, append = TRUE, col_names = TRUE)



tTable <- as.data.frame(LFCtableUTC[2:5])
row.names(tTable) <- LFCtableUTC[[1]]

## For each row, determine which had the maximum expression, 12hpc, 24hpc, 36hpc, 48hpc 
maxtTable <- sapply(rownames(tTable), function(row) which(tTable[row,] == max(tTable[row,])))

## Subet tTable for genes with max at 24hpc
maxtTable <- tTable[names(maxtTable[which(maxtTable == 2)]),] %>% 
  tibble::rownames_to_column("V3name")
## Log2 of max table
#maxLogtTable <- bind_cols(maxtTable[1], log(maxtTable[2:6],2))

# For each row, determine which had the minimum expression, 12hpc, 24hpc, 36hpc, 48hpc 
mintTable <- sapply(rownames(tTable), function(row) which(tTable[row,] == min(tTable[row,])))
## Subet tTable for genes with min at 24hpc
mintTable <- tTable[names(mintTable[which(mintTable == 2)]),] %>% 
  tibble::rownames_to_column("V3name")
## Log2 of min table
#minLogtTable <- bind_cols(minTable[1], log(mintTable[2:6],2))

## Generate gene set where either max or min expression occured at 24hpc
extreme24hpcGenes <- bind_rows(maxtTable, mintTable)

tRealExtremeG <- right_join(realGeneName, extreme24hpcGenes , by = c("V3annotation" = "V3name")) #%>%
# #   #left_join(Fennall_Grape_data, by = c("realGenename" = "Unique_ID")) 
summary(tRealExtremeG)
write.csv(tRealExtremeG, file.path(output_dir, paste(outLabel, "_Extreme24hpc_DESeqresults.csv", sep = "")), row.names = FALSE)


cat("\n\n\n", file=logFile,append=TRUE)

### Output log2 fold changes to txt file
### Update filename above based on samples included
sink(logFile, append = TRUE)
  print(LOGsampleFiles, quote = FALSE)
  #print(paste("Outlier Samples Removed: ", all_outliers, collapse = " "), quote = FALSE)
  print(paste("Number of Samples: ", length(colnames(assay(ddsCondition))), sep = ""), quote = FALSE)
  print(table(ddsCondition$condition))
  print(summary(LFCtableUTC))
sink()

## Create new table with columns showing LFC relative to previous treatments
LFCtableTimeSeries <- tibble("V3name" = row.names(tTableUnionCondition), 
                             "L2FC12vsUTC" = log(tTableUnionCondition$qMean12hpc/tTableUnionCondition$qMeanUTC, 2), 
                             "L2FC24vs12" = log(tTableUnionCondition$qMean24hpc/tTableUnionCondition$qMean12hpc, 2), 
                             "L2FC36vs24" = log(tTableUnionCondition$qMean36hpc/tTableUnionCondition$qMean24hpc, 2), 
                             "L2FC48vs36" = log(tTableUnionCondition$qMean48hpc/tTableUnionCondition$qMean36hpc, 2))

# t <- LFCtableUTC[which(LFCtableUTC$L2FC48vsUTC > 2), "V3name"]
# # ## Adds real gene names and Fennall data to "t" and writes to CSV file
tRealG <- right_join(realGeneName, LFCtableUTC , by = c("V3annotation" = "V3name")) #%>%
# #   #left_join(Fennall_Grape_data, by = c("realGenename" = "Unique_ID")) 
summary(tRealG)
write.csv(tRealG, file.path(output_dir, paste(outLabel, "_DESeqresults.csv", sep = "")), row.names = FALSE)



## Create new table from "LFCtableUTC" for plotting and plot it; deals with LFC compared to UTC
tGather <- gather(LFCtableUTC, 
                  "L2FC12vsUTC", 
                  "L2FC24vsUTC", 
                  "L2FC36vsUTC", 
                  "L2FC48vsUTC", 
                  key = "sumGroup", 
                  value = "Log2FoldChange")
tGather$sumGroup <- factor(tGather$sumGroup, 
                           levels = c("L2FC12vsUTC", 
                                      "L2FC24vsUTC", 
                                      "L2FC36vsUTC", 
                                      "L2FC48vsUTC"))


outfile <- file.path(output_dir, paste(LOGsampleFiles, "_3Models_Expression_Profiles_TreatmentVSUTC.pdf", sep = ""))
pdf(outfile)
ggplot(tGather, 
       aes(x = sumGroup, 
           y = Log2FoldChange, 
           group = V3name)) + 
  geom_line() + 
  labs(title = "Expression Profiles for Significant Genes Intersect from 3 Models", 
       subtitle = paste("3 Models: %Sucept, Treatment, Time^2     Dataset: ", 
                        LOGsampleFiles, 
                        "   Sig Genes: ", 
                        dim(LFCtableUTC)[1], 
                        sep = ""), 
       x= "Pair-Wise Comparison (Relative to UTC)")
dev.off()


ggplot(tGather, 
       aes(x = sumGroup, 
           y = Log2FoldChange, 
           group = V3name)) + 
  geom_line() + 
  labs(title = "Expression Profiles for Significant Genes Union from 3 Models", 
       subtitle = paste("3 Models: %Sucept, Treatment, Time^2     Dataset: ",
                        LOGsampleFiles,
                        "   Sig Genes: ", 
                        dim(LFCtableUTC)[1], 
                        sep = ""), 
       x= "Pair-Wise Comparison")

outfile <- file.path(output_dir, paste(LOGsampleFiles, "_3Models_Expression_Profiles_TreatmentVSUTC.png", sep = ""))
ggsave(outfile, device = "png")


## Create new table from "LFCtableUTCTimeSeries" for plotting and plot it; deals with LFC compared to the precious treatment group
tTimeGather <- gather(LFCtableTimeSeries, 
                      "L2FC12vsUTC", 
                      "L2FC24vs12", 
                      "L2FC36vs24", 
                      "L2FC48vs36", 
                      key = "sumGroup", 
                      value = "Log2FoldChange")
tTimeGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vs12", "L2FC36vs24", "L2FC48vs36"))
ggplot(tTimeGather, aes(x = sumGroup, y = Log2FoldChange, group = V3name)) + geom_line()


## Generate heatmap showing the correlation between the treatments group
pheatmap(cor(tTableUnionCondition), 
         cluster_cols=F, 
         cluster_rows = F,
         cexCol=0.75, 
         cexRow=0.75, 
         show_rownames = TRUE, 
         main = "Correlation Between Treatment Types")

## Treatments Relative To UTC
outfile <- file.path(output_dir, paste(LOGsampleFiles, "_3Models_L2FC_Vs_UTC_Heatmap.png", sep = ""))
pheatmap(LFCtableUTC[2:5],
         cluster_cols=F, 
         cluster_rows = T, 
         cexCol=0.75, 
         cexRow=0.75, 
         show_rownames = FALSE, 
         main = "Log2 Fold Change of Treatments Relative To UTC",
         filename = outfile)

## Treatments Relative To Previous time point
outfile <- file.path(output_dir, paste(LOGsampleFiles, "_3Models_L2FC_Vs_Previous_Heatmap.png", sep = ""))
pheatmap(LFCtableTimeSeries[2:5], 
         cluster_cols=F, 
         cluster_rows = T,
         cexCol=0.75, 
         cexRow=0.75, 
         show_rownames = FALSE, 
         main = "Log2 Fold Change of Treatments Relative To Treatment", 
         filename = outfile)


## Extremes Analysis
outfile <- file.path(output_dir, paste(LOGsampleFiles,"_3Models_Extreme24hpc_L2FC_Vs_UTC_Heatmap.pdf", sep = ""))
pheatmap(extreme24hpcGenes[2:5], 
         cluster_cols=F, 
         cluster_rows = T,
         cexCol=0.75, 
         cexRow=0.75, 
         show_rownames = FALSE, 
         main = "Extreme 24hpc genes: Log2 Fold Change of Treatments Relative To UTC", 
         filename = outfile)


tGather <- gather(LFCtableUTC, 
                  "L2FC12vsUTC", 
                  "L2FC24vsUTC", 
                  "L2FC36vsUTC", 
                  "L2FC48vsUTC", 
                  key = "sumGroup", 
                  value = "Log2FoldChange")
tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", 
                                                        "L2FC24vsUTC", 
                                                        "L2FC36vsUTC", 
                                                        "L2FC48vsUTC"))

outfile <- file.path(output_dir, paste(LOGsampleFiles, "_3Models_Expression_Profiles_Extreme24hpc_TreatmentVSUTC.pdf", sep = ""))
pdf(outfile)
ggplot(tGather, 
       aes(x = sumGroup, 
           y = Log2FoldChange, 
           group = V3name)) + 
  geom_line() + 
  labs(title = "Expression Profiles for Extreme24hpc Genes Intersect from 3 Models", 
       subtitle = paste("3 Models: %Sucept, Treatment, Time^2     Dataset: ", 
                        LOGsampleFiles,
                        "   Sig Genes: ",
                        dim(extreme24hpcGenes)[1],
                        sep = ""),
       x= "Pair-Wise Comparison")
dev.off()



# Cluster Analysis --------------------------------------------------------

set.seed(1) ## Necessary for reproducibiliy of clusters; DOUBLE CHECK THIS WORKS

maxCluster <- 13
minCluster <- 4

averageClusterDisp <- data.frame(cluster = c(minCluster:maxCluster), 
                                 averageDisp = c(minCluster:maxCluster))
averageClusterDisp$averageDisp <- NA
clusterNameLog <- c()

for (testCluster in c(minCluster:maxCluster))
{

  ex.m <- as.matrix(LFCtableUTC[2:5])
  eset <- new('ExpressionSet', exprs=ex.m)
  kl <- kmeans2(eset ,k=testCluster, iter.max=10000)
  
  test <- bind_cols(LFCtableUTC, cluster=kl$cluster)

  tGather <- gather(test, "L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC", key = "sumGroup", value = "Log2FoldChange")
  
  tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))
  levels(tGather$sumGroup) <- c("12hpc", "24hpc", "36hpc", "48hpc")
  
  clusterNameVar <- paste('cluster', testCluster, '_LFCtableClusters', sep = "")
  assign(clusterNameVar, kl$cluster)
  clusterNameLog <- c(clusterNameVar, clusterNameLog)
  
  averageClusterDisp[which(averageClusterDisp$cluster == testCluster), "averageDisp"] <- tGather %>% group_by(cluster, sumGroup) %>%
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

optClusterLFC <- bind_cols(LFCtableUTC, "cluster" = get(clusterKeepName))


tGather <- gather(optClusterLFC, "L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC", 
                  key = "sumGroup", 
                  value = "Log2FoldChange", -"cluster")
tGather$sumGroup <- factor(tGather$sumGroup, levels = c("L2FC12vsUTC", "L2FC24vsUTC", "L2FC36vsUTC", "L2FC48vsUTC"))
levels(tGather$sumGroup) <- c("12hpc", "24hpc", "36hpc", "48hpc")



clustersByABSMean <- tGather %>% group_by(cluster, sumGroup) %>%
  summarise(MeanL2FC=mean(Log2FoldChange)) %>%
  filter(sumGroup == "24hpc") %>% 
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

write_tsv(clustersByABSMean, file.path(output_Cluster_dir, paste(LOGsampleFiles, "_3Models_Full_Cluster_Analysis_ClustersbyABSMean.tsv", sep = "")))


outfile <- file.path(output_Cluster_dir, 
                     paste(LOGsampleFiles, "_3Models_Full_Cluster_Analysis_", optCluster, "_Clusters.pdf", sep = ""))

# Output all clusters in a single pdf file
pdf(outfile, onefile = TRUE)
for (cluster in clustersByABSMean$cluster){
  
  plotSubtitle <- paste("3 Models, Dataset: ",
                        LOGsampleFiles, ",
                        Sig Genes: ",
                        dim(optClusterLFC[which(optClusterLFC$cluster == cluster),])[1],
                        "/",
                        dim(optClusterLFC)[1],
                        ", ABS Mean L2FC @24hrs: ",
                        round(clustersByABSMean[which(clustersByABSMean$cluster == cluster), "AbsMeanLFC"], 2),
                        sep = "")
  print(ggplot(tGather[which(tGather$cluster == cluster),],
               aes(x = sumGroup,
                   y = Log2FoldChange,
                   group = V3name)) +
          geom_line() +
          ylim(-4.2, 4.2) +
          labs(title = paste("Expression Profiles Cluster #", cluster, sep = ""),
               subtitle = plotSubtitle,
               x= "Pair-Wise Comparison") +
          theme(plot.subtitle=element_text(size=9, 
                                           hjust=0.5, 
                                           face="italic",
                                           color="black")))

  optClusterLFC[which(optClusterLFC$cluster == cluster),] %>% 
    left_join(realGeneName, by = c("V3name" = "V3annotation")) %>%
  write_csv(path = file.path(output_Cluster_dir, 
                             paste(outLabel, "_Full_Cluster_Analysis_cluster_", cluster, "_genes.csv", sep = "")))

}
dev.off()


# Output data from analysis so analysis does not need to be re-run --------

outDATA <- file.path(output_dir, paste("Analysis_RData_", outLabel, "_", outDate, ".RData", sep = ''))
save.image(file = outDATA)


