
rm(list=ls())

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

directory <- "Vitis_2018_HTseq/Vitis_Cold_Count_Data_Files_No_Hold/"

filesList <- list.files(directory)
i<-1
allFileCounts <- data.frame()
for (file in filesList){

  # Read in file and rename columns
  countTest <- read.table(paste(directory, file, sep =""))
  names(countTest) <- c("Gene", "Counts")

  # Add up all counts in count column, including "no_feature", "ambiguous", and "alignment_not_unique" summaries at bottom
  countSum <- countTest %>%
    summarise(countSum = sum(Counts))
  

  # Insert file name and read count sum
  countTest <- countTest %>%
    add_row(Gene="File_Name", Counts=file, .after=(length(countTest$Counts)-5)) %>%
    add_row(Gene="Total_Reads", Counts=countSum[[1]], .after=(length(countTest$Counts)-4)) %>%
    tail(n=7)

  countTest$Entry <- i

  allFileCounts <- bind_rows(allFileCounts, countTest)

  i <- i+1

}


# Clean up count table
allFileCounts <- spread(allFileCounts, Gene, Counts)
colnames(allFileCounts) <- gsub("_", "", colnames(allFileCounts))
allFileCounts[,-7] <- sapply(allFileCounts[,-7], as.numeric)
allFileCounts$FileName <- sub("_counts.txt", "", allFileCounts$FileName)


# Create counted and Uncounted features
allFileCounts <- mutate(allFileCounts, TotalCounted= TotalReads - (alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual))
allFileCounts <- mutate(allFileCounts, TotalCounted2 = TotalReads - (alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual))
allFileCounts <- mutate(allFileCounts, Uncounted=alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual)
allFileCounts <- mutate(allFileCounts, Uncounted2=alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual)

# Calculate feature percentages
allFileCounts <- mutate(allFileCounts, AmbigousCount = (ambiguous / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, AlignmentNotUnique = (alignmentnotunique / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, AlignmentNotUnique2 = (alignmentnotunique / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, NoFeature = (nofeature / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, UncountPercent = (Uncounted2 / TotalReads) * 100)


# Plot % of reads counted and Uncounted -----------------------------------


gatheredallAmbigCounts <- gather(allFileCounts, 
                                 NoFeature, 
                                 AlignmentNotUnique, 
                                 AmbigousCount, 
                                 key = "AmbigTotals", value = "PercentTotalReads")


ggplot(gatheredallAmbigCounts, aes(x = reorder(FileName, AlignmentNotUnique2), 
                                   y = PercentTotalReads, 
                                   group = FileName)) + 
  geom_col(aes(fill = AmbigTotals), 
           position = "dodge") + 
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) + 
  labs(x = "Sample", 
       y = "% of Total Reads", 
       title = "Vitis 2018 HTseq Uncounted Reads")



# Plot total reads, total counted, and Uncounted --------------------------


gatheredallFileCounts <- gather(allFileCounts,
                                TotalReads,
                                TotalCounted,
                                Uncounted,
                                key = "Totals", value = "Reads")

gatheredallFileCounts$Totals <- factor(gatheredallFileCounts$Totals, levels = c("TotalReads", "TotalCounted", "Uncounted"))

# Summary stats for total counted
mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"])
min(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"])
max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"])

# Summary stats for Uncounted
mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"])
min(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"])
max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"])


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
            axis.text.y = element_text(size = 8),
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




# Output with total counts and mean count lines  --------------------------



p <- ggplot(gatheredallFileCounts, aes(x = reorder(FileName, TotalCounted2), y = Reads, group = FileName)) + 
  geom_col(aes(fill = Totals), position = "dodge") + 
  geom_hline(yintercept = mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"]), color="blue") +
  geom_hline(yintercept = mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"]), color="green") +
  geom_hline(yintercept = mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalReads"), "Reads"]), color="red") +
  scale_y_continuous(breaks = seq(0, max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalReads"), "Reads"]) + 10000, by = 1000000), labels = unit_format(unit = "m", scale = 1e-6)) + 
  coord_flip() +
  theme(legend.position="bottom") + 
  labs(x="Sample", y="Reads (millions)", title=expression(paste(italic("V. vinifera"), " HTseq Read Count Breakdown")))


p <- p + scale_colour_Publication() + theme_Publication()
p
plotOutFile <- file.path(paste(getwd(), "/Vitis_Outliers", sep = ""), "Vvinifera_2018_HTseq_Read_Count_Breakdown.tiff")
#ggsave(filename = plotOutFile, plot = p, device = "png", width = 7.5, height = 9.5, units = "in")
ggsave(filename = plotOutFile, plot = p, device = "tiff")


# Stack plot of counted and Uncounted -------------------------------------

gatheredallFileCounts <- gather(allFileCounts,
                                TotalReads,
                                TotalCounted,
                                Uncounted,
                                key = "Totals", value = "Reads")

gatheredallFileCounts <- filter(gatheredallFileCounts, Totals %in% c("TotalCounted", "Uncounted"))

gatheredallFileCounts$Totals <- factor(gatheredallFileCounts$Totals, levels = c("TotalCounted", "Uncounted"))

# Summary stats for Uncounted
mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"])
min(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"])
max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "Uncounted"), "Reads"])



# Output only counted and uncounted distributions-------------------------------------------------------------


p <- ggplot(gatheredallFileCounts, aes(x = reorder(FileName, TotalCounted2 + Uncounted2), y = Reads)) + 
  geom_col(aes(fill = Totals), position = "stack") + 
  #scale_fill_manual(values=c("red", "green")) +
  scale_y_continuous(breaks = seq(0, 9000000, by = 1000000), labels = unit_format(unit = "m", scale = 1e-6)) + 
  coord_flip() +
  theme(legend.position="bottom") + 
  scale_fill_discrete(labels=c("Total Counted", "Uncounted")) +
  labs(x="Sample", y="Reads (millions)", title=expression(paste(italic("V. vinifera"), " HTseq Read Count Breakdown")))


p <- p + scale_colour_Publication() + theme_Publication()
plotOutFile <- file.path(paste(getwd(), "/Vitis_Outliers", sep = ""), "Frontier_2018_HTseq_Read_Count_Breakdown.tiff")
ggsave(filename = plotOutFile, plot = p, device = "tiff", width = 7.5, height = 9.5, units = "in", compression = "lzw")
p

#View(gatheredallFileCounts)
