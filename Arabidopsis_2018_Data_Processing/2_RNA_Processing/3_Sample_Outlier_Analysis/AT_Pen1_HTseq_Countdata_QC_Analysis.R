
rm(list=ls())

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

directory <- "Athaliana_Pen1_2018_HTseq/"

filesList <- list.files(directory)
i<-1
allFileCounts <- data.frame()
for (file in filesList){

countTest <- read.table(paste(directory, file, sep =""))
names(countTest) <- c("Gene", "Counts")

countSum <- countTest %>%
  summarise(countSum = sum(Counts))
  


countTest <- countTest %>%
  add_row(Gene="File_Name", Counts=file, .after=(length(countTest$Counts)-5)) %>%
  add_row(Gene="Total_Reads", Counts=countSum[[1]], .after=(length(countTest$Counts)-4)) %>%
  tail(n=7)

countTest$Entry <- i

allFileCounts <- bind_rows(allFileCounts, countTest)

i <- i+1

}

allFileCounts <- spread(allFileCounts, Gene, Counts)

colnames(allFileCounts) <- gsub("_", "", colnames(allFileCounts))

allFileCounts[,-7] <- sapply(allFileCounts[,-7], as.numeric)

allFileCounts$FileName <- sub("_counts.txt", "", allFileCounts$FileName)



allFileCounts <- mutate(allFileCounts, TotalCounted= TotalReads - (alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual))
allFileCounts <- mutate(allFileCounts, TotalCounted2 = TotalReads - (alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual))


allFileCounts <- mutate(allFileCounts, UnCounted=alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual)
allFileCounts <- mutate(allFileCounts, UnCounted2=alignmentnotunique+ambiguous+nofeature+notaligned+toolowaQual)

# Ambiguous percentages
allFileCounts <- mutate(allFileCounts, AmbigousCount = (ambiguous / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, AlignmentNotUnique = (alignmentnotunique / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, AlignmentNotUnique2 = (alignmentnotunique / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, NoFeature = (nofeature / TotalReads) * 100)
allFileCounts <- mutate(allFileCounts, UncountPercent = (UnCounted2 / TotalReads) * 100)


gatheredallAmbigCounts <- gather(allFileCounts,  NoFeature, AlignmentNotUnique, AmbigousCount, key = "AmbigTotals", value = "PercentTotalReads")



ggplot(gatheredallAmbigCounts, aes(x = reorder(FileName, AlignmentNotUnique2), y = PercentTotalReads, group = FileName)) + geom_col(aes(fill = AmbigTotals), position = "dodge") + theme(axis.text.x = element_text(angle = 55, hjust = 1)) + labs(x="Sample", y="% of Total Reads", title="Vitis 2018 HTseq Uncounted Reads")



gatheredallFileCounts <- gather(allFileCounts, TotalReads, TotalCounted, UnCounted, key = "Totals", value = "Reads")


gatheredallFileCounts$Totals <- factor(gatheredallFileCounts$Totals, levels = c("TotalReads", "TotalCounted", "UnCounted"))


# Summary stats for total counted
mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"])
min(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"])
max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"])

# Summary stats for uncounted
mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "UnCounted"), "Reads"])
min(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "UnCounted"), "Reads"])
max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "UnCounted"), "Reads"])


ggplot(gatheredallFileCounts, aes(x = reorder(FileName, TotalCounted2), y = Reads, group = FileName)) + 
  geom_col(aes(fill = Totals), position = "dodge") + 
  geom_hline(yintercept = mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "UnCounted"), "Reads"]), color="blue") +
  geom_hline(yintercept = mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalCounted"), "Reads"]), color="green") +
  geom_hline(yintercept = mean(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalReads"), "Reads"]), color="red") +
  scale_y_continuous(breaks = seq(0, max(gatheredallFileCounts[which(gatheredallFileCounts$Totals == "TotalReads"), "Reads"]) + 10000, by = 1000000), labels = unit_format(unit = "m", scale = 1e-6)) + 
  coord_flip() +
  theme(legend.position="bottom") + 
  labs(x="Sample", y="Reads (million)", title="A. thaliana Pen-1 2018 HTseq Read Breakdown")



plotOutFile <- file.path(getwd(), "Athaliana_Pen1_2018_HTseq_Read_Breakdown.png")
ggsave(plotOutFile, device = "png", width = 6.5, height = 7.5, units = "in")







#View(gatheredallFileCounts)
