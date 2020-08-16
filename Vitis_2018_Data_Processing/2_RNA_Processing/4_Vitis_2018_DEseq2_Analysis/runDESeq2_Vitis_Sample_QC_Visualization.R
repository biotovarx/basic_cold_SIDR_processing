# This script will generate various plots using a DESeqDataSet (dds) object
# This script should be sourced following calling the function DESeq the input matrix object (ddsHTSeq)


# Plots will be saved to "output_dir" which is initiated in the analysis script


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
            axis.text.x = element_text(angle = 90, vjust = 0.3),
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



res <- results(ddsCondition)

vsd <- vst(ddsHTSeq, blind = FALSE)
#rld <- rlog(dds, blind = FALSE)

data <- plotPCA(vsd, intgroup=c("condition", "experiment", "trueDiscNumber", "discNumber"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

outfile <- file.path(output_dir, paste('Vitis_No_Hold_', outLabel,'_PCA_Plot_by_condition_discNumber.pdf', sep =""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, shape=experiment, color = condition )) +
  #geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #facet_wrap(~experiment)+
  geom_text(aes(label=discNumber),vjust="inward",hjust="inward")
#  coord_cartesian(xlim = c(-20,16)) 
dev.off()

outfile <- file.path(output_dir, paste('Vitis_No_Hold_', outLabel,'_PCA_Plot_by_condition_discNumber_experiment.pdf', sep =""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, shape=experiment, color = condition )) +
  #geom_point(size = 2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  facet_wrap(~experiment)+
  geom_text(aes(label=trueDiscNumber),vjust="inward",hjust="inward") +
  theme(legend.position="bottom")
#  coord_cartesian(xlim = c(-20,16)) 
dev.off()


# Color Treatment ---------------------------------------------------------

outfile <- file.path(output_dir, paste('Vitis_No_Hold_', outLabel, '_PCA_Plot_by_condition.pdf', sep = ""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, color=condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()



# Color Experiment --------------------------------------------------------

outfile <- file.path(output_dir, paste('Vitis_No_Hold_', outLabel, '_PCA_Plot_by_experiment.tiff', sep = ""))

p <- ggplot(data, aes(PC1, PC2, color=experiment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title="Vitis 2018 - All Samples - PCA Plot", 
       subtitle="First 2 principal components, colored by experiment")
 

p <- p + theme_Publication() +
  scale_colour_Publication() + scale_color_discrete(name = "Experiment", labels = c("E1", "E2", "E3"))
ggsave(filename = outfile, plot = p, device = "tiff", compression = "lzw")


# Color Disk Number -------------------------------------------------------


outfile <- file.path(output_dir, paste('Vitis_No_Hold_', outLabel, '_PCA_Plot_by_experiment_discNumber.pdf', sep = ""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, shape=experiment, color = discNumber)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

varStabData <- getVarianceStabilizedData(dds)

# Create dendrogram of samples based on similarity
#d <- dist(cor(assay(dds)), method = "euclidean")

sampleDists<- dist(t(assay(vsd)), method = "euclidean")

# Generates heatmap that scales to PDF page
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$experiment, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
outfile <- file.path(output_dir,'Vitis_No_Hold_heatmap.pdf')
pdf(outfile, width=12, height=20)
heatmap(cor(varStabData), cexCol=0.75, cexRow=0.75, scale = "column", margins=c(12,12))
#heatmap.2(cor(vsd), cexCol=0.75, cexRow=0.75, scale = "column", margins=c(12,12))
#pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, height = 12, width = 20)
dev.off()

fit <- hclust(sampleDists, method = "single")

outfile <- file.path(output_dir, paste('Vitis_No_hold_', outLabel, '_dendrogram.pdf', sep = ""))
pdf(outfile, width=18, height=14)
par(mar=c(4,4,4,16))
plot(as.dendrogram(fit), horiz=T)
dev.off()

res <- res[order(res$padj),]
dim(res[which(res$padj < 0.05),])

#png('_DEseq.png')
#plotMA(dds, ylim=c(-2,2), main="DESeq2; experiment + condition")
#dev.off()


# Generate MA plot
# A scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
plotMA(dds, ylim=c(-4,4), main="DESeq2") 

# Plot independent filtering of results
metadata(res)$alpha
metadata(res)$filterThreshold
plot(metadata(res)$filterNumRej,
    type="b", ylab="number of rejections",
    xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)




plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))



