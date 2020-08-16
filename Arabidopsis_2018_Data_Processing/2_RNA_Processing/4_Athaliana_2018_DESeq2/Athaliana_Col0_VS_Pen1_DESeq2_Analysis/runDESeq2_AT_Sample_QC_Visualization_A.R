results(dds)
res <- results(dds)

vsd <- vst(dds, blind = FALSE)

data <- plotPCA(vsd, intgroup=c("condition", "experiment", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

outfile <- file.path(output_dir, paste('Athaliana_Col0_', LOGsampleFiles,'_PCA_Plot_by_condition.pdf', sep =""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, color=condition, shape = genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

outfile <- file.path(output_dir, paste('Athaliana_Col0_', LOGsampleFiles, '_PCA_Plot_by_experiment.pdf', sep = ""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, color=experiment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

outfile <- file.path(output_dir, paste('Athaliana_Col0_', LOGsampleFiles, '_PCA_Plot_by_genotype_condition_experiment.pdf', sep = ""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, color=genotype, shape = experiment)) +
  geom_point(size=3) +
  geom_text(aes(label=condition),vjust="inward",hjust="inward")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

outfile <- file.path(output_dir, paste('Athaliana_Col0_', LOGsampleFiles, '_PCA_Plot_by_genotype_condition.pdf', sep = ""))
pdf(outfile)
ggplot(data, aes(PC1, PC2, color=genotype)) +
  #fgeom_point(size=3) +
  geom_text(aes(label=condition),vjust="inward",hjust="inward")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

outfile <- file.path(output_dir, paste0(outLabel, '_PCA_Plot_by_facet_genotype_condition.pdf'))
p <- ggplot(data, aes(PC1, PC2, color=condition, shape = experiment)) +
  geom_point(size = 4) + 
  facet_grid(~genotype) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(strip.text.x = element_text(size = 12))
ggsave(outfile, plot = p, device = "pdf")

outfile <- file.path(output_dir, paste0(outLabel, '_PCA_Plot_by_facet_genotype_condition.png'))
ggsave(outfile, plot = p, device = "png")


outfile <- file.path(output_dir, paste0(outLabel, '_PCA_Plot_by_facet_experiment_genotype_condition.pdf'))
p <- ggplot(data, aes(PC1, PC2, color=condition, shape = genotype)) +
  geom_point(size = 4) + 
  facet_grid(~experiment) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(strip.text.x = element_text(size = 12))
p
ggsave(outfile, plot = p, device = "pdf")

outfile <- file.path(output_dir, paste0(outLabel, '_PCA_Plot_by_facet_experiment_genotype_condition.png'))
ggsave(outfile, plot = p, device = "png")



varStabData <- getVarianceStabilizedData(dds)

# Generates heatmap that scales to PDF page
outfile <- file.path(output_dir,paste('Athaliana_Col0_', LOGsampleFiles, '_heatmap.pdf', sep = ""))
pdf(outfile, width=12, height=14)
heatmap(cor(varStabData), cexCol=0.75, cexRow=0.75, scale = "column", margins=c(12,12))
#heatmap.2(cor(vsd), cexCol=0.75, cexRow=0.75, scale = "column", margins=c(12,12))
dev.off()



# Create dendrogram of samples based on similarity
#d <- dist(cor(assay(dds)), method = "euclidean")

sampleDists<- dist(t(assay(vsd)), method = "euclidean")

# Generates heatmap that scales to PDF page
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$experiment, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# outfile <- file.path(output_dir,'Vitis_No_Hold_heatmap.pdf')
pdf(outfile, width=12, height=20)
heatmap(cor(varStabData), cexCol=0.75, cexRow=0.75, scale = "column", margins=c(12,12))
#heatmap.2(cor(vsd), cexCol=0.75, cexRow=0.75, scale = "column", margins=c(12,12))
#pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, height = 12, width = 20)

dev.off()

d <- dist(cor(varStabData), method = "euclidean")

sampleDists<- dist(t(assay(vsd)), method = "euclidean")

fit <- hclust(sampleDists, method = "single")


outfile <- file.path(output_dir,'Athaliana_Col0_dendrogram.pdf')
pdf(outfile, width=18, height=14)
par(mar=c(4,4,4,16))
plot(as.dendrogram(fit), horiz=T)
dev.off()

res <- res[order(res$padj),]
dim(res[which(res$padj < 0.05),])

#png('_DEseq.png')
#plotMA(dds, ylim=c(-2,2), main="DESeq2; experiment + condition")
#dev.off()

plotMA(dds, ylim=c(-2,2), main="DESeq2")
