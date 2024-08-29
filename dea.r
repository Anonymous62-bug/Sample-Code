#Run differential expression analysis between groups (pairwise)
#Usage: Rscript dea.r comparisons.tsv counts.tsv sampleinfo.tsv
#All files need to be TSV

library(tidyr)
library(edgeR)
library(readr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#Plotting PCA with batch and group information
plot_pca <- function(dge,dir_path){
	counts <- data.frame(cpm(dge, normalized.lib.sizes = T, log = T), check.names = F)
	pc <- prcomp(counts)
	print(summary(counts))
	pca_summary <- summary(counts)
	prop_variance <- summary(test1)$importance[2,]*100
	coord <- data.frame(pca_summary$rotation)[1:2]
	coord$group <- dge$samples[rownames(coord),"group"]
	coord$batch <- dge$samples[rownames(coord),"batch"]

	pdf(paste0(dir_path,"/pca.pdf"))
	print(ggplot(coord, aes(x=PC1, y=PC2, shape=batch, color=group, fill=group)) +
		geom_point(size = 1.5) +
		scale_color_brewer(palette="Dark2") +
		geom_text_repel(label = rownames(coord), vjust="inward",hjust="inward") +
		xlab(paste0("PC1 captures ",prop_variance[1],"% of variance")) +
		ylab(paste0("PC2 captures ",prop_variance[2],"% of variance")) +
		theme_light())
	dev.off()

#Volcano Plot of Genes Colouring Significant Genes
plot_volcano <- function(dea,dir_path){
	dea$type <- "Not Significant"
	pval <- 0.05
	fc <- [0.5,-0.5,1,-1]
	dea[dea$FDR < pval & dea$logFC > fc[1],"type"] <- "logFC > as.character(fc[1])"
	dea[dea$FDR < pval & dea$logFC < fc[2],"type"] <- "logFC < as.character(fc[2])"
	dea[dea$FDR < pval & dea$logFC > fc[3],"type"] <- "logFC > as.character(fc[3])"
	dea[dea$FDR < pval & dea$logFC < fc[4],"type"] <- "logFC < as.character(fc[4])"
	dea$logFDR <- -log10(dea$FDR)
	pdf(paste0(dir_path,"/volcano.pdf"))
	print(ggplot(dea, aes(x=logFC,y=logFDR, color=type)) +
	geom_point(size=0.5) +
	scale_color_brewer(palette="Dark2") +
	ylab("-log10(FDR)"))
	dev.off()
}

#File1: args[1]; Info: Provide case & control samples to be compared
#Format:
#case	control
#sample1	sample2

#File2: args[2]; Info: Raw counts with column 1 containing gene id/name
#Format:
#geneid	sample1	sample2	sample3
#geneA	26	32	1

#File3: args[2]; Info: Containing any sample metadata available
#Format:
#sample	group	batch (add other columns as needed)
#sample1	GroupA	batch1
#sample2	GroupB	batch1

comparison <- read_tsv(args[1], col_names=T) %>%
	data.frame()
counts <- read_tsv(args[2], col_names=T) %>%
	column_to_rownames(var = "geneid") %>%
	as.matrix()
counts <- tail(counts, -4)
samples <- read_tsv(args[3], col_names=T) %>%
	data.frame()
rownames(samples) <- samples$sample

for(i in c(1:nrow(comparison))){
	comp_samples <- samples[samples$group %in% comparison[i,],]
	comp_counts <- counts[,comp_samples$sample]
	comp_group <- comp_samples$group
	has_replicates <- ifelse(min(table(comp_samples$group)) > 1,TRUE,FALSE)

	#Checking if both case and control have multiple batches
	has_batch <- data.frame(name=unique(comp_samples$batch))
	has_batch$case <- unique(comp_samples$batch) %in% comp_samples$batch[comp_samples$group == comparison[i,"case"]]
	has_batch$control <- unique(comp_samples$batch) %in% comp_samples$batch[comp_samples$group == comparison[i,"control"]]
	has_batch$both <- has_batch$case && has_batch$control
	has_batches <- ifelse(sum(has_batch$both) >=2 ,TRUE,FALSE)
	dge <- DGEList(counts = comp_counts, group = comp_group,samples = comp_samples)
	dir_path <- paste0(comparison[i,"case"],"-vs-",comparison[i,"control"])
	print(dir_path)
	dir.create(dir_path, showWarnings = FALSE)

	#Filter & Normalize counts matrix
	keep <- filterByExpr(dge)
	dge <- dge[keep,,keep.lib.sizes=FALSE]
	dge <- calcNormFactors(dge)
	plot_pca(dge,dir_path)

	#Experimental design & DEA model
	if(!has_replicates){
		model <- exactTest(dge, dispersion=0.2^2)
	}else{
		if(has_batches){
			design <- model.matrix(~0 + comp_group + comp_samples$batch)
		}else{
			design <- model.matrix(~0 + comp_group)
		}
		write_tsv(as.data.frame(design),paste0(dir_path,"/design.tsv"), col_names=TRUE)
		index <- which(colnames(design) == paste0("comp_group",comparison[i,"control"]))
		if(index == 1){
			contrasts <- c(c(-1,1),rep(x=0, times=ncol(design)-2))
		}else{
			contrasts <- c(c(1,-1),rep(x=0, times=ncol(design)-2))
		}
		dge <- estimateDisp(dge, design)
		pdf(paste0(dir_path,"/bcv.pdf"))
		plotBCV(dge)
		dev.off()
		fit <- glmQLFit(dge, design)
		model <- glmQLFTest(fit, contrast=contrasts)
		pdf(paste0(dir_path,"/MD.pdf"))
		plotMD(model)
		abline(h=c(-1,1), col="blue")
		dev.off()
	}

	#Write Output
	dea <- rownames_to_column(as.data.frame(topTags(model,n=model$model$table)),"geneid")
	plot_volcano(dea,dir_path)
	gtf <- read_tsv("/path/to/id2gene.tsv", col_names=TRUE)
	dea <- merge(gtf,dea, by="geneid")
	dea <- dea[!duplicated(dea),]
	dea <- dea[order(dea$FDR),]
	write_tsv(dea,paste0(dir_path,"/DEG.tsv"), col_names=TRUE)
}
