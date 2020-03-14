smoc2_rawcounts <- fibrosis #GENE-ID row

# Create genotype vector
genotype <- c('smoc2_oe', 'smoc2_oe', 'smoc2_oe', 'smoc2_oe', 'smoc2_oe' ,'smoc2_oe', 'smoc2_oe')

# Create condition vector
condition <- c('fibrosis', 'fibrosis','fibrosis', 'fibrosis', 'normal', 'normal', 'normal')

# Create data frame
smoc2_metadata <- data.frame(genotype, condition)

# Assign the row names of the data frame
rownames(smoc2_metadata) <- c('smoc2_fibrosis1', 'smoc2_fibrosis2', 
                              'smoc2_fibrosis3', 'smoc2_fibrosis4', 'smoc2_normal1',
                              'smoc2_normal4', 'smoc2_normal3')


'''
DESeq2 and edgeR are the most popular model for differential expression anaysis
both use negative bionomial model. DESeq2 is part of bioconductor package. 

'''
# Use the match() function to reorder the columns of the raw counts
# Use the match() function to reorder the columns of the raw counts
match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))

# Reorder the columns of the count data
reordered_smoc2_rawcounts <- smoc2_rawcounts[ , match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))]

#first we have to create DESeq2 object

library(DESeq2)
# Create a DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                                    colData = smoc2_metadata,
                                    design = ~ condition)


#determine the size factors to use normalization
dds_smoc2 <- estimateSizeFactors(dds_smoc2)

#extract the normalized counts
smoc2_normalized_counts <- counts(dds_smoc2, normalize = TRUE)


# Transform the normalized counts 
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)

# Compute the correlation values between samples
vsd_cor_smoc2 <- cor(vsd_mat_smoc2) 

library(dplyr)
library(pheatmap)
# Plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(smoc2_metadata, condition))

# Plot dispersions
plotDispEsts(dds_smoc2)


# Create DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                                    colData = smoc2_metadata,
                                    design = ~ condition)

# Run the DESeq2 analysis
dds_smoc2 <- DESeq(dds_smoc2)

# Explore the results() function
?results

# Extract results
smoc2_res <- results(dds_smoc2, 
                     contrast = c("condition", "fibrosis", "normal"), 
                     alpha = 0.05, 
                     lfcThreshold = 0.32)

# Shrink the log2 fold changes
smoc2_res <- lfcShrink(dds_smoc2, 
                       contrast = c("condition", "fibrosis", "normal"), 
                       res = smoc2_res)

# Get an overview of the results    
summary(smoc2_res)

# Save results as a data frame
smoc2_res_all <- data.frame(smoc2_res)

load("C:/Users/laayt/Downloads/grcm38.rda")

new_data1 <- cbind(new_data$X, as.data.frame(smoc2_normalized_counts))
smoc2_res_all <- cbind(new_data1, smoc2_res_all$padj)



# Subset the results to only return the significant genes with p-adjusted values less than 0.05
smoc2_res_sig <- subset(smoc2_res_all, `smoc2_res_all$padj` < 0.05)

top_20 <- data.frame(smoc2_res_sig)[1:20, ] %>%
  rownames_to_column(var = "ensgene")

top_20 <- inner_join(top_20,
                     rownames_to_column(smoc2_metadata, var = "samplename"),
                     by = "samplename")

top_20 <- gather(top_20, 
                 key = "samplename", 
                 value = "normalized_counts", 
                 2:8)

library(ggplot2)
ggplot(top_20) +
  geom_point(aes(x = ensgene, y = normalized_counts, color = condition)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))


library(RCurl)
## paste URL to make it easier to read code (cosmetic!)
dat_url <- paste0("https://github.com/stephenturner/annotables/blob/master/data/grcm38.rda")
f <- getBinaryURL()
L <- load(rawConnection(f))
