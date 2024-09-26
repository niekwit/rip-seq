# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)
library(viridisLite)
library(rtracklayer)
library(openxlsx)

# Load Snakemake parameters
bed.files <- snakemake@input[["bed"]]
gtf <- snakemake@input[["gtf"]]

# Load GTF file
db <- rtracklayer::import(gtf)

# Extract relevant information
edb <- data.frame(geneId = db$gene_id, 
                  geneName = db$gene_name, 
                  geneBiotype = db$gene_biotype) %>%
  distinct()

# Load annotation database
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)

# Load sample information
sample_info <- read.csv("config/samples.csv", header = TRUE)

# Add sample names to bed files
samples <- sub(".*\\/([^\\/]+)\\.bed", "\\1", bed.files)
names(bed.files) <- samples

# Get sample condition from sample_info
conditions <- unique(str_replace(sample_info$sample, "_[0-9]+$", ""))

# Annotate bed files
peakAnnoList <- lapply(bed.files,
                       annotatePeak,
                       TxDb = txdb,
                       tssRegion = c(-3000, 3000)
                        )

# Plot binding relative to TSS
pdf(snakemake@output[["dt"]],
    width = 10,
    height = length(bed.files) * 2.5)
plotDistToTSS(peakAnnoList,
              title =  "Distribution of binding sites relative to TSS") +
  theme(axis.line.y = element_line(linewidth = 0),
        axis.line.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = viridis(6),
                    name = "Distance to TSS")
dev.off()

# Plot annotation bar
pdf(snakemake@output[["fd"]])
plotAnnoBar(peakAnnoList,
            title =  "Binding site distribution") +
  theme(axis.line.y = element_line(linewidth = 0),
        axis.line.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))
dev.off()

# Annotate peaks with gene information
annotation_list <- list()
for (i in seq_along(peakAnnoList)) {
  df <- as.data.frame(peakAnnoList[[i]])
  df <- left_join(df, edb, by = "geneId")
  # rename V4 column to peak_id
  colnames(df)[which(names(df) == "V4")] <- "peak_id"
  annotation_list[[i]] <- df
  names(annotation_list)[[i]] <- names(peakAnnoList[i])
}

write.xlsx(annotation_list,
           snakemake@output[["xlsx"]],
           colNames = TRUE)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")