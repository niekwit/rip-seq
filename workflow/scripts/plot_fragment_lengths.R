# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(RColorBrewer)

data <- read.delim(snakemake@input[[1]], 
                   header = TRUE, 
                   skip = 1,
                   stringsAsFactors = TRUE)

# Get max number of replicate samples
samples <- data$Sample
replicate_number <- str_split_fixed(samples, "_", 2)[,2] %>% 
  as.numeric() %>% 
  max()

# Add replicate to column (replicate_n)
data$replicate <- paste0("Replicate_", str_split_fixed(samples, "_", 2)[,2]) %>% 
  as.factor()

# Create vector for colours
unique_samples <- str_split_fixed(data$Sample, "_", 2)[,1] %>% 
  unique() %>%
  length()
colours_all <- brewer.pal(8,"Dark2")[1:unique_samples]
colours <- vector()
for (i in 1:replicate_number) {
  colours <- c(colours, rep(colours_all[i], replicate_number))
}

# Plot histogram with facet for each replicate
p <- ggplot(data %>% tidyr::uncount(Occurrences), 
            aes(x = Size,
                fill = Sample)) +
  geom_histogram(binwidth = 10,
                 position = "identity",
                 alpha = 0.65) +
  facet_wrap(~replicate, ncol = replicate_number) +
  scale_fill_manual(values = colours) +
  theme_cowplot(18)

ggsave(snakemake@output[[1]], 
       p, 
       width = 10, 
       height = 6)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")