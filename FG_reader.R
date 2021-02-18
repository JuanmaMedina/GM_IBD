library("tidyverse")

# Set path for FG (i.e. Fermentative Genes) database
setwd("~/ANALYSIS/IBD_cohort/results_new/RPOB/")

# Read metadata and remove samples with less than 1kk reads (these output abnormally high values in the normalization)
samples <- read.table("~/ANALYSIS/IBD_cohort/metadata_prep/metadata_MGX_reads.txt")[c(1,2,4,5)]
colnames(samples) <- c("sample","patient","status","kkreads")
samples <- subset(samples, samples[ , 4] > 1)

# N of each group
table(samples$status)


read_plus <- function(f) {
  # This function reads a .tsv creating a new variable with the sample name
  read_tsv(f, col_names = F) %>%
    mutate(filename = f)
}


# Read files from folder (https://community.rstudio.com/t/64629)
d <- list.files(path = ".", pattern = "matches") %>%
  
  # Rbind them as DFs
  map_dfr(~read_plus(.)) %>%                                                      
  setNames(., c("genus","strain","gene_ID","evalue","pid","sample"))


# Remove "_matches.m8" extension
d$sample <- substr(d$sample,1,nchar(d$sample)-11)

d %>%
  select("gene_ID","sample") %>%

  # Number of reads aligning to the gene in each sample
  add_count(sample, name = "reads") %>%

  # Keep one record for sample
  distinct(sample, .keep_all = T) %>%
  {. ->> d_all}


d %>%
  select("genus","strain","gene_ID","sample") %>%
  
  # Number of reads aligning to the gene in each bacterial genus
  add_count(genus,sample, name = "reads") %>%
  
  # Keep one record for genus and sample
  distinct(genus,sample, .keep_all = T) %>%                                       
  {. ->> d_genus}


# Merge with metadata
d_all <- merge(d_all,samples,by = "sample")
d_genus <- merge(d_genus,samples,by = "sample")

# Normalize gene counts by library size and remove column
d_all$reads <- round(d_all$reads / d_all$kkreads)
d_all$kkreads <- NULL

d_genus$reads <- round(d_genus$reads / d_genus$kkreads)
d_genus$kkreads <- NULL

# Writing results
write.table(d_all, file="~/FGreader_out/rpob_all.tsv", quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(d_genus, file="~/FGreader_out/rpob_genus.tsv", quote=FALSE, sep='\t', col.names = F, row.names = F)
