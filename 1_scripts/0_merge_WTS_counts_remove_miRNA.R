

# Dependencies ------------------------------------------------------------
library(purrr)
library(tidyverse)
library(biomaRt)

# Setup -------------------------------------------------------------------
# List 44 mice files
fc_files1 <-
  list.files(
    "0_data/featureCounts",
    pattern = ".txt",
    full.names = T
  )

# Concatenate vectors of files
#fc_files_final <- c(fc_files1, fc_files2)
fc_files_final <- fc_files1

# Function for reading the data in
read_in_feature_counts <- function(file) {
  cnt <- read_tsv(file, col_names = T, comment = "#")
  cnt <- cnt %>% dplyr::select(-Chr,-Start,-End,-Strand,-Length)
  return(cnt)
}

# Call the mapping parsing function and map records into list
raw_counts <- map(fc_files_final, read_in_feature_counts)

# Reduce to dataframe
raw_counts_df <- purrr::reduce(raw_counts, inner_join) 

# Remove reads mapping to miRNAs
gtf <- rtracklayer::import('0_data/Homo_sapiens.GRCh38.109.chr.gtf.gz')
gtf_df <- as.data.frame(gtf)
gtf_df_filtered <- gtf_df %>% 
  dplyr::select(gene_id, gene_biotype) %>% 
  filter(gene_biotype == 'miRNA')

# Filter out miRNAs 
`%nin%` = Negate(`%in%`) # negate %in%
raw_counts_df_no_miRNAs <- raw_counts_df %>% 
  filter(Geneid %nin% gtf_df_filtered$gene_id) %>% 
  arrange(Geneid)

# Add biological gene names form ENSEMBLIDs
# Allows to choose dataset and databse in one step
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
raw_counts_df_no_miRNAs_bio <-
  biomaRt::getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name'),
    filters = 'ensembl_gene_id',
    values =  raw_counts_df_no_miRNAs$Geneid,
    mart = ensembl) %>%
  dplyr::rename('Geneid' = 'ensembl_gene_id',
                'Gene' = 'external_gene_name') %>%
  full_join(raw_counts_df_no_miRNAs)

# Export data
# remove X from first colnames starting with numbers...
colnames(raw_counts_df_no_miRNAs_bio) <- gsub('^X', '', colnames(raw_counts_df_no_miRNAs_bio))
openxlsx::write.xlsx(raw_counts_df_no_miRNAs_bio,
                     '0_data/raw_counts_no-miRNAs.xlsx')













