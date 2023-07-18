
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(readxl)

# Setup -------------------------------------------------------------------
clinic_gruczo <- read_excel('0_data/clinic_revised_full_rrbs_wts.xlsx',
                            sheet = 'gruczolak_komplet') %>% 
  filter(WTS_Tn == 1 &
           WTS_Tp == 1 &
           WTS_Tc == 1 &
           rrbs_Tn == 1 &
           rrbs_Tp == 1 &
           rrbs_Tc == 1)

clinic_plasko <- read_excel('0_data/clinic_revised_full_rrbs_wts.xlsx',
                            sheet = 'plasko_komplet') %>% 
  filter(WTS_Tn == 1 &
           WTS_Tp == 1 &
           WTS_Tc == 1 &
           rrbs_Tn == 1 &
           rrbs_Tp == 1 &
           rrbs_Tc == 1)

# Read in sample names which had size < 100mb
mb100 <- read_excel('0_data/samples_100mb.xlsx')
mb100$samples <- str_replace(mb100$samples, '_R2_001.fastq.gz', '')
mb100$samples <- str_replace(mb100$samples, '_R1_001.fastq.gz', '')

# Read raw counts
raw_counts_wts <-
  read_excel(
    '0_data/raw_counts_no-miRNAs.xlsx'
  )

# Trim the col names
colnames(raw_counts_wts) <-
  stringr::str_replace(colnames(raw_counts_wts),
                       'Aligned.sortedByCoord.out.bam',
                       '')

# Filter out samples which size was less than 100mb
raw_counts_wts_100mb <- raw_counts_wts %>%
  dplyr::select(!contains(mb100$samples))


# Plasko
raw_counts_wts_100mb_plasko_tc <-
  raw_counts_wts_100mb[str_extract(colnames(raw_counts_wts_100mb), "[^_]+") %in% c('Geneid', clinic_plasko$pacjent)] %>% 
  dplyr::select(Geneid, contains('Tc') |
                  contains('Tn')) %>%
  column_to_rownames(var = 'Geneid')

raw_counts_wts_100mb_plasko_tp <-
  raw_counts_wts_100mb[str_extract(colnames(raw_counts_wts_100mb), "[^_]+") %in% c('Geneid', clinic_plasko$pacjent)] %>% 
  dplyr::select(Geneid, contains('Tp') |
                  contains('Tn')) %>%
  column_to_rownames(var = 'Geneid')

# Gruczo
raw_counts_wts_100mb_gruczo_tc <-
  raw_counts_wts_100mb[str_extract(colnames(raw_counts_wts_100mb), "[^_]+") %in% c('Geneid', clinic_gruczo$pacjent)] %>% 
  dplyr::select(Geneid, contains('Tc') |
                  contains('Tn')) %>%
  column_to_rownames(var = 'Geneid')

raw_counts_wts_100mb_gruczo_tp <-
  raw_counts_wts_100mb[str_extract(colnames(raw_counts_wts_100mb), "[^_]+") %in% c('Geneid', clinic_gruczo$pacjent)] %>% 
  dplyr::select(Geneid, contains('Tp') |
                  contains('Tn')) %>%
  column_to_rownames(var = 'Geneid')

# meta
meta_plasko_tc <- data.frame(samples = colnames(raw_counts_wts_100mb_plasko_tc),
                             # Extract match between first two underscores
                             condition = factor(
                               gsub(
                                 "(?:[^_]+_){1}([^_]+).*",
                                 "\\1",
                                 colnames(raw_counts_wts_100mb_plasko_tc)
                               ),
                               levels = c('Tn', 'Tc')
                             ))

meta_plasko_tp <- data.frame(samples = colnames(raw_counts_wts_100mb_plasko_tp),
                             # Extract match between first two underscores
                             condition = factor(
                               gsub(
                                 "(?:[^_]+_){1}([^_]+).*",
                                 "\\1",
                                 colnames(raw_counts_wts_100mb_plasko_tp)
                               ),
                               levels = c('Tn', 'Tp')
                             ))

meta_gruczo_tc <- data.frame(samples = colnames(raw_counts_wts_100mb_gruczo_tc),
                             # Extract match between first two underscores
                             condition = factor(
                               gsub(
                                 "(?:[^_]+_){1}([^_]+).*",
                                 "\\1",
                                 colnames(raw_counts_wts_100mb_gruczo_tc)
                               ),
                               levels = c('Tn', 'Tc')
                             ))

meta_gruczo_tp <- data.frame(samples = colnames(raw_counts_wts_100mb_gruczo_tp),
                             # Extract match between first two underscores
                             condition = factor(
                               gsub(
                                 "(?:[^_]+_){1}([^_]+).*",
                                 "\\1",
                                 colnames(raw_counts_wts_100mb_gruczo_tp)
                               ),
                               levels = c('Tn', 'Tp')
                             ))

# DE Analysis -------------------------------------------------------------
library(DESeq2)
# Create DESeq2 object
dds_plasko_tc <- DESeqDataSetFromMatrix(countData = raw_counts_wts_100mb_plasko_tc,
                                        colData = meta_plasko_tc,
                                        design = ~ condition)
if (all(colnames(raw_counts_wts_100mb_plasko_tc) == meta_plasko_tc$samples)) {
  print("OK")
} else {
  stop("samples not in proper order for DESeq2")
}
keep_plasko_tc <- rowSums(counts(dds_plasko_tc)) >= 20
dds_f_plasko_tc <- dds_plasko_tc[keep_plasko_tc, ]
dds_f_deseq_plasko_tc <- DESeq(dds_f_plasko_tc)
res_plasko_tc <- results(dds_f_deseq_plasko_tc)
summary(res_plasko_tc)

dds_plasko_tp <- DESeqDataSetFromMatrix(countData = raw_counts_wts_100mb_plasko_tp,
                                        colData = meta_plasko_tp,
                                        design = ~ condition)
if (all(colnames(raw_counts_wts_100mb_plasko_tp) == meta_plasko_tp$samples)) {
  print("OK")
} else {
  stop("samples not in proper order for DESeq2")
}
keep_plasko_tp <- rowSums(counts(dds_plasko_tp)) >= 20
dds_f_plasko_tp <- dds_plasko_tp[keep_plasko_tp, ]
dds_f_deseq_plasko_tp <- DESeq(dds_f_plasko_tp)
res_plasko_tp <- results(dds_f_deseq_plasko_tp)
summary(res_plasko_tp)

dds_gruczo_tc <- DESeqDataSetFromMatrix(countData = raw_counts_wts_100mb_gruczo_tc,
                                        colData = meta_gruczo_tc,
                                        design = ~ condition)
if (all(colnames(raw_counts_wts_100mb_gruczo_tc) == meta_gruczo_tc$samples)) {
  print("OK")
} else {
  stop("samples not in proper order for DESeq2")
}
keep_gruczo_tc <- rowSums(counts(dds_gruczo_tc)) >= 20
dds_f_gruczo_tc <- dds_gruczo_tc[keep_gruczo_tc, ]
dds_f_deseq_gruczo_tc <- DESeq(dds_f_gruczo_tc)
res_gruczo_tc <- results(dds_f_deseq_gruczo_tc)
summary(res_gruczo_tc)

dds_gruczo_tp <- DESeqDataSetFromMatrix(countData = raw_counts_wts_100mb_gruczo_tp,
                                        colData = meta_gruczo_tp,
                                        design = ~ condition)
if (all(colnames(raw_counts_wts_100mb_gruczo_tp) == meta_gruczo_tp$samples)) {
  print("OK")
} else {
  stop("samples not in proper order for DESeq2")
}
keep_gruczo_tp <- rowSums(counts(dds_gruczo_tp)) >= 20
dds_f_gruczo_tp <- dds_gruczo_tp[keep_gruczo_tp, ]
dds_f_deseq_gruczo_tp <- DESeq(dds_f_gruczo_tp)
res_gruczo_tp <- results(dds_f_deseq_gruczo_tp)
summary(res_gruczo_tp)


# Get DEGs ----------------------------------------------------------------

degs_plasko_tc <- res_plasko_tc %>% 
  data.frame() %>% 
  rownames_to_column(var = 'Geneid') %>% 
  dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05)
degs_plasko_tc$metadata <- rep(res_plasko_tc@elementMetadata$description, length.out = nrow(degs_plasko_tc))

degs_plasko_tp <- res_plasko_tp %>% 
  data.frame() %>% 
  rownames_to_column(var = 'Geneid') %>% 
  dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05)
degs_plasko_tp$metadata <- rep(res_plasko_tp@elementMetadata$description, length.out = nrow(degs_plasko_tp))

degs_gruczo_tc <- res_gruczo_tc %>% 
  data.frame() %>% 
  rownames_to_column(var = 'Geneid') %>% 
  dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05)
degs_gruczo_tc$metadata <- rep(res_gruczo_tc@elementMetadata$description, length.out = nrow(degs_gruczo_tc))

degs_gruczo_tp <- res_gruczo_tp %>% 
  data.frame() %>% 
  rownames_to_column(var = 'Geneid') %>% 
  dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05)
degs_gruczo_tp$metadata <- rep(res_gruczo_tp@elementMetadata$description, length.out = nrow(degs_gruczo_tp))

# Heatmap -----------------------------------------------------------------
all_samples <- c(
  colnames(raw_counts_wts_100mb_plasko_tc),
  colnames(raw_counts_wts_100mb_plasko_tp),
  colnames(raw_counts_wts_100mb_gruczo_tc),
  colnames(raw_counts_wts_100mb_gruczo_tp)
)

all_samples_degs <- c(
  degs_plasko_tc$Geneid,
  degs_plasko_tp$Geneid,
  degs_gruczo_tc$Geneid,
  degs_gruczo_tp$Geneid
)

# Get DEGs from count matrix
raw_counts_wts_100mb2 <- raw_counts_wts_100mb %>% 
  filter(Geneid %in% all_samples_degs) %>% 
  column_to_rownames(var = 'Geneid')

# Make sure only samples from DE are selected
all_samples_counts <- raw_counts_wts_100mb2[colnames(raw_counts_wts_100mb2) %in% all_samples]

meta_heatmap <- data.frame(
  sample = colnames(all_samples_counts),
  condition = 1
)

dds_heatmap <- DESeqDataSetFromMatrix(countData = all_samples_counts, colData = meta_heatmap, design =~1)
dds_heatmap2 <- estimateSizeFactors(dds_heatmap)
normalized_counts_heatmap <- counts(dds_heatmap2, normalized=TRUE)


data_TPM_HGNC_heatmap <- normalized_counts_heatmap
# To properly calculate the z-score, the data has to be will be approximately 
# normally distributed and suitable for calculating z-scores. Z-score are 
# used to enhance the visualisation of the data (it helps to even the values out).
TPM_heatmap_in_log2 <- data_TPM_HGNC_heatmap %>% 
  as.matrix() %>% 
  log2()

TPM_heatmap_in_log2[is.infinite(TPM_heatmap_in_log2)] <- 0

column_names <- colnames(TPM_heatmap_in_log2)

TPM_heatmap_in_log2_scaled <- base::apply(TPM_heatmap_in_log2, 1, scale) %>% 
  t()

colnames(TPM_heatmap_in_log2_scaled) <- column_names

TPM_heatmap_in_log2_scaled[is.nan(TPM_heatmap_in_log2_scaled)] <- 0

col_fun2 = circlize::colorRamp2(c(-2, 0, 2),
                                c('mediumblue', 'white', 'red2'), space = 'sRGB')



colnames(TPM_heatmap_in_log2_scaled) <- paste0(str_replace(str_replace(str_replace(sapply(str_split(colnames(TPM_heatmap_in_log2_scaled), '_'), "[[" , 2), 'Tn', '1_Tn'), 'Tp', '2_Tp'), 'Tc', '3_Tc'),'_', colnames(TPM_heatmap_in_log2_scaled))
TPM_heatmap_in_log2_scaled <- TPM_heatmap_in_log2_scaled[, order(colnames(TPM_heatmap_in_log2_scaled))]

TPM_heatmap_in_log2_scaled_gruczo <- TPM_heatmap_in_log2_scaled[,sapply(strsplit(colnames(TPM_heatmap_in_log2_scaled), '_'), '[', 3) %in% clinic_gruczo$pacjent]
TPM_heatmap_in_log2_scaled_plasko <- TPM_heatmap_in_log2_scaled[,sapply(strsplit(colnames(TPM_heatmap_in_log2_scaled), '_'), '[', 3) %in% clinic_plasko$pacjent]


library(ComplexHeatmap)
ha_Gruczo = HeatmapAnnotation(bar = sapply(strsplit(colnames(TPM_heatmap_in_log2_scaled_gruczo),"_"), '[', 2),
                       col = list(bar = c("Tn" = "#ebc999", 
                                          "Tp" = "#cd7700", 
                                          "Tc" = "#4d3227")))

ha_Plasko = HeatmapAnnotation(bar = sapply(strsplit(colnames(TPM_heatmap_in_log2_scaled_plasko),"_"), '[', 2),
                              col = list(bar = c("Tn" = "#37a987", 
                                                 "Tp" = "#b7b1d2", 
                                                 "Tc" = "#4b3d8f")))


temp2_gruczo <- ComplexHeatmap::Heatmap(TPM_heatmap_in_log2_scaled_gruczo, name = 'z-score', col = col_fun2,
                                 cluster_rows = T, cluster_columns = FALSE, show_row_names = F,
                                 border = F, show_column_names = F,column_title = '',
                                 row_dend_width = unit(1, 'cm'), heatmap_legend_param = list(legend_height = unit(6, "cm")),
                                 top_annotation = ha_Gruczo
)

temp2_plasko <- ComplexHeatmap::Heatmap(TPM_heatmap_in_log2_scaled_plasko, name = 'z-score', col = col_fun2,
                                        cluster_rows = T, cluster_columns = FALSE, show_row_names = F,
                                        border = F, show_column_names = F,column_title = '',
                                        row_dend_width = unit(1, 'cm'), heatmap_legend_param = list(legend_height = unit(6, "cm")),
                                        top_annotation = ha_Plasko
)

szer_gru <- ncol(TPM_heatmap_in_log2_scaled_gruczo)/12

plas_gru <- ncol(TPM_heatmap_in_log2_scaled_plasko)/12

tiff("2_results/DEGs_heatmaps/2top_ann_gruczo.tiff", width = szer_gru, height = 6, units = 'in', res = 300)
ComplexHeatmap::draw(temp2_gruczo, heatmap_legend_side="left")
dev.off()

tiff("2_results/DEGs_heatmaps/2top_ann_plasko.tiff", width = plas_gru, height = 6, units = 'in', res = 300)
ComplexHeatmap::draw(temp2_plasko, heatmap_legend_side="left")
dev.off()

# Venn --------------------------------------------------------------------
venn_list <- list(
  p_tc = degs_plasko_tc$Geneid,
  p_tp = degs_plasko_tp$Geneid,
  g_tc = degs_gruczo_tc$Geneid,
  g_tp = degs_gruczo_tp$Geneid
)
library(VennDiagram)
pdf("2_results/VennDiagram/VennDiagram.pdf", width = 5, height = 5)
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(
  venn_list,
  category.names = c("p_tc" , "p_tp" , "g_tc", "g_tp"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#4b3d8f",
           "#b7b1d2",
           "#4d3227",
           "#cd7700"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1)
)

dev.off()


# GSEA  -------------------------------------------------------------
library(clusterProfiler)
library(DOSE)
library(ReactomePA)

# Data processing for GSEA Algorithm to GO and KEGG
res_plasko_tc_entrez <- bitr(rownames(res_plasko_tc), 
                             fromType = "ENSEMBL",
                             toType = c("ENTREZID"),
                             OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
  rename(entrezgene_id = ENTREZID,
         id = ENSEMBL) %>% 
  full_join(data.frame(res_plasko_tc) %>% rownames_to_column(var = 'id'))

geneList_res_plasko_tc_entrez <- res_plasko_tc_entrez$log2FoldChange
names(geneList_res_plasko_tc_entrez) <- res_plasko_tc_entrez$entrezgene_id
geneList_res_plasko_tc_entrez <- sort(geneList_res_plasko_tc_entrez, decreasing = TRUE)

res_plasko_tp_entrez <- bitr(rownames(res_plasko_tp), 
                             fromType = "ENSEMBL",
                             toType = c("ENTREZID"),
                             OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
  rename(entrezgene_id = ENTREZID,
         id = ENSEMBL) %>% 
  full_join(data.frame(res_plasko_tp) %>% rownames_to_column(var = 'id'))

geneList_res_plasko_tp_entrez <- res_plasko_tp_entrez$log2FoldChange
names(geneList_res_plasko_tp_entrez) <- res_plasko_tp_entrez$entrezgene_id
geneList_res_plasko_tp_entrez <- sort(geneList_res_plasko_tp_entrez, decreasing = TRUE)

res_gruczo_tc_entrez <- bitr(rownames(res_gruczo_tc), 
                             fromType = "ENSEMBL",
                             toType = c("ENTREZID"),
                             OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
  rename(entrezgene_id = ENTREZID,
         id = ENSEMBL) %>% 
  full_join(data.frame(res_gruczo_tc) %>% rownames_to_column(var = 'id'))

geneList_res_gruczo_tc_entrez <- res_gruczo_tc_entrez$log2FoldChange
names(geneList_res_gruczo_tc_entrez) <- res_gruczo_tc_entrez$entrezgene_id
geneList_res_gruczo_tc_entrez <- sort(geneList_res_gruczo_tc_entrez, decreasing = TRUE)

res_gruczo_tp_entrez <- bitr(rownames(res_gruczo_tp), 
                             fromType = "ENSEMBL",
                             toType = c("ENTREZID"),
                             OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
  rename(entrezgene_id = ENTREZID,
         id = ENSEMBL) %>% 
  full_join(data.frame(res_gruczo_tp) %>% rownames_to_column(var = 'id'))

geneList_res_gruczo_tp_entrez <- res_gruczo_tp_entrez$log2FoldChange
names(geneList_res_gruczo_tp_entrez) <- res_gruczo_tp_entrez$entrezgene_id
geneList_res_gruczo_tp_entrez <- sort(geneList_res_gruczo_tp_entrez, decreasing = TRUE)

# create element as an Sys.time() object - 1648568306
(the_seed = 1648568306 %% 100000) # 50700

set.seed(the_seed)
plasko_tc_go_ora_wts <- gseGO(geneList          = geneList_res_plasko_tc_entrez,
                              OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1)
set.seed(the_seed)
plasko_tp_go_ora_wts <- gseGO(geneList          = geneList_res_plasko_tp_entrez,
                              OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1)
set.seed(the_seed)
gruczo_tc_go_ora_wts <- gseGO(geneList          = geneList_res_gruczo_tc_entrez,
                              OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1)
set.seed(the_seed)
gruczo_tp_go_ora_wts <- gseGO(geneList          = geneList_res_gruczo_tp_entrez,
                              OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1)

# KEGG over-rep -----------------------------------------------------------
set.seed(the_seed)
plasko_tc_kegg_ora_wts <- gseKEGG(geneList = geneList_res_plasko_tc_entrez,
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  organism     = 'hsa',
                                  pvalueCutoff = 1)
set.seed(the_seed)
plasko_tp_kegg_ora_wts <- gseKEGG(geneList = geneList_res_plasko_tp_entrez,
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  organism     = 'hsa',
                                  pvalueCutoff = 1)
set.seed(the_seed)
gruczo_tc_kegg_ora_wts <- gseKEGG(geneList = geneList_res_gruczo_tc_entrez,
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  organism     = 'hsa',
                                  pvalueCutoff = 1
)
set.seed(the_seed)
gruczo_tp_kegg_ora_wts <- gseKEGG(geneList = geneList_res_gruczo_tp_entrez,
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  organism     = 'hsa',
                                  pvalueCutoff = 1
)

# Wiki-Pathways Ora -------------------------------------------------------
set.seed(the_seed)
plasko_tc_wp_ora_wts <- gseWP(geneList =geneList_res_plasko_tc_entrez, 
                              organism      = "Homo sapiens",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1
) 
set.seed(the_seed)
plasko_tp_wp_ora_wts <- gseWP(geneList =geneList_res_plasko_tp_entrez, 
                              organism      = "Homo sapiens",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1
) 
set.seed(the_seed)
gruczo_tc_wp_ora_wts <- gseWP(geneList =geneList_res_gruczo_tc_entrez, 
                              organism      = "Homo sapiens",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1
) 
set.seed(the_seed)
gruczo_tp_wp_ora_wts <- gseWP(geneList =geneList_res_gruczo_tp_entrez, 
                              organism      = "Homo sapiens",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              pvalueCutoff  = 1
) 


# Reactome ORA ------------------------------------------------------------

set.seed(the_seed)
plasko_tc_reactome_ora_wts <- gsePathway(geneList =geneList_res_plasko_tc_entrez, 
                                         organism = "human",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = "BH",
                                         minGSSize = 10,
                                         maxGSSize = 500)

set.seed(the_seed)
plasko_tp_reactome_ora_wts <- gsePathway(geneList =geneList_res_plasko_tp_entrez, 
                                         organism = "human",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = "BH",
                                         minGSSize = 10,
                                         maxGSSize = 500)

set.seed(the_seed)
gruczo_tc_reactome_ora_wts <- gsePathway(geneList =geneList_res_gruczo_tc_entrez, 
                                         organism = "human",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = "BH",
                                         minGSSize = 10,
                                         maxGSSize = 500)

set.seed(the_seed)
gruczo_tp_reactome_ora_wts <- gsePathway(geneList =geneList_res_gruczo_tp_entrez, 
                                         organism = "human",
                                         pvalueCutoff = 1,
                                         pAdjustMethod = "BH",
                                         minGSSize = 10,
                                         maxGSSize = 500)


# Disease Ontology --------------------------------------------------------
set.seed(the_seed)
plasko_tc_do_ora_wts <- gseDO(geneList =geneList_res_plasko_tc_entrez, 
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              verbose      = TRUE)

set.seed(the_seed)
plasko_tp_do_ora_wts <- gseDO(geneList =geneList_res_plasko_tp_entrez, 
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              verbose      = TRUE)

set.seed(the_seed)
gruczo_tc_do_ora_wts <- gseDO(geneList =geneList_res_gruczo_tc_entrez, 
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              verbose      = TRUE)

set.seed(the_seed)
gruczo_tp_do_ora_wts <- gseDO(geneList =geneList_res_gruczo_tp_entrez, 
                              pvalueCutoff  = 1,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500,
                              verbose      = TRUE)


# Network of Cancer Gene --------------------------------------------------
set.seed(the_seed)
plasko_tc_ncg_ora2_wts <- gseNCG(geneList =geneList_res_plasko_tc_entrez, 
                                 pvalueCutoff  = 1,
                                 pAdjustMethod = "BH",
                                 minGSSize     = 10,
                                 maxGSSize     = 500)

set.seed(the_seed)
plasko_tp_ncg_ora_wts <- gseNCG(geneList =geneList_res_plasko_tp_entrez, 
                                
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

set.seed(the_seed)
gruczo_tc_ncg_ora_wts <- gseNCG(geneList =geneList_res_gruczo_tc_entrez, 
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

set.seed(the_seed)
gruczo_tp_ncg_ora_wts <- gseNCG(geneList =geneList_res_gruczo_tp_entrez, 
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

# Network of Disease Gene Network -----------------------------------------
set.seed(the_seed)
plasko_tc_dgn_ora_wts <- gseDGN(geneList =geneList_res_plasko_tc_entrez, 
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

set.seed(the_seed)
plasko_tp_dgn_ora_wts <- gseDGN(geneList =geneList_res_plasko_tp_entrez, 
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

set.seed(the_seed)
gruczo_tc_dgn_ora_wts <- gseDGN(geneList =geneList_res_gruczo_tc_entrez, 
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

set.seed(the_seed)
gruczo_tp_dgn_ora_wts <- gseDGN(geneList =geneList_res_gruczo_tp_entrez, 
                                pvalueCutoff  = 1,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)

# Export ------------------------------------------------------------------

# Exported at 13 Jul 2023 11:00 AM
save.image(file='2_results/Rdata_objects/1_WTS.RData')

