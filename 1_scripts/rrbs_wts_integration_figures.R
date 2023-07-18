

# THIS SCRIPTS CONVERTS GENE IDS TO HUGO AND CHECK WHETER DMGs OVERLAPS WITH DEGs
# spoiler; they DO NOT
# only LUAD/LUSC genes overlap 

# Overlap of genes --------------------------------------------------------
# plasko_tc

overlap_plasko_tc <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values =  rownames(res_plasko_tc),
  mart = biomaRt::useDataset(
    "hsapiens_gene_ensembl",
    mart = biomaRt::useMart("ensembl", host = "uswest.ensembl.org")
  )
) %>%
  dplyr::rename('id' = 'ensembl_gene_id',
                'hgnc' = 'external_gene_name') %>%
  full_join(data.frame(res_plasko_tc) %>% 
              rownames_to_column(var = 'id')) %>% 
 # dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% 
  filter(!is.na(hgnc)) %>% 
  filter(hgnc != "") %>% 
  distinct(hgnc, .keep_all = T)

degs_plasko_tc$hgnc

myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol

intersect(degs_plasko_tc$hgnc, myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol)

phyper(
  length(intersect(degs_plasko_tc$hgnc, myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol))-1,
  length(myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol),
  length(unique(c(data.frame(myDiff_plasko__tc_all_annotated)$annot.symbol, degs_plasko_tc$hgnc))) - length(myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol),
  length(degs_plasko_tc$hgnc),
  lower.tail = FALSE,
  log.p = FALSE
)



set.seed(the_seed)
overlap_go_plasko_tc <- enrichGO(gene          = bitr(intersect(degs_plasko_tc$hgnc, myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol), 
                                                      fromType = "SYMBOL",
                                                      toType = c("ENTREZID"),
                                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                  universe      = bitr(unique(c(unique(overlap_plasko_tc$hgnc),data.frame(myDiff_plasko__tc_all_annotated)$annot.symbol)), 
                                                       fromType = "SYMBOL",
                                                       toType = c("ENTREZID"),
                                                       OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)

# plasko_tp

overlap_plasko_tp <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values =  rownames(res_plasko_tp),
  mart = biomaRt::useDataset(
    "hsapiens_gene_ensembl",
    mart = biomaRt::useMart("ensembl", host = "uswest.ensembl.org")
  )
) %>%
  dplyr::rename('id' = 'ensembl_gene_id',
                'hgnc' = 'external_gene_name') %>%
  full_join(data.frame(res_plasko_tp) %>% 
              rownames_to_column(var = 'id')) %>% 
  # dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% 
  filter(!is.na(hgnc)) %>% 
  filter(hgnc != "") %>% 
  distinct(hgnc, .keep_all = T)

degs_plasko_tp$hgnc

myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol

intersect(degs_plasko_tp$hgnc, myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol)

phyper(
  length(intersect(degs_plasko_tp$hgnc, myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol))-1,
  length(myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol),
  length(unique(c(data.frame(myDiff_plasko__tp_all_annotated)$annot.symbol, degs_plasko_tp$hgnc))) - length(myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol),
  length(degs_plasko_tp$hgnc),
  lower.tail = FALSE,
  log.p = FALSE
)



set.seed(the_seed)
overlap_go_plasko_tp <- enrichGO(gene          = bitr(intersect(degs_plasko_tp$hgnc, myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol), 
                                                      fromType = "SYMBOL",
                                                      toType = c("ENTREZID"),
                                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                    universe      = bitr(unique(c(unique(overlap_plasko_tp$hgnc),data.frame(myDiff_plasko__tp_all_annotated)$annot.symbol)), 
                                                         fromType = "SYMBOL",
                                                         toType = c("ENTREZID"),
                                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "ALL",
                                 pAdjustMethod = "BH",
                                 minGSSize     = 10,
                                 maxGSSize     = 500,
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)

# gruczo_tc

overlap_gruczo_tc <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values =  rownames(res_gruczo_tc),
  mart = biomaRt::useDataset(
    "hsapiens_gene_ensembl",
    mart = biomaRt::useMart("ensembl", host = "uswest.ensembl.org")
  )
) %>%
  dplyr::rename('id' = 'ensembl_gene_id',
                'hgnc' = 'external_gene_name') %>%
  full_join(data.frame(res_gruczo_tc) %>% 
              rownames_to_column(var = 'id')) %>% 
  # dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% 
  filter(!is.na(hgnc)) %>% 
  filter(hgnc != "") %>% 
  distinct(hgnc, .keep_all = T)

degs_gruczo_tc$hgnc

myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol

intersect(degs_gruczo_tc$hgnc, myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol)

phyper(
  length(intersect(degs_gruczo_tc$hgnc, myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol))-1,
  length(myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol),
  length(unique(c(data.frame(myDiff_gruczo__tc_all_annotated)$annot.symbol, degs_gruczo_tc$hgnc))) - length(myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol),
  length(degs_gruczo_tc$hgnc),
  lower.tail = FALSE,
  log.p = FALSE
)




set.seed(the_seed)
overlap_go_gruczo_tc <- enrichGO(gene          = bitr(intersect(degs_gruczo_tc$hgnc, myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol), 
                                                      fromType = "SYMBOL",
                                                      toType = c("ENTREZID"),
                                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                    universe      = bitr(unique(c(unique(overlap_gruczo_tc$hgnc),data.frame(myDiff_gruczo__tc_all_annotated)$annot.symbol)), 
                                                         fromType = "SYMBOL",
                                                         toType = c("ENTREZID"),
                                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "ALL",
                                 pAdjustMethod = "BH",
                                 minGSSize     = 10,
                                 maxGSSize     = 500,
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)


# gruczo_tp

overlap_gruczo_tp <- biomaRt::getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = 'ensembl_gene_id',
  values =  rownames(res_gruczo_tp),
  mart = biomaRt::useDataset(
    "hsapiens_gene_ensembl",
    mart = biomaRt::useMart("ensembl", host = "uswest.ensembl.org")
  )
) %>%
  dplyr::rename('id' = 'ensembl_gene_id',
                'hgnc' = 'external_gene_name') %>%
  full_join(data.frame(res_gruczo_tp) %>% 
              rownames_to_column(var = 'id')) %>% 
  # dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% 
  filter(!is.na(hgnc)) %>% 
  filter(hgnc != "") %>% 
  distinct(hgnc, .keep_all = T)

degs_gruczo_tp$hgnc

myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol

intersect(degs_gruczo_tp$hgnc, myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol)

phyper(
  length(intersect(degs_gruczo_tp$hgnc, myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol))-1,
  length(myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol),
  length(unique(c(data.frame(myDiff_gruczo__tp_all_annotated)$annot.symbol, degs_gruczo_tp$hgnc))) - length(myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol),
  length(degs_gruczo_tp$hgnc),
  lower.tail = FALSE,
  log.p = FALSE
)



set.seed(the_seed)
overlap_go_gruczo_tp <- enrichGO(gene          = bitr(intersect(degs_gruczo_tp$hgnc, myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol), 
                                                      fromType = "SYMBOL",
                                                      toType = c("ENTREZID"),
                                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                    universe      = bitr(unique(c(unique(overlap_gruczo_tp$hgnc),data.frame(myDiff_gruczo__tp_all_annotated)$annot.symbol)), 
                                                         fromType = "SYMBOL",
                                                         toType = c("ENTREZID"),
                                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db) %>% 
                                   dplyr::pull(ENTREZID),
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "ALL",
                                 pAdjustMethod = "BH",
                                 minGSSize     = 10,
                                 maxGSSize     = 500,
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)


# -------------------------------------------------------------------------

phyper(
  length(intersect(degs_plasko_tc$hgnc, degs_plasko_tp$hgnc))-1,
  length(unique(degs_plasko_tc$hgnc)),
  length(unique(c(overlap_plasko_tc$hgnc, overlap_plasko_tp$hgnc))) - length(unique(degs_plasko_tc$hgnc)),
  length(unique(degs_plasko_tp$hgnc)),
  lower.tail = FALSE,
  log.p = FALSE
)


phyper(
  length(intersect(degs_gruczo_tc$hgnc, degs_gruczo_tp$hgnc))-1,
  length(unique(degs_gruczo_tc$hgnc)),
  length(unique(c(overlap_gruczo_tc$hgnc, overlap_gruczo_tp$hgnc))) - length(unique(degs_gruczo_tc$hgnc)),
  length(unique(degs_gruczo_tp$hgnc)),
  lower.tail = FALSE,
  log.p = FALSE
)




