
# Dependencies ------------------------------------------------------------
library(methylKit)
library(magrittr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggsignif)
library(ggpubr)

# Read methylation dbs from methylkit -------------------------------------
# normalization by coverage was also applied (check previous script for more)
plasko_tc <- readMethylDB('/Users/andrzejeljaszewicz/Adrian/Omics/MOBIT/mobit_rrbs/methylDB_plasko_Tc_vs_Tn/methylBase_4cc0236cb87.txt.bgz')
plasko_tp <- readMethylDB('/Users/andrzejeljaszewicz/Adrian/Omics/MOBIT/mobit_rrbs/methylDB_plasko_Tp_vs_Tn/04042022_tp_vs_tn/methylBase_5642f005095.txt.bgz')

gruczo_tc <- readMethylDB('/Users/andrzejeljaszewicz/Adrian/Omics/MOBIT/mobit_rrbs/methylDB_Tc_Tn/methylBase_14bbb14c77394.txt.bgz')
gruczo_tp <- readMethylDB('/Users/andrzejeljaszewicz/Adrian/Omics/MOBIT/mobit_rrbs/methylDB_Tp_Tn/methylBase_15127190b7fe8.txt.bgz')

# Kepp only 1:22 and X chromosomes ----------------------------------------
plasko_tc1 <- getData(plasko_tc) 
plasko_tp1 <- getData(plasko_tp) 
v <- c(1:22,'X')
plasko_tc2 <- plasko_tc[which(plasko_tc1$chr %in% v),]
plasko_tp2 <- plasko_tp[which(plasko_tp1$chr %in% v),]
#plasko_tp2_temp <- plasko_tp[which(plasko_tp1$chr %in% v),]
# plasko_tp2_temp@sample.ids[14] <- 'KT19266e_Tp.1'
# plasko_tp2 <- reorganize(plasko_tp2_temp,
#            sample.ids = plasko_tp2_temp@sample.ids[-14],
#            treatment = plasko_tp2_temp@treatment[-14])
plasko_tc3 <- getData(plasko_tc2)
plasko_tp3 <- getData(plasko_tp2)

gruczo_tc1 <- getData(gruczo_tc) 
gruczo_tp1 <- getData(gruczo_tp) 
v <- c(1:22,'X')
gruczo_tc2 <- gruczo_tc[which(gruczo_tc1$chr %in% v),]
gruczo_tp2 <- gruczo_tp[which(gruczo_tp1$chr %in% v),]
gruczo_tc3 <- getData(gruczo_tc2)
gruczo_tp3 <- getData(gruczo_tp2)

# Diff meth ---------------------------------------------------------------
myDiff_plasko__tc <- calculateDiffMeth(plasko_tc2, mc.cores = 6)
myDiff_plasko__tp <- calculateDiffMeth(plasko_tp2, mc.cores = 6)

myDiff_gruczo__tc <- calculateDiffMeth(gruczo_tc2, mc.cores = 6)
myDiff_gruczo__tp <- calculateDiffMeth(gruczo_tp2, mc.cores = 6)

# Percentage matrix -------------------------------------------------------
plasko_tc_pm <- percMethylation(plasko_tc2)
plasko_tp_pm <- percMethylation(plasko_tp2)

gruczo_tc_pm <- percMethylation(gruczo_tc2)
gruczo_tp_pm <- percMethylation(gruczo_tp2)

# Calculate global methylation --------------------------------------------
plasko_tn_pm_sub <- data.frame(plasko_tc_pm) %>% 
  select(grep('_Tn', colnames(plasko_tc_pm)))
plasko_tp_pm_sub <- data.frame(plasko_tp_pm) %>% 
  select(grep('_Tp', colnames(plasko_tp_pm)))
plasko_tc_pm_sub <- data.frame(plasko_tc_pm) %>% 
  select(grep('_Tc', colnames(plasko_tc_pm)))

gruczo_tn_pm_sub <- data.frame(gruczo_tc_pm) %>% 
  select(grep('_Tn', colnames(gruczo_tc_pm)))
gruczo_tp_pm_sub <- data.frame(gruczo_tp_pm) %>% 
  select(grep('_Tp', colnames(gruczo_tp_pm)))
gruczo_tc_pm_sub <- data.frame(gruczo_tc_pm) %>% 
  select(grep('_Tc', colnames(gruczo_tc_pm)))

# Prepare perc matrix for plotting ----------------------------------------
global_pm <- data.frame(
  mean_perc_meth = 
    c(Tn_plasko = apply(plasko_tn_pm_sub, 2, mean),
      Tp_plasko = apply(plasko_tp_pm_sub, 2, mean),
      Tc_plasko = apply(plasko_tc_pm_sub, 2, mean),
      Tn_gruczo = apply(gruczo_tn_pm_sub, 2, mean),
      Tp_gruczo = apply(gruczo_tp_pm_sub, 2, mean),
      Tc_gruczo = apply(gruczo_tc_pm_sub, 2, mean)),
  sample_name = names(c(Tn_plasko = apply(plasko_tn_pm_sub, 2, mean),
                        Tp_plasko = apply(plasko_tp_pm_sub, 2, mean),
                        Tc_plasko = apply(plasko_tc_pm_sub, 2, mean),
                        Tn_gruczo = apply(gruczo_tn_pm_sub, 2, mean),
                        Tp_gruczo = apply(gruczo_tp_pm_sub, 2, mean),
                        Tc_gruczo = apply(gruczo_tc_pm_sub, 2, mean)))
)

global_pm2 <- global_pm %>% 
  mutate(group = sub('\\..*', '', global_pm$sample_name)) %>%
  mutate(pheno = sub('*.._', '', group))
#filter(group %in% c('Tn_plasko', 'Tp_plasko', 'Tc_plasko'))

global_pm2$group <- factor(global_pm2$group,
                           levels = c('Tn_gruczo', 'Tp_gruczo', 'Tc_gruczo',
                                      'Tn_plasko', 'Tp_plasko', 'Tc_plasko'))

global_pm2$pheno <- factor(global_pm2$pheno,
                           levels = c('gruczo', 'plasko'))

global_pm2_kruskal_gruczo <- kruskal_test(global_pm2[global_pm2$pheno=='gruczo',],
                                          mean_perc_meth ~ group)
global_pm2_kruskal_plasko <- kruskal_test(global_pm2[global_pm2$pheno=='plasko',],
                                          mean_perc_meth ~ group)

global_pm2_dunn_gruczo <- dunn_test(data = global_pm2[global_pm2$pheno=='gruczo',], mean_perc_meth ~ group, p.adjust.method = "holm")
global_pm2_dunn_plasko <- dunn_test(data = global_pm2[global_pm2$pheno=='plasko',], mean_perc_meth ~ group, p.adjust.method = "holm")

global_pm2_dunn_gruczo2 <- add_xy_position(global_pm2_dunn_gruczo, x = "group", step.increase = 0.075)
global_pm2_dunn_gruczo2$y.position <- c(54.8, 58.7, 62.8)
pdf('gruczo.pdf', 3,5)
(p1 <- ggplot(global_pm2[global_pm2$pheno=='gruczo',], aes(x=group, y=mean_perc_meth, color = group)) + 
    #geom_violin(trim = TRUE) +
    stat_summary(fun = "median", fun.min = "median", fun.max= "median", size= 0.3, width = 0.75, geom = "crossbar", col = 'black') +
    geom_jitter(position=position_jitter(0.35), size = 1.5) +
    scale_color_manual(values=c("#ebc999",
                                "#cd7700", 
                                "#4d3227"
    )) +
    lemon::facet_rep_grid(~pheno, repeat.tick.labels = T, space = 'free', scales = 'free')+
    ggthemes::theme_hc() + 
    stat_pvalue_manual(global_pm2_dunn_gruczo2, hide.ns = TRUE)+
    labs(
      subtitle = get_test_label(global_pm2_kruskal_gruczo, detailed = TRUE),
      caption = get_pwc_label(global_pm2_dunn_gruczo2)
    )+
    xlab('Phenotype')+
    ylab('Mean Methylation [%]')+
    scale_x_discrete(labels= c('Tn','Tp','Tc'))+
    scale_y_continuous(limits = c(10, 65),
                       breaks = c(0,10,20,30,40,50,60)))
dev.off()

global_pm2_dunn_plasko2 <- add_xy_position(global_pm2_dunn_plasko, x = "group", step.increase = 1)
global_pm2_dunn_plasko2$xmin <- c(1,1,2)
global_pm2_dunn_plasko2$xmax <- c(2,3,3)
global_pm2_dunn_plasko2$y.position <- c(45.6, 49.5, 53.5)
pdf('plasko.pdf', 3,5)
(p2 <- ggplot(global_pm2[global_pm2$pheno=='plasko',], aes(x=group, y=mean_perc_meth, color = group)) + 
    #geom_violin(trim = TRUE) +
    stat_summary(fun = "median", fun.min = "median", fun.max= "median", size= 0.3, width = 0.75, geom = "crossbar", col = 'black') +
    geom_jitter(position=position_jitter(0.4), size = 1.5) + 
    scale_color_manual(values=c("#37a987",
                                "#b7b1d2", 
                                "#4b3d8f"
    )) +
    lemon::facet_rep_grid(~pheno, repeat.tick.labels = T, space = 'free', scales = 'free')+
    ggthemes::theme_hc() + 
    stat_pvalue_manual(global_pm2_dunn_plasko2, hide.ns = TRUE)+
    labs(
      subtitle = get_test_label(global_pm2_kruskal_plasko, detailed = TRUE),
      caption = get_pwc_label(global_pm2_dunn_plasko2),
    )+
    xlab('Phenotype')+
    ylab('Mean Methylation [%]')+
    scale_x_discrete(labels= c('Tn','Tp','Tc'))+
    scale_y_continuous(limits = c(10, 65),
                       breaks = c(0,10,20,30,40,50,60)))
dev.off()

# Prepare for circular manhatann ------------------------------------------
myDiff_plasko__tp_data <- getData(myDiff_plasko__tp) %>% 
  dplyr::select(end, chr, start, qvalue) %>% 
  mutate(qvalue_tp = qvalue) %>% 
  arrange(factor(chr, levels = c(1:22,'X')), chr, start) %>% 
  mutate(SNP = factor(paste0('EndNo_',end)),
         Chromosome = factor(chr),
         Position = factor(start),
         # qvalue_tp = runif(length(myDiff_plasko__tp_data$SNP), 0, 1) # testing values whether plotting works
  ) %>% 
  dplyr::select(SNP, Chromosome, Position, qvalue_tp)

myDiff_plasko__tc_data <- getData(myDiff_plasko__tc) %>% 
  dplyr::select(end, chr, start, qvalue) %>% 
  mutate(qvalue_tc = qvalue) %>% 
  arrange(factor(chr, levels = c(1:22,'X')), chr, start) %>% 
  mutate(SNP = factor(paste0('EndNo_',end)),
         Chromosome = factor(chr),
         Position = factor(start),
         # qvalue_tp = runif(length(myDiff_plasko__tp_data$SNP), 0, 1) # testing values whether plotting works
  ) %>% 
  dplyr::select(SNP, Chromosome, Position, qvalue_tc)

cm_plasko_merged <- full_join(myDiff_plasko__tp_data, myDiff_plasko__tc_data,
                              by = c('SNP', 'Chromosome', 'Position'))

cm_plasko_merged$qvalue_tp[is.na(cm_plasko_merged$qvalue_tp)] <- 1
cm_plasko_merged$qvalue_tc[is.na(cm_plasko_merged$qvalue_tc)] <- 1

cm_plasko_merged$qvalue_tc[cm_plasko_merged$qvalue_tc <= 0] <- 1e-300
cm_plasko_merged$qvalue_tp[cm_plasko_merged$qvalue_tp <= 0] <- 1e-300

sum(cm_plasko_merged$qvalue_tc <= 0)
sum(cm_plasko_merged$qvalue_tc > 1)

sum(cm_plasko_merged$qvalue_tp <= 0)
sum(cm_plasko_merged$qvalue_tp > 1)

# Plot
library(CMplot)
CMplot(
  cm_plasko_merged,
  type = "p",
  plot.type = "c",
  r = 1,
  col = matrix(c("#b7b1d2", 
                 "#4b3d8f"), 2, 1),
  chr.labels = paste("Chr", c(1:22, "X"), sep = ""),
  threshold = c(0.05),
  cir.chr.h = 1,
  amplify = FALSE,
  threshold.lty = c(1, 2),
  threshold.lwd = 2,
  threshold.col = c("red",
                    "blue"),
  LOG10 = TRUE,
  main = 'CpG Methylations in squamous cancer',
  signal.line = 1,
  signal.col = c("red", "green"),
 # chr.den.col = c("darkgreen", "yellow", "red"),
  chr.den.col="black",
  bin.size = 1e6,
  outward = FALSE,
  file = "pdf",
  memo = "plasko_whole",
  dpi = 300,
  ylim = c(0, 320),
  file.output = TRUE,
  verbose = TRUE,
  width = 10,
  height = 10
)


myDiff_gruczo__tp_data <- getData(myDiff_gruczo__tp) %>% 
  dplyr::select(end, chr, start, qvalue) %>% 
  mutate(qvalue_tp = qvalue) %>% 
  arrange(factor(chr, levels = c(1:22,'X')), chr, start) %>% 
  mutate(SNP = factor(paste0('EndNo_',end)),
         Chromosome = factor(chr),
         Position = factor(start),
         # qvalue_tp = runif(length(myDiff_gruczo__tp_data$SNP), 0, 1) # testing values whether plotting works
  ) %>% 
  dplyr::select(SNP, Chromosome, Position, qvalue_tp)

myDiff_gruczo__tc_data <- getData(myDiff_gruczo__tc) %>% 
  dplyr::select(end, chr, start, qvalue) %>% 
  mutate(qvalue_tc = qvalue) %>% 
  arrange(factor(chr, levels = c(1:22,'X')), chr, start) %>% 
  mutate(SNP = factor(paste0('EndNo_',end)),
         Chromosome = factor(chr),
         Position = factor(start),
         # qvalue_tp = runif(length(myDiff_gruczo__tp_data$SNP), 0, 1) # testing values whether plotting works
  ) %>% 
  dplyr::select(SNP, Chromosome, Position, qvalue_tc)

cm_gruczo_merged <- full_join(myDiff_gruczo__tp_data, myDiff_gruczo__tc_data,
                              by = c('SNP', 'Chromosome', 'Position'))

cm_gruczo_merged$qvalue_tp[is.na(cm_gruczo_merged$qvalue_tp)] <- 1
cm_gruczo_merged$qvalue_tc[is.na(cm_gruczo_merged$qvalue_tc)] <- 1

cm_gruczo_merged$qvalue_tc[cm_gruczo_merged$qvalue_tc <= 0] <- 1e-300
cm_gruczo_merged$qvalue_tp[cm_gruczo_merged$qvalue_tp <= 0] <- 1e-300

sum(cm_gruczo_merged$qvalue_tc <= 0)
sum(cm_gruczo_merged$qvalue_tc > 1)

sum(cm_gruczo_merged$qvalue_tp <= 0)
sum(cm_gruczo_merged$qvalue_tp > 1)

# Plot
library(CMplot)
CMplot(
  cm_gruczo_merged,
  type = "p",
  plot.type = "c",
  r = 1,
  col = matrix(c("#cd7700", 
                 "#4d3227"), 2, 1),
  chr.labels = paste("Chr", c(1:22, "X"), sep = ""),
  threshold = c(0.05),
  cir.chr.h = 1,
  amplify = FALSE,
  threshold.lty = c(1, 2),
  threshold.lwd = 2,
  threshold.col = c("red",
                    "blue"),
  LOG10 = TRUE,
  main = 'CpG Methylations in squamous cancer',
  signal.line = 1,
  signal.col = c("red", "green"),
  # chr.den.col = c("darkgreen", "yellow", "red"),
  chr.den.col="black",
  bin.size = 1e6,
  outward = FALSE,
  file = "pdf",
  memo = "gruczo_whole",
  dpi = 300,
  ylim = c(0, 320),
  file.output = TRUE,
  verbose = TRUE,
  width = 10,
  height = 10
)

# Difference between significant CpG --------------------------------------
# apply cutoffs
myDiff_plasko__tc_25p <- getMethylDiff(myDiff_plasko__tc, 
                                       difference=25, 
                                       qvalue=0.05)
myDiff_plasko__tp_25p <- getMethylDiff(myDiff_plasko__tp, 
                                       difference=25, 
                                       qvalue=0.05)
myDiff_gruczo__tc_25p <- getMethylDiff(myDiff_gruczo__tc, 
                                       difference=25, 
                                       qvalue=0.05)
myDiff_gruczo__tp_25p <- getMethylDiff(myDiff_gruczo__tp, 
                                       difference=25, 
                                       qvalue=0.05)


# check it for chromosomes
myDiff_plasko__tc_25p_chrom <- diffMethPerChr(
  myDiff_plasko__tc_25p,
  plot = T,
  qvalue.cutoff = 0.05,
  meth.cutoff = 25
)
myDiff_plasko__tp_25p_chrom <- diffMethPerChr(
  myDiff_plasko__tp_25p,
  plot = F,
  qvalue.cutoff = 0.05,
  meth.cutoff = 25
)
myDiff_gruczo__tc_25p_chrom <- diffMethPerChr(
  myDiff_gruczo__tc_25p,
  plot = F,
  qvalue.cutoff = 0.05,
  meth.cutoff = 25
)
myDiff_gruczo__tp_25p_chrom <- diffMethPerChr(
  myDiff_gruczo__tp_25p,
  plot = F,
  qvalue.cutoff = 0.05,
  meth.cutoff = 25
)


# Venn Diagram ------------------------------------------------------------
myDiff_plasko__tc_25p_venn <- getData(myDiff_plasko__tc_25p) %>% 
  dplyr::mutate(ID = paste0(chr,'_',start,'_',end)) %>% 
  dplyr::select(ID, meth.diff) %>% 
  dplyr::rename(logFC = meth.diff)

myDiff_plasko__tp_25p_venn <- getData(myDiff_plasko__tp_25p) %>% 
  dplyr::mutate(ID = paste0(chr,'_',start,'_',end)) %>% 
  dplyr::select(ID, meth.diff) %>% 
  dplyr::rename(logFC = meth.diff)

myDiff_gruczo__tc_25p_venn <- getData(myDiff_gruczo__tc_25p) %>% 
  dplyr::mutate(ID = paste0(chr,'_',start,'_',end)) %>% 
  dplyr::select(ID, meth.diff) %>% 
  dplyr::rename(logFC = meth.diff)

myDiff_gruczo__tp_25p_venn <- getData(myDiff_gruczo__tp_25p) %>% 
  dplyr::mutate(ID = paste0(chr,'_',start,'_',end)) %>% 
  dplyr::select(ID, meth.diff) %>% 
  dplyr::rename(logFC = meth.diff)

venn_list <- list(
  p_tc = myDiff_plasko__tc_25p_venn$ID, 
  p_tp = myDiff_plasko__tp_25p_venn$ID, 
  g_tc = myDiff_gruczo__tc_25p_venn$ID,
  g_tp = myDiff_gruczo__tp_25p_venn$ID
)
library(VennDiagram)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

pdf('v2_venn.pdf',10,10)

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

# Chromo map --------------------------------------------------------------
setwd('/Users/andrzejeljaszewicz/Adrian/Omics/MOBIT/mobit_rrbs')

myDiff_plasko__tc_25p_data <- getData(myDiff_plasko__tc_25p) %>% 
  mutate(ID = paste0(chr,'_',start)) %>% 
  dplyr::select(ID, chr, start, end, meth.diff)
readr::write_delim(myDiff_plasko__tc_25p_data, 'myDiff_plasko__tc_25p_data.csv',
                   col_names = F, delim = "\t")

myDiff_plasko__tp_25p_data <- getData(myDiff_plasko__tp_25p) %>% 
  mutate(ID = paste0(chr,'_',start)) %>% 
  dplyr::select(ID, chr, start, end, meth.diff)
readr::write_delim(myDiff_plasko__tp_25p_data, 'myDiff_plasko__tp_25p_data.csv',
                   col_names = F, delim = "\t")

myDiff_gruczo__tc_25p_data <- getData(myDiff_gruczo__tc_25p) %>% 
  mutate(ID = paste0(chr,'_',start)) %>% 
  dplyr::select(ID, chr, start, end, meth.diff)
readr::write_delim(myDiff_gruczo__tc_25p_data, 'myDiff_gruczo__tc_25p_data.csv',
                   col_names = F, delim = "\t")

myDiff_gruczo__tp_25p_data <- getData(myDiff_gruczo__tp_25p) %>% 
  mutate(ID = paste0(chr,'_',start)) %>% 
  dplyr::select(ID, chr, start, end, meth.diff)
readr::write_delim(myDiff_gruczo__tp_25p_data, 'myDiff_gruczo__tp_25p_data.csv',
                   col_names = F, delim = "\t")
# chromosome files
chr_file = "/Users/andrzejeljaszewicz/Adrian/Omics/MOBIT/mobit_rrbs/chromosomes.txt"
# annotation files
myDiff_plasko__tc_25p_data_anno_file = "myDiff_plasko__tc_25p_data.csv"
myDiff_plasko__tp_25p_data_anno_file = "myDiff_plasko__tp_25p_data.csv"
myDiff_gruczo__tc_25p_data_anno_file = "myDiff_gruczo__tc_25p_data.csv"
myDiff_gruczo__tp_25p_data_anno_file = "myDiff_gruczo__tp_25p_data.csv"

data_colors <- list(c('blue','white','red'))
library(chromoMap)

chromoMap(chr_file,
          myDiff_plasko__tc_25p_data_anno_file,
          data_based_color_map = T,
          legend = c(TRUE),
          data_colors = data_colors,
          export.options = T,
          numeric.domain = c(-50,50))

chromoMap(chr_file,
          myDiff_plasko__tp_25p_data_anno_file,
          data_based_color_map = T,
          legend = c(TRUE),
          data_colors = data_colors,
          export.options = T,
          numeric.domain = c(-50,50))

chromoMap(chr_file,
          myDiff_gruczo__tc_25p_data_anno_file,
          data_based_color_map = T,
          legend = c(TRUE),
          data_colors = data_colors,
          export.options = T,
          numeric.domain = c(-50,50))

chromoMap(chr_file,
          myDiff_gruczo__tp_25p_data_anno_file,
          data_based_color_map = T,
          legend = c(TRUE),
          data_colors = data_colors,
          export.options = T,
          numeric.domain = c(-50,50))

# Annotation --------------------------------------------------------------
library(annotatr)
annots = c('hg38_basicgenes', 'hg38_cpgs', 'hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)

myDiff_plasko__tc_25p_dm_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_plasko__tc_25p) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_plasko__tc_25p))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)
myDiff_plasko__tc_25p_dm_annotated_f <- data.frame(myDiff_plasko__tc_25p_dm_annotated) %>% distinct(SNP_ID, .keep_all = T) %>% 
  dplyr::filter(!is.na(annot.symbol)) %>% 
  group_by(annot.symbol) %>% 
  filter(n()>= 4) %>% 
  ungroup() 

myDiff_plasko__tp_25p_dm_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_plasko__tp_25p) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_plasko__tp_25p))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)
myDiff_plasko__tp_25p_dm_annotated_f <- data.frame(myDiff_plasko__tp_25p_dm_annotated) %>% distinct(SNP_ID, .keep_all = T) %>% 
  dplyr::filter(!is.na(annot.symbol)) %>% 
  group_by(annot.symbol) %>% 
  filter(n()>= 4) %>% 
  ungroup() 

myDiff_gruczo__tc_25p_dm_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_gruczo__tc_25p) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_gruczo__tc_25p))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)
myDiff_gruczo__tc_25p_dm_annotated_f <- data.frame(myDiff_gruczo__tc_25p_dm_annotated) %>% distinct(SNP_ID, .keep_all = T) %>% 
  dplyr::filter(!is.na(annot.symbol)) %>% 
  group_by(annot.symbol) %>% 
  filter(n()>= 4) %>% 
  ungroup() 

myDiff_gruczo__tp_25p_dm_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_gruczo__tp_25p) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_gruczo__tp_25p))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)
myDiff_gruczo__tp_25p_dm_annotated_f <- data.frame(myDiff_gruczo__tp_25p_dm_annotated) %>% distinct(SNP_ID, .keep_all = T) %>% 
  dplyr::filter(!is.na(annot.symbol)) %>% 
  group_by(annot.symbol) %>% 
  filter(n()>= 4) %>% 
  ungroup() 


myDiff_plasko__tc_25p_dm_annotated_summary = summarize_annotations(
  annotated_regions = myDiff_plasko__tc_25p_dm_annotated,
  quiet = TRUE) %>% 
  mutate(group = 'plasko_Tc')

myDiff_plasko__tp_25p_dm_annotated_summary = summarize_annotations(
  annotated_regions = myDiff_plasko__tp_25p_dm_annotated,
  quiet = TRUE) %>% 
  mutate(group = 'plasko_Tp')

myDiff_gruczo__tc_25p_dm_annotated_summary = summarize_annotations(
  annotated_regions = myDiff_gruczo__tc_25p_dm_annotated,
  quiet = TRUE) %>% 
  mutate(group = 'gruczo_Tc')

myDiff_gruczo__tp_25p_dm_annotated_summary = summarize_annotations(
  annotated_regions = myDiff_gruczo__tp_25p_dm_annotated,
  quiet = TRUE) %>% 
  mutate(group = 'gruczo_Tp')

summary_merged <- myDiff_plasko__tc_25p_dm_annotated_summary %>% 
  full_join(myDiff_plasko__tp_25p_dm_annotated_summary) %>% 
  full_join(myDiff_gruczo__tc_25p_dm_annotated_summary) %>% 
  full_join(myDiff_gruczo__tp_25p_dm_annotated_summary) %>% 
  arrange(n)
summary_merged$annot.type <- factor(summary_merged$annot.type,
                                    levels = unique(summary_merged$annot.type))

summary_merged2 <- summary_merged %>% 
  group_by(annot.type, group) %>% 
  mutate(perc = n/sum(n))
pdf('annotations.pdf', 7,7)
ggplot(summary_merged2,
       aes(x = annot.type,
           y = n,
           fill = group)) +
  geom_bar(stat = "identity", width = 0.65) +
  ggthemes::theme_hc() +
  scale_fill_manual(values=c(
                             "#4d3227",
                             "#cd7700","#4b3d8f",
                             "#b7b1d2"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()
dev.off()
ggplot(summary_merged,
       aes(x = annot.type,
           y = n,
           fill = group)) +
  geom_bar(stat = "identity", width = 0.65) +
  ggthemes::theme_hc() +
  scale_fill_manual(values=c("#B38B59", 
                             "#5a4b29",
                             "#748cac",
                             "#2c3c4c"
  )) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()

# For all 

myDiff_plasko__tc_all_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_plasko__tc) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_plasko__tc))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)

myDiff_plasko__tp_all_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_plasko__tp) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_plasko__tp))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)

myDiff_gruczo__tc_all_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_gruczo__tc) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_gruczo__tc))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)

myDiff_gruczo__tp_all_annotated = annotate_regions(
  regions = regioneR::toGRanges(getData(myDiff_gruczo__tp) %>% dplyr::mutate(SNP_ID = paste0('SNP_', 1:nrow(getData(myDiff_gruczo__tp))))),
  annotations = diffloop::rmchr(annotations),
  ignore.strand = TRUE,
  quiet = FALSE)


# GO over-rep -------------------------------------------------------------
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
# create element as an Sys.time() object - 1648568306
(the_seed = 1648568306 %% 100000) # 50700

# Run and set seed each time for rep
plasko_tc_signif_cpg <- myDiff_plasko__tc_25p_dm_annotated_f$annot.gene_id
plasko_tc_signif_cpg_all <- data.frame(myDiff_plasko__tc_all_annotated)$annot.gene_id

plasko_tp_signif_cpg <- myDiff_plasko__tp_25p_dm_annotated_f$annot.gene_id
plasko_tp_signif_cpg_all <- data.frame(myDiff_plasko__tp_all_annotated)$annot.gene_id

gruczo_tc_signif_cpg <- myDiff_gruczo__tc_25p_dm_annotated_f$annot.gene_id
gruczo_tc_signif_cpg_all <- data.frame(myDiff_gruczo__tc_all_annotated)$annot.gene_id

gruczo_tp_signif_cpg <- myDiff_gruczo__tp_25p_dm_annotated_f$annot.gene_id
gruczo_tp_signif_cpg_all <- data.frame(myDiff_gruczo__tp_all_annotated)$annot.gene_id


set.seed(the_seed)
plasko_tc_go_ora_rrbs <- enrichGO(gene          = plasko_tc_signif_cpg,
                                  universe      = plasko_tc_signif_cpg_all,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)
set.seed(the_seed)
plasko_tp_go_ora_rrbs <- enrichGO(gene          = plasko_tp_signif_cpg,
                                  universe      = plasko_tp_signif_cpg_all,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)
set.seed(the_seed)
gruczo_tc_go_ora_rrbs <- enrichGO(gene          = gruczo_tc_signif_cpg,
                                  universe      = gruczo_tc_signif_cpg_all,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)
set.seed(the_seed)
gruczo_tp_go_ora_rrbs <- enrichGO(gene          = gruczo_tp_signif_cpg,
                                  universe      = gruczo_tp_signif_cpg_all,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)

# KEGG over-rep -----------------------------------------------------------
set.seed(the_seed)
plasko_tc_kegg_ora_rrbs <- enrichKEGG(gene          = plasko_tc_signif_cpg,
                                      universe      = plasko_tc_signif_cpg_all,
                                      minGSSize     = 10,
                                      maxGSSize     = 500,
                                      organism     = 'hsa',
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
set.seed(the_seed)
plasko_tp_kegg_ora_rrbs <- enrichKEGG(gene          = plasko_tp_signif_cpg,
                                      universe      = plasko_tp_signif_cpg_all,
                                      minGSSize     = 10,
                                      maxGSSize     = 500,
                                      organism     = 'hsa',
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
set.seed(the_seed)
gruczo_tc_kegg_ora_rrbs <- enrichKEGG(gene          = gruczo_tc_signif_cpg,
                                      universe      = gruczo_tc_signif_cpg_all,
                                      minGSSize     = 10,
                                      maxGSSize     = 500,
                                      organism     = 'hsa',
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)
set.seed(the_seed)
gruczo_tp_kegg_ora_rrbs <- enrichKEGG(gene          = gruczo_tp_signif_cpg,
                                      universe      = gruczo_tp_signif_cpg_all,
                                      minGSSize     = 10,
                                      maxGSSize     = 500,
                                      organism     = 'hsa',
                                      pvalueCutoff = 1,
                                      qvalueCutoff = 1)

# Wiki-Pathways Ora -------------------------------------------------------
set.seed(the_seed)
plasko_tc_wp_ora_rrbs <- enrichWP(gene =        plasko_tc_signif_cpg, 
                                  universe      = plasko_tc_signif_cpg_all,
                                  organism      = "Homo sapiens",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1) 
set.seed(the_seed)
plasko_tp_wp_ora_rrbs <- enrichWP(gene =        plasko_tp_signif_cpg, 
                                  universe      = plasko_tp_signif_cpg_all,
                                  organism      = "Homo sapiens",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1) 
set.seed(the_seed)
gruczo_tc_wp_ora_rrbs <- enrichWP(gene =        gruczo_tc_signif_cpg, 
                                  universe      = gruczo_tc_signif_cpg_all,
                                  organism      = "Homo sapiens",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1) 
set.seed(the_seed)
gruczo_tp_wp_ora_rrbs <- enrichWP(gene =        gruczo_tp_signif_cpg, 
                                  universe      = gruczo_tp_signif_cpg_all,
                                  organism      = "Homo sapiens",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  pvalueCutoff  = 1,
                                  qvalueCutoff  = 1) 


# Reactome ORA ------------------------------------------------------------

set.seed(the_seed)
plasko_tc_reactome_ora_rrbs <- enrichPathway(gene =        plasko_tc_signif_cpg, 
                                             universe      = plasko_tc_signif_cpg_all,
                                             organism = "human",
                                             pvalueCutoff = 1,
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 1,
                                             minGSSize = 10,
                                             maxGSSize = 500,
                                             readable = TRUE)

set.seed(the_seed)
plasko_tp_reactome_ora_rrbs <- enrichPathway(gene =        plasko_tp_signif_cpg, 
                                             universe      = plasko_tp_signif_cpg_all,
                                             organism = "human",
                                             pvalueCutoff = 1,
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 1,
                                             minGSSize = 10,
                                             maxGSSize = 500,
                                             readable = TRUE)

set.seed(the_seed)
gruczo_tc_reactome_ora_rrbs <- enrichPathway(gene =        gruczo_tc_signif_cpg, 
                                             universe      = gruczo_tc_signif_cpg_all,
                                             organism = "human",
                                             pvalueCutoff = 1,
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 1,
                                             minGSSize = 10,
                                             maxGSSize = 500,
                                             readable = TRUE)

set.seed(the_seed)
gruczo_tp_reactome_ora_rrbs <- enrichPathway(gene =        gruczo_tp_signif_cpg, 
                                             universe      = gruczo_tp_signif_cpg_all,
                                             organism = "human",
                                             pvalueCutoff = 1,
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 1,
                                             minGSSize = 10,
                                             maxGSSize = 500,
                                             readable = TRUE)


# Disease Ontology --------------------------------------------------------
set.seed(the_seed)
plasko_tc_do_ora_rrbs <- enrichDO(gene =        plasko_tc_signif_cpg, 
                                  universe      = plasko_tc_signif_cpg_all,
                                  ont           = "DO",
                                  pvalueCutoff  = 1,
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)

set.seed(the_seed)
plasko_tp_do_ora_rrbs <- enrichDO(gene =        plasko_tp_signif_cpg, 
                                  universe      = plasko_tp_signif_cpg_all,
                                  ont           = "DO",
                                  pvalueCutoff  = 1,
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)

set.seed(the_seed)
gruczo_tc_do_ora_rrbs <- enrichDO(gene =        gruczo_tc_signif_cpg, 
                                  universe      = gruczo_tc_signif_cpg_all,
                                  ont           = "DO",
                                  pvalueCutoff  = 1,
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)

set.seed(the_seed)
gruczo_tp_do_ora_rrbs <- enrichDO(gene =        gruczo_tp_signif_cpg, 
                                  universe      = gruczo_tp_signif_cpg_all,
                                  ont           = "DO",
                                  pvalueCutoff  = 1,
                                  pAdjustMethod = "BH",
                                  minGSSize     = 10,
                                  maxGSSize     = 500,
                                  qvalueCutoff  = 1,
                                  readable      = TRUE)


# Network of Cancer Gene --------------------------------------------------
set.seed(the_seed)
plasko_tc_ncg_ora_rrbs <- enrichNCG(gene =        plasko_tc_signif_cpg, 
                                    universe      = plasko_tc_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

set.seed(the_seed)
plasko_tp_ncg_ora_rrbs <- enrichNCG(gene =        plasko_tp_signif_cpg, 
                                    universe      = plasko_tp_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

set.seed(the_seed)
gruczo_tc_ncg_ora_rrbs <- enrichNCG(gene =        gruczo_tc_signif_cpg, 
                                    universe      = gruczo_tc_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

set.seed(the_seed)
gruczo_tp_ncg_ora_rrbs <- enrichNCG(gene =        gruczo_tp_signif_cpg, 
                                    universe      = gruczo_tp_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

# Network of Disease Gene Network -----------------------------------------
set.seed(the_seed)
plasko_tc_dgn_ora_rrbs <- enrichDGN(gene =        plasko_tc_signif_cpg, 
                                    universe      = plasko_tc_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

set.seed(the_seed)
plasko_tp_dgn_ora_rrbs <- enrichDGN(gene =        plasko_tp_signif_cpg, 
                                    universe      = plasko_tp_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

set.seed(the_seed)
gruczo_tc_dgn_ora_rrbs <- enrichDGN(gene =        gruczo_tc_signif_cpg, 
                                    universe      = gruczo_tc_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)

set.seed(the_seed)
gruczo_tp_dgn_ora_rrbs <- enrichDGN(gene =        gruczo_tp_signif_cpg, 
                                    universe      = gruczo_tp_signif_cpg_all,
                                    pvalueCutoff  = 1,
                                    pAdjustMethod = "BH",
                                    minGSSize     = 10,
                                    maxGSSize     = 500,
                                    qvalueCutoff  = 1,
                                    readable      = TRUE)





# phyper ------------------------------------------------------------------

phyper(
  length(intersect(unique(myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol), unique(myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol))) - 1,
       length(unique(myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol)),
       length(unique(c(data.frame(myDiff_gruczo__tc_all_annotated)$annot.symbol, data.frame(myDiff_gruczo__tp_all_annotated)$annot.symbol)))-length(unique(myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol)),
       length(unique(myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol)),
       lower.tail = FALSE,
       log.p = FALSE)


phyper(length(intersect(unique(myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol),unique(myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol)))-1,
       length(unique(myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol)),
       length(unique(c(data.frame(myDiff_plasko__tc_all_annotated)$annot.symbol,data.frame(myDiff_plasko__tp_all_annotated)$annot.symbol)))-length(unique(myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol)),
       length(unique(myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol)),
       lower.tail = FALSE,
       log.p = FALSE)

openxlsx::write.xlsx(data.frame(intersect(unique(myDiff_plasko__tc_25p_dm_annotated_f$annot.symbol), unique(myDiff_plasko__tp_25p_dm_annotated_f$annot.symbol))), 'common_plasko_fixed.xlsx')
openxlsx::write.xlsx(data.frame(intersect(unique(myDiff_gruczo__tc_25p_dm_annotated_f$annot.symbol), unique(myDiff_gruczo__tp_25p_dm_annotated_f$annot.symbol))), 'common_gruczo_fixed.xlsx')


