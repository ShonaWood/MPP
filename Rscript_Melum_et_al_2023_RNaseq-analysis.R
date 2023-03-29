#### Script to attach to submission Melum et al. 2023

## Performed with:
# R version 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting"
# Copyright (C) 2022 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)

# in RStudio 2023.03.0+386 "Cherry Blossom" Release (3c53477afb13ab959aeb5b34df1f10c237b256c3, 2023-03-09) for Windows

###########################

## Required packages 
library(factoextra)
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(VennDiagram)
library(EnhancedVolcano)
library(readxl)
library(ggpubr)
library(janitor)
library(ggsci)
library(scales)
library(edgeR)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(matrixStats)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(patchwork) 
library(forcats)
library(viridis)

#### Data filtering and statistical analysis 

# read count data:

df <- read_delim("Preproccessed_data/merged_table_counts_2_20220314.txt")

names(df)<-gsub("\\_","",names(df))

genes_mmus <- read_delim("External_data/Genes_mmus_info.txt", delim = " ")


df_geneID <- df[,1]

df_full<- df_geneID %>% 
  left_join(genes_mmus,by=c("Geneid"="ensembl_gene_id"),copy=F,keep=F)

df_full <- df_full %>% 
  distinct(Geneid, .keep_all = TRUE)

# Convert column 1 to rownames
df2 <- df %>% 
  remove_rownames() %>% 
  column_to_rownames(var="Geneid")

#change order, so LP, LPIP,  SPIP, SP to keep it as want to display it.  

df3 <- df2 %>% 
  relocate(6:11,1:5,12:16,17:21)

# Remove SPIP19 due to poor sequencing quality

df4 <- df3 %>% 
  dplyr::select(-SPIP19)

sample_names <- colnames(df4)

mt2 <- as.matrix(df4)

sample_names <- colnames(df4)

#design matrix/ group 
group2<- factor(c(rep("LP", 6), rep("LPIP", 5),rep("SPIP",4),rep("SP",5)))

DGE_glm<-DGEList(mt2, group = group2, samples = sample_names)

# for looking at all groups referenced to LP (intercept)
design2<-model.matrix(~group2)

#for being able to contrast all groups "real" values:

design3 <- model.matrix(~ 0 + group2, data = DGE_glm$samples)

colnames(design3) <- levels(DGE_glm$samples$group)

cpm_unfiltered <- cpm(DGE_glm, log=FALSE, normalized.lib.sizes = TRUE)

tb_cpm_unfiltered <- as_tibble(cpm_unfiltered, rownames="ENSEMBLID")

keep <- rowMedians(cpm(DGE_glm, log=FALSE, normalized.lib.sizes = TRUE))>=1

DGE_glm <- DGE_glm[keep, , keep.lib.sizes=FALSE]

# Keep only genes that have median cpm >=1 - 13345 

# Make cpm - log2cpm dataframes 

### make cpm, with normalized lib sizes to take into account the differences between sequencing depths:

log2cpm <- cpm(DGE_glm, log = TRUE, normalized.lib.sizes = TRUE)

logcounts_gr <- cpmByGroup(DGE_glm, log=TRUE, normalized.lib.sizes = TRUE)

cpm_gr <- cpmByGroup(DGE_glm, log=FALSE, normalized.lib.sizes = TRUE)

cpm <- cpm(DGE_glm, log=FALSE, normalized.lib.sizes = TRUE)


#Calculate normality factors
DGE_glm <-calcNormFactors(DGE_glm, method = "TMM")

# estimate dispersion
DGE_glm <- estimateDisp(DGE_glm,design3)


fit<- glmQLFit(DGE_glm, design = design2)

fit2 <- glmQLFit(DGE_glm, design = design3)


# After glmQLFit, can perform glmWLFTest to run the quasi-likelihood (QL) F-test
#Try to run for whole dataset, ie LPIP vs LP, SP vs LP, SPIP vs LP, this will indicate which genes that differ most in all groups contrasted with LP. 



#test all groups vs LP, since I now have made ~group to the design matrix I have made an intercept, this means that I can specify coef=2:4, this will compare LPIP, SP and SPIP to the LP baseline

qlf2_all<- glmQLFTest(fit, coef = 2:4)


Toptabs_all_vs_LP <-topTags(qlf2_all, n=Inf, adjust.method = "fdr", sort.by = "none", p.value = Inf)

tb_qlf_all_vs_LP_2 <-Toptabs_all_vs_LP$table

tb_qlf_all_vs_LP <- tb_qlf_all_vs_LP_2 %>% 
  rownames_to_column(var="ENSEMBLID")%>% 
  left_join(ID_conv, by=c("ENSEMBLID"="ENSEMBL")) %>% 
  relocate(1,9:11) %>% 
  distinct(ENSEMBLID, .keep_all = TRUE)

tb_qlf_all_vs_LP
write.table(tb_qlf_all_vs_LP, file="output/glm_test_All_vs_LP_mediancpm1.txt")

summary(decideTestsDGE(qlf2_all))

# Genelist to use for subset of genes in heatmap, logCPM >= 0, FDR<0.001, logFC 0.5 in all compared to LP. 

tb_all_vs_LP_fdr <-topTags(qlf2_all, n=Inf, adjust.method = "fdr", sort.by = "none", p.value = 0.001)

tb_all_vs_LP_fdr <- tb_all_vs_LP_fdr$table

tb_all_vs_LP_fdr_logCPM <- tb_all_vs_LP_fdr %>% 
  rownames_to_column(var="Geneid") %>% 
  filter(logCPM>=0)

## continue with other contrasts:

qlf.3vs1_volcano <- glmQLFTest(fit2, contrast=c(-1,0,1,0), aveLogCPM(fit2))

Toptabs_SPvsLP <-topTags(qlf.3vs1_volcano,n=Inf, adjust.method = "fdr", sort.by = "logFC")

tb_qlf_SP_vs_LP <-Toptabs_SPvsLP$table

tb_qlf_SP_LP <- tb_qlf_SP_vs_LP %>% 
  rownames_to_column(var="ENSEMBLID")%>% 
  left_join(ID_conv, by=c("ENSEMBLID"="ENSEMBL")) %>% 
  relocate(1,7:9) %>% 
  distinct(ENSEMBLID, .keep_all = TRUE)

write.table(tb_qlf_SP_LP, file="output/glm_test_SP_vs_LP_anno_median_cpm1.txt")

tb_qlf_SP_vs_LP_vol <- tb_qlf_SP_LP %>% 
  filter(logCPM>0)

# continue with LPIP vs LP

qlf.2vs1_volcano <- glmQLFTest(fit2, contrast=c(-1,1,0,0), aveLogCPM(fit2))

Toptabs_LPIPvsLP <-topTags(qlf.2vs1_volcano,n=Inf, adjust.method = "fdr", sort.by = "logFC")

tb_qlf_LPIPvsLP <-Toptabs_LPIPvsLP$table

tb_qlf_LPIP_LP_an <- tb_qlf_LPIPvsLP %>% 
  rownames_to_column(var="ENSEMBLID")%>% 
  left_join(ID_conv, by=c("ENSEMBLID"="ENSEMBL")) %>% 
  relocate(1,7:9) %>% 
  distinct(ENSEMBLID, .keep_all = TRUE)

write.table(tb_qlf_LPIP_LP_an, file="output/glm_test_LPIP_LP_anno_median_cpm1.txt")

tb_qlf_LPIP_LP_vol <- tb_qlf_LPIP_LP_an %>% 
  filter(logCPM>0)

# SPIP vs SP

qlf.4vs3_volcano <- glmQLFTest(fit2, contrast=c(0,0,-1,1), aveLogCPM(fit2))

Toptabs_SPIPvsSP <-topTags(qlf.4vs3_volcano,n=Inf, adjust.method = "fdr", sort.by = "logFC")

tb_qlf_SPIPvsSP  <-Toptabs_SPIPvsSP$table

tb_qlf_SPIP_SP  <- tb_qlf_SPIPvsSP  %>% 
  rownames_to_column(var="ENSEMBLID")%>% 
  left_join(ID_conv, by=c("ENSEMBLID"="ENSEMBL")) %>% 
  relocate(1,7:9) %>% 
  distinct(ENSEMBLID, .keep_all = TRUE)

write.table(tb_qlf_SPIP_SP , file="output/glm_test_SPIP_SP_anno_median_cpm1.txt")

tb_qlf_SPIP_SP_vol <- tb_qlf_SPIP_SP %>% 
  filter(logCPM>0)

######

### Figure 1
# 1A - Study design
# 1B- Testis weight 

plot_testis2<- df_phys2 %>%
  mutate(Group= fct_relevel(Group,"LP", "LPIP", "SPIP", "SP")) %>% 
  ggplot(aes(x=Group, y= `Testisweight(mg)`))+ 
  geom_dotplot(binaxis= "y",
               stackdir = "center",
               dotsize = 1,
               fill = 1)+
  scale_shape_manual(values = c(LP= 2, LPIP= 17, SPIP=16, SP= 16))+
  theme_pubr()+
  stat_summary( fun.y=mean, geom="crossbar", size=0.5, width=0.5,  color="gray")+
  xlab("Groups") +
  ylab("mg")+
  ylim(0,400)+
  theme(axis.line = element_line(colour = 'black', size = 1.75),
        text = element_text(size = 30),
        axis.ticks = element_line(colour = 'black', size = 1.75),
        axis.ticks.length = unit(.25, "cm"))


plot_testis2


# 1C- Z-score body weight
Created in Graphpad prism 8

# 1F - Average CPM per cluster

plot_all_cpm_clust_abs <- ggerrorplot(data= df_tot_cpm_clust_pres_abs ,
                                      x= "Clustervis",
                                      y= "total_cpm",
                                      color= "Clustercol",
                                      desc_stat = "mean_ci",
                                      add= "mean",
                                      error.plot = "errorbar",
                                      combine= FALSE,
                                      merge= FALSE)



plot_all_cpm_clust_abs_ed2 <-ggpar(plot_all_cpm_clust_abs,
                                   xlab= "Cell specific clusters and remaining dataset",
                                   ylab= "Average total cpm",
                                   legend = "none",
                                   legend.title = "Cell types & dataset",
                                   x.text.angle = 45,
                                   axis.line = element_line(colour = 'black', size = 1.5),
                                   text = element_text(size = 20))

plot_all_cpm_clust_abs_ed2


ggsave("output/Fig1_All_tot_cpm_221010.pdf", plot= plot_all_cpm_clust_abs_ed2, dpi=300, width = 28, height = 14, units=c("cm") )

eds3 <- plot_all_cpm_clust_abs_ed2+
  theme(axis.line = element_line(colour = 'black', size = 1.5),
        text = element_text(size = 20))




ggsave("output/Fig1_All_tot_cpm_default_linethick_221010.pdf", plot=eds3, dpi=300, width = 28, height = 14, units=c("cm") )

# 1G- barchart 
Created in graphpad prism 8

# 2A - PCA analysis

mt_cpm_filt <- tb_cpm %>% 
  column_to_rownames(var="ENSEMBLID") %>% 
  as.matrix()

metadata<-data.frame(row.names = colnames(mt_cpm_filt))

metadata$group <- group2


PCA_mt_filt <-PCAtools::pca(mt_cpm_filt, metadata=metadata, removeVar = 0.1)


## BW and no individual labels:

biplot_pca_filt <-PCAtools::biplot(PCA_mt_filt,
                                   colby = "group", 
                                   colkey = c(
                                     LPIP="black",
                                     SPIP="grey", 
                                     SP="black"),
                                   shape= "group", 
                                   shapekey = c(LP= 2, LPIP= 17, SPIP=16, SP= 16),
                                   pointSize = 5,
                                   lab = NULL)
biplot_pca_filt

ggsave("output/biplot_pca_BW_nolabel_on_cpm_filt_data_221014.pdf", 
       plot=biplot_pca_filt,
       dpi=300)


# 2B - Volcano plot SPvsLP

keyvals.color2 <- ifelse(
  tb_qlf_SP_vs_LP_vol$logFC< 0 & tb_qlf_SP_vs_LP_vol$FDR < 0.001 , 'dodgerblue3',
  ifelse(tb_qlf_SP_vs_LP_vol$logFC > 0 & tb_qlf_SP_vs_LP_vol$FDR < 0.001, 'tomato3',
         'darkgreen'))

keyvals.color2[is.na(keyvals.color2)] <- 'darkgreen'
names(keyvals.color2)[keyvals.color2 == 'darkgreen'] <- 'FDR > 0.001'
names(keyvals.color2)[keyvals.color2== "dodgerblue3"] <- "Upregulated LP"
names(keyvals.color2)[keyvals.color2== "tomato3"] <- "Upregulated SP"




SP_LP_vol <-EnhancedVolcano(tb_qlf_SP_vs_LP_vol,
                            lab = tb_qlf_SP_vs_LP_vol$SYMBOL,
                            selectLab = c("Dio3", "Dio2", "Gpr50"),
                            labSize = 3,
                            x = 'logFC',
                            y = 'FDR',
                            xlab = bquote(~Log[2]~ 'fold change'),
                            ylab= bquote(~-Log[10]~ 'FDR'),
                            title = 'SP vs LP',
                            pCutoff = 0.001,
                            FCcutoff = 0,
                            pointSize = 2,
                            colCustom = keyvals.color2,
                            colAlpha = 0.9,
                            legendPosition = 'top',
                            legendLabSize = 15,
                            legendIconSize = 5.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.75,
                            colConnectors = 'black',
                            arrowheads = FALSE,
                            gridlines.major = TRUE,
                            gridlines.minor = FALSE,
                            border = 'partial',
                            borderWidth = 1.5,
                            borderColour = 'black',
                            caption = bquote(~Log[2]~ "fold change cutoff, 0; FDR cutoff, 0.001; total = 13345 variables"),
                            captionLabSize = 10,
                            xlim = c(-10,10),
                            ylim =c(0,17))



ggsave("output/Volcano_SP_LP.pdf", plot=SP_LP_vol, dpi = 300, height = 7, width = 7 )
# 2C - Go term enrichment LP
Up_LP_GO <- read.csv("up_lp_go-fdr0.01-manual-reduncancy-2.csv", header =TRUE)
Up_LP_GO %>% sample_n(5)
Up_LP_GO %>%
  ggplot(aes(x= Generatio, y= reorder(Description, Generatio),  size =Enrichment, colour = log10FDR ))+
  geom_point()+
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(2,10), oob = scales::squish, name = '-Log10FDR')+
  scale_size(limits = c(1,6), range = c(1,20))+
  theme_bw()+
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 
  
# 2D - Dot plot of expression for cilia genes

gene_cilia <- read_csv("cilia-final.csv")
gene_cilia %>% sample_n(5)
markers <- gene_cilia$Gene %>% unique()
gene_cilia %>% filter(Gene %in% markers) %>%
  ggplot(aes(x=group, y= Gene, color = log2counts, size = normcounts))+
  geom_point()+
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,10), oob = scales::squish, name = 'log2count')+
  scale_size(range = c(-1,5))+
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

# 3A - Barchart 

# 3B - Barchart 

# 3C - Lollipop plot enriched SP
Up_SP_GO <- read.csv("up_SP_go-fdr0.01-manual-reduncancy.csv", header =TRUE)
Up_SP_GO %>% sample_n(5)
Up_SP_GO %>%
  ggplot(aes(x= Generatio, y= reorder(Description, Generatio),  size =Enrichment, colour = log10FDR ))+
  geom_point()+
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(2,10), oob = scales::squish, name = '-Log10FDR')+
  scale_size(limits = c(1,6), range = c(1,10))+
  theme_bw()+
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

# 3D - dotplot Ascl1 

plot_Ascl1 <- tb_cpm_long_gr %>% 
  filter(SYMBOL == "Ascl1") %>% 
  ggdotplot(x= "Group",
            y= "Cpm",
            color= "black",
            add= "mean_ci",
            combine= FALSE,
            merge= FALSE)+
  ylim(0,40)

plot_Ascl1

ggsave("output/Ascl1_BW_TMM_cpm_20221017.pdf", plot_Ascl1, dpi=300)

# 3E - Done in graphpad prism 8

# 3F - Barchart Done in graphpad prism 8

# 4A - Volcano plot SPIP vs SP

keyvals.color3 <- ifelse(
  tb_qlf_SP_SPIP_vol$logFC< 0 & tb_qlf_SP_SPIP_vol$FDR < 0.001 , 'darkgreen',
  ifelse(tb_qlf_SP_SPIP_vol$logFC > 0 & tb_qlf_SP_SPIP_vol$FDR < 0.001, 'dodgerblue3',
         'darkgreen'))

keyvals.color3[is.na(keyvals.color3)] <- 'grey'
names(keyvals.color3)[keyvals.color3 == 'grey'] <- 'FDR > 0.001'
names(keyvals.color3)[keyvals.color3== "darkgreen"] <- "Upregulated SPIP"
names(keyvals.color3)[keyvals.color3== "dodgerblue3"] <- "Upregulated SP"


SP_SPIP_vol_2 <-EnhancedVolcano(tb_qlf_SP_SPIP_vol,
                                lab = tb_qlf_SP_SPIP_vol$SYMBOL,
                                selectLab = c("Dio3", "Dio2", "Gpr50"),
                                labSize = 3,
                                x = 'logFC',
                                y = 'FDR',
                                xlab = bquote(~Log[2]~ 'fold change'),
                                ylab= bquote(~-Log[10]~ 'FDR'),
                                title = 'SP vs SPIP',
                                pCutoff = 0.001,
                                FCcutoff = 0,
                                pointSize = 2,
                                colCustom = keyvals.color3,
                                colAlpha = 0.9,
                                legendPosition = 'top',
                                legendLabSize = 15,
                                legendIconSize = 5.0,
                                drawConnectors = TRUE,
                                widthConnectors = 0.75,
                                colConnectors = 'black',
                                arrowheads = FALSE,
                                gridlines.major = TRUE,
                                gridlines.minor = FALSE,
                                border = 'partial',
                                borderWidth = 1.5,
                                borderColour = 'black',
                                caption = bquote(~Log[2]~ "fold change cutoff, 0; FDR cutoff, 0.001; total = 13345 variables"),
                                captionLabSize = 10,
                                xlim = c(-10,10),
                                ylim =c(0,17))



ggsave("output/Volcano_SP_SPIP.pdf", plot=SP_SPIP_vol_2, dpi = 300, height = 7, width = 7)


# 4A - Volcano plot LPIP vs LP

keyvals.color <- ifelse(
  tb_qlf_LP_LPIP_vol$logFC< 0 & tb_qlf_LP_LPIP_vol$FDR < 0.001 , 'dodgerblue3',
  ifelse(tb_qlf_LP_LPIP_vol$logFC > 0 & tb_qlf_LP_LPIP_vol$FDR < 0.001, 'tomato3',
         'grey'))

keyvals.color[is.na(keyvals.color)] <- 'grey'
names(keyvals.color)[keyvals.color == 'grey'] <- 'FDR > 0.001'
names(keyvals.color)[keyvals.color== "dodgerblue3"] <- "Upregulated LPIP"
names(keyvals.color)[keyvals.color== "tomato3"] <- "Upregulated LP"



# LP vs LPIP 

LP_LPIP_vol <-EnhancedVolcano(tb_qlf_LP_LPIP_vol,
                              lab = tb_qlf_LP_LPIP_vol$SYMBOL,
                              selectLab = c("Dio3", "Dio2", "Gpr50"),
                              labSize = 3,
                              x = 'logFC',
                              y = 'FDR',
                              xlab = bquote(~Log[2]~ 'fold change'),
                              ylab= bquote(~-Log[10]~ 'FDR'),
                              title = 'LP vs LPIP',
                              pCutoff = 0.001,
                              FCcutoff = 0,
                              pointSize = 2,
                              colCustom = keyvals.color,
                              colAlpha = 0.9,
                              legendPosition = 'top',
                              legendLabSize = 15,
                              legendIconSize = 5.0,
                              drawConnectors = TRUE,
                              widthConnectors = 0.75,
                              colConnectors = 'black',
                              arrowheads = FALSE,
                              gridlines.major = TRUE,
                              gridlines.minor = FALSE,
                              border = 'partial',
                              borderWidth = 1.5,
                              borderColour = 'black',
                              caption = bquote(~Log[2]~ "fold change cutoff, 0; FDR cutoff, 0.001; total = 13345 variables"),
                              captionLabSize = 10,
                              xlim = c(-10,10),
                              ylim =c(0,17))






ggsave("output/Volcano_LP_LPIP.pdf", plot= LP_LPIP_vol, dpi = 300, height = 7, width = 7)


# 4C - done in Graphpad prism 8

# 4D - Venndiagram done in inkscape

# 4E - heatmap SPIP uDEG

SPIP_transition <- read_delim("109-SPIP-transistion-genes-for-heatmap.txt")

SPIP_transition <- SPIP_transition %>% 
  dplyr::select(ENSEMBLID)

df_SPIP_trans <- tb_cpm_anno %>%
  right_join(SPIP_transition)


mt_SPIP_trans <- df_SPIP_trans %>% 
  column_to_rownames(var="SYMBOL") %>% 
  dplyr::select(4:23)



mt_SPIP_trans_rowscale <- t( scale(t(mt_SPIP_trans), center = TRUE, scale = TRUE) )

hm_transition_symb <- Heatmap( mt_SPIP_trans_rowscale,
                               name = "z-score",
                               col = viridis::viridis(100),            
                               show_row_names = TRUE,
                               cluster_columns = FALSE,
                               cluster_rows = TRUE,
                               show_row_dend = FALSE,
                               row_names_max_width = unit(6, "cm"),
                               row_names_gp = gpar(fontsize = 10),
                               column_names_side = c("bottom"),
                               column_names_rot = 45,
                               column_names_centered = TRUE, 
                               row_title_rot = 0)



pdf(file="output/heatmap_SPIP_transition_genes_symbols_ind_cpm_221013.pdf", width = 14, height=21)
draw( hm_transition_symb )
dev.off()

# 4F - enriched TFBS
TFBS_odds <- read.csv("TFBS-odds-all.csv", header =TRUE)
TFBS_odds %>% sample_n(5)
TFBS_odds %>%
  ggplot(aes(x= odds.ratio, y= reorder(motif.name, odds.ratio), fill = log10fdr))+
  geom_col()+
  scale_fill_gradientn(colours = viridis::viridis(20), limits = c(1.36,2.5), oob = scales::squish, name = '-Log10FDR')+
  scale_x_continuous(limits = c(-15, 15))+
  cowplot::theme_cowplot() +
  theme_bw()+
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) 

