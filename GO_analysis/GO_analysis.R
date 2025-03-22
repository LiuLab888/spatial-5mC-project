#### R packages ####
library(DOSE)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggnewscale)
library(ggupset)
library(europepmc)
library(enrichplot)
library(biomaRt)
library(RColorBrewer)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ReactomePA)
library(openxlsx)

####  Nutrient-supplier Stem Cell ####

library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

rm(list=ls())

#### raw data ####

region1 <- read.table( 
  ".../data/mE65spM_4_2_split/Cluster3/mE65spM_4_2_split_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


region2 <- read.table( 
  ".../data/mE75spM_3_4_split/Cluster1/mE75spM_3_4_split_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


region3 <- read.table( 
  ".../data/mE85spM_1_4_split/Cluster6/mE85spM_1_4_split_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


ref_region <- read.table(
  ".../data/mm39_10kb_bin_order_tab.txt",
  sep = "\t", header = F)


#### Peak extracting #### 

region_merge1 <- merge(ref_region, region1, by.x="V4", by.y="V5", sort=F)

colnames(region_merge1) <- c("bin","chr","start","end","level","CpG","C","T")

peak1 <- data.frame(chr=region_merge1$chr,
                    start=region_merge1$start,
                    end=region_merge1$end)

write.table(peak1, ".../data/peak/mE65spM_4_2_split_peak.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")



region_merge2 <- merge(ref_region, region2, by.x="V4", by.y="V5", sort=F)

colnames(region_merge2) <- c("bin","chr","start","end","level","CpG","C","T")

peak2 <- data.frame(chr=region_merge2$chr,
                    start=region_merge2$start,
                    end=region_merge2$end)

write.table(peak2, ".../data/peak/mE75spM_3_4_split_peak.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


  
region_merge3 <- merge(ref_region, region3, by.x="V4", by.y="V5", sort=F)

colnames(region_merge3) <- c("bin","chr","start","end","level","CpG","C","T")

peak3 <- data.frame(chr=region_merge3$chr,
                    start=region_merge3$start,
                    end=region_merge3$end)

write.table(peak3, ".../data/peak/mE85spM_1_4_split_peak.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


#### annotation ####

file_path1 <- ".../data/peak/mE65spM_4_2_split_peak.txt"

file_path2 <- ".../data/peak/mE75spM_3_4_split_peak.txt"

file_path3 <- ".../data/peak/mE85spM_1_4_split_peak.txt"


peakAnno1 <- annotatePeak(file_path1, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno1)



peakAnno2 <- annotatePeak(file_path2, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno2)



peakAnno3 <- annotatePeak(file_path3, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

plotAnnoPie(peakAnno3)



peakAnnoList <- lapply(list("mE6.5"=file_path1, "mE7.5"=file_path2, "mE8.5"=file_path3), 
                       annotatePeak, 
                       TxDb =txdb,
                       tssRegion=c(-2000, 500))

plotAnnoBar(peakAnnoList)


setwd(".../data/peak")


png(file="peak_distribution.png",width=500,height=300,res=100)

plotAnnoBar(peakAnnoList)

dev.off()


pdf(file="peak_distribution.pdf", width=5, height=3)

plotAnnoBar(peakAnnoList)

dev.off()


#### 300bp ####

library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#### raw data ####

region1 <- read.table( 
  ".../data/mE65spM_4_2_split/Cluster3/mE65spM_4_2_split_300bp_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


region2 <- read.table( 
  ".../data/mE75spM_3_4_split/Cluster1/mE75spM_3_4_split_300bp_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


region3 <- read.table( 
  ".../data/mE85spM_1_4_split/Cluster6/mE85spM_1_4_split_300bp_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


ref_region <- read.table(
  ".../data/mm39_300bp_bin_order_tab.txt",
  sep = "\t", header = F)


#### Peak #### 

region_merge1 <- merge(ref_region, region1, by.x="V4", by.y="V5", sort=F)

colnames(region_merge1) <- c("bin","chr","start","end","level","CpG","C","T")

peak1 <- data.frame(chr=region_merge1$chr,
                    start=region_merge1$start,
                    end=region_merge1$end)

write.table(peak1, ".../data/peak/mE65spM_4_2_split_peak_300bp.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



region_merge2 <- merge(ref_region, region2, by.x="V4", by.y="V5", sort=F)

colnames(region_merge2) <- c("bin","chr","start","end","level","CpG","C","T")

peak2 <- data.frame(chr=region_merge2$chr,
                    start=region_merge2$start,
                    end=region_merge2$end)

write.table(peak2, ".../data/peak/mE75spM_3_4_split_peak_300bp.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



region_merge3 <- merge(ref_region, region3, by.x="V4", by.y="V5", sort=F)

colnames(region_merge3) <- c("bin","chr","start","end","level","CpG","C","T")

peak3 <- data.frame(chr=region_merge3$chr,
                    start=region_merge3$start,
                    end=region_merge3$end)

write.table(peak3, ".../data/peak/mE85spM_1_4_split_peak_300bp.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")


#### annotation ####

file_path1 <- ".../data/peak/mE65spM_4_2_split_peak_300bp.txt"

file_path2 <- ".../data/peak/mE75spM_3_4_split_peak_300bp.txt"

file_path3 <- ".../data/peak/mE85spM_1_4_split_peak_300bp.txt"


peakAnno1 <- annotatePeak(file_path1, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno1)



peakAnno2 <- annotatePeak(file_path2, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno2)



peakAnno3 <- annotatePeak(file_path3, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

plotAnnoPie(peakAnno3)



peakAnnoList <- lapply(list("mE6.5"=file_path1, "mE7.5"=file_path2, "mE8.5"=file_path3), 
                       annotatePeak, 
                       TxDb =txdb,
                       tssRegion=c(-2000, 500))

plotAnnoBar(peakAnnoList)


setwd(".../data/peak")


png(file="peak_distribution_300bp.png",width=500,height=300,res=100)

plotAnnoBar(peakAnnoList)

dev.off()


pdf(file="peak_distribution_300bp.pdf", width=5, height=3)

plotAnnoBar(peakAnnoList)

dev.off()



#### overlapped peaks -- GO ####

bins_merge <-  region1$V5[which(region1$V5 %in% region2$V5)]

bins_merge1 <- bins_merge[which(bins_merge %in% region3$V5)]

ref_region_filter <- ref_region[which(ref_region$V4 %in% bins_merge1),]

write.table(ref_region_filter[,1:3], ".../data/peak/overlapped_peaks_300bp.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



file_path <- ".../data/peak/overlapped_peaks_300bp.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

#gene_list <- peakAnno@anno$geneId[which(abs(as.numeric(peakAnno@anno$distanceToTSS)) <= 10000)]

gene_list <- peakAnno@anno$geneId[which(peakAnno@anno$annotation == "Promoter (<=1kb)") ]


setwd(".../data/peak")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_hypo_methyl_merge_peaks_300bp_promoter.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]


setwd(".../data/peak")


geneGO_enrich_All_results_sorted <- read.xlsx("GO_hypo_methyl_merge_peaks_300bp_promoter.xlsx", sheet=1)

Description <- c("positive regulation of cell cycle", "actin binding", "regulation of mRNA processing" ,
                 "cell leading edge", "nucleocytoplasmic transport",
                 "regulation of organelle assembly", "autophagy" ,"lamellipodium", 
                 "regulation of supramolecular fiber organization", "phospholipid binding" )

geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[
  which(geneGO_enrich_All_results_sorted$Description %in% Description), ]



geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdYlBu",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot




png(file="GO_hypo_methyl_merge_peaks_300bp_promoter_used.png", width=800, height=400)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_merge_peaks_300bp_promoter_used.pdf", width=8, height=4)
geneGO_enrich_All_results_subset_plot

dev.off()



png(file="GO_hypo_methyl_merge_peaks_300bp_promoter.png", width=800, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_merge_peaks_300bp_promoter.pdf", width=8, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()


####************************************************************************####
####************************************************************************####
#### Nutrient-supplier Stem Cell vs inner EPC ####
library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

rm(list=ls())

#### raw data ####

region1 <- read.table( 
  ".../data/mE85spM_1_4_split/Cluster1_EM/mE85spM_1_4_split_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


region2 <- read.table( 
  ".../data/mE85spM_1_4_split/Cluster6/mE85spM_1_4_split_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


ref_region <- read.table(
  ".../data/mm39_10kb_bin_order_tab.txt",
  sep = "\t", header = F)


#### Peak extracting #### 

region_merge1 <- merge(ref_region, region1, by.x="V4", by.y="V5", sort=F)

colnames(region_merge1) <- c("bin","chr","start","end","level","CpG","C","T")

peak1 <- data.frame(chr=region_merge1$chr,
                    start=region_merge1$start,
                    end=region_merge1$end)

write.table(peak1, ".../data/peak_PD/Cluster1_EM_peak.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")



region_merge2 <- merge(ref_region, region2, by.x="V4", by.y="V5", sort=F)

colnames(region_merge2) <- c("bin","chr","start","end","level","CpG","C","T")

peak2 <- data.frame(chr=region_merge2$chr,
                    start=region_merge2$start,
                    end=region_merge2$end)

write.table(peak2, ".../data/peak_PD/Cluster6_peak.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


#### annotation ####

file_path1 <- ".../data/peak_PD/Cluster1_EM_peak.txt"

file_path2 <- ".../data/peak_PD/Cluster6_peak.txt"



peakAnno1 <- annotatePeak(file_path1, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno1)



peakAnno2 <- annotatePeak(file_path2, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno2)



peakAnnoList <- lapply(list("Placenta"=file_path1, "Decidua"=file_path2), 
                       annotatePeak, 
                       TxDb =txdb,
                       tssRegion=c(-2000, 500))

plotAnnoBar(peakAnnoList)


setwd(".../data/peak_PD")


png(file="peak_distribution.png",width=500,height=300,res=100)

plotAnnoBar(peakAnnoList)

dev.off()


pdf(file="peak_distribution.pdf", width=5, height=3)

plotAnnoBar(peakAnnoList)

dev.off()



#### overlapped peaks -- GO  ####

bins_merge <-  region1$V5[which(region1$V5 %in% region2$V5)]

ref_region_filter <- ref_region[which(ref_region$V4 %in% bins_merge),]

write.table(ref_region_filter[,1:3], ".../data/peak_PD/overlapped_peaks.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")



file_path <- ".../data/peak_PD/overlapped_peaks.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-2000, 500), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId



setwd(".../data/peak_PD")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_hypo_methyl_merge_peaks.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]




#先提取富集结果表
geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[1:20,]


#指定绘图顺序（转换为因子）：
geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

#Top20富集条形图：
mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdBu",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_hypo_methyl_merge_peaks.png", width=800, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_merge_peaks.pdf", width=8, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()


#### 300bp ####
library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rm(list=ls())

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

#### raw data ####

region1 <- read.table( 
  ".../data/mE85spM_1_4_split/Cluster1_EM/mE85spM_1_4_split_300bp_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


region2 <- read.table( 
  ".../data/mE85spM_1_4_split/Cluster6/mE85spM_1_4_split_300bp_bin_methylation_level_new_filter2.txt",
  sep = "\t", header = F)


ref_region <- read.table(
  ".../data/mm39_300bp_bin_order_tab.txt",
  sep = "\t", header = F)


#### Peak extracting #### 

region_merge1 <- merge(ref_region, region1, by.x="V4", by.y="V5", sort=F)

colnames(region_merge1) <- c("bin","chr","start","end","level","CpG","C","T")

peak1 <- data.frame(chr=region_merge1$chr,
                    start=region_merge1$start,
                    end=region_merge1$end)

write.table(peak1, ".../data/peak_PD/Cluster1_EM_peak_300bp.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



region_merge2 <- merge(ref_region, region2, by.x="V4", by.y="V5", sort=F)

colnames(region_merge2) <- c("bin","chr","start","end","level","CpG","C","T")

peak2 <- data.frame(chr=region_merge2$chr,
                    start=region_merge2$start,
                    end=region_merge2$end)

write.table(peak2, ".../data/peak_PD/Cluster6_peak_300bp.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")


#### annotation ####

file_path1 <- ".../data/peak_PD/Cluster1_EM_peak_300bp.txt"

file_path2 <- ".../data/peak_PD/Cluster6_peak_300bp.txt"



peakAnno1 <- annotatePeak(file_path1, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno1)



peakAnno2 <- annotatePeak(file_path2, tssRegion=c(-2000, 500), TxDb=txdb, 
                          annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)
plotAnnoPie(peakAnno2)



peakAnnoList <- lapply(list("Placenta"=file_path1, "Decidua"=file_path2), 
                       annotatePeak, 
                       TxDb =txdb,
                       tssRegion=c(-2000, 500))

plotAnnoBar(peakAnnoList)


setwd(".../data/peak_PD")


png(file="peak_distribution_300bp.png",width=500,height=300,res=100)

plotAnnoBar(peakAnnoList)

dev.off()


pdf(file="peak_distribution_300bp.pdf", width=5, height=3)

plotAnnoBar(peakAnnoList)

dev.off()


#### overlapped peaks ####

#### promoter ####

bins_merge <-  region1$V5[which(region1$V5 %in% region2$V5)]

ref_region_filter <- ref_region[which(ref_region$V4 %in% bins_merge),]

write.table(ref_region_filter[,1:3], ".../data/peak_PD/overlapped_peaks_300bp.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")



file_path <- ".../data/peak_PD/overlapped_peaks_300bp.txt"


peakAnno <- annotatePeak(file_path, tssRegion=c(-2000, 500), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which(peakAnno@anno$annotation == "Promoter (<=1kb)" & abs(peakAnno@anno$distanceToTSS) == 0 ) ]



#### Go ####

setwd(".../data/peak_PD")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_hypo_methyl_merge_peaks_300bp.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]


setwd(".../data/peak_PD")

geneGO_enrich_All_results_sorted <- read.xlsx("GO_hypo_methyl_merge_peaks_300bp.xlsx", sheet=1)

Description <- c("regulation of cellular component size" ,"regulation of epithelial cell proliferation" ,
                 "regulation of actin filament organization" , "regulation of epithelial cell migration",
                 "regulation of endothelial cell migration" , "positive regulation of cell cycle",
                 "regulation of protein kinase B signaling", "protein serine/threonine kinase activity" )

geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[
  which(geneGO_enrich_All_results_subset$Description %in% Description),]


geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "BrBG",direction = -1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_hypo_methyl_merge_peaks_300bp_used.png", width=900, height=400)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_merge_peaks_300bp_used.pdf", width=9, height=4)

geneGO_enrich_All_results_subset_plot

dev.off()



png(file="GO_hypo_methyl_merge_peaks_300bp.png", width=900, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_merge_peaks_300bp.pdf", width=9, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()


####************************************************************************####
####************************************************************************####

#### Nutrient-supplier Stem Cell only -- GO ####

#### promoter ####

bins_merge <-  region2$V5[-which(region2$V5 %in% region1$V5)]

ref_region_filter <- ref_region[which(ref_region$V4 %in% bins_merge),]

write.table(ref_region_filter[,1:3], ".../data/peak_PD/Decidua_specific_peaks_300bp.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


file_path <- ".../data/peak_PD/Decidua_specific_peaks_300bp.txt"

peakAnno <- annotatePeak(file_path, tssRegion=c(-2000, 500), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which( peakAnno@anno$annotation == "Promoter (<=1kb)" & abs(peakAnno@anno$distanceToTSS) == 0 )  ]


#### Go ####

setwd(".../data/peak_PD")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_hypo_methyl_Decidua_specific_peaks_300bp.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]



setwd(".../data/peak_PD")

geneGO_enrich_All_results_sorted <- read.xlsx("GO_hypo_methyl_Decidua_specific_peaks_300bp.xlsx")

Description <- c("axonogenesis", "regulation of actin filament-based process",
                 "regulation of mRNA metabolic process", "regulation of protein catabolic process", 
                 "exocytosis",  "columnar/cuboidal epithelial cell differentiation")

geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[
  which(geneGO_enrich_All_results_sorted$Description %in% Description),]


geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "OrRd",direction = 1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot



png(file="GO_hypo_methyl_Decidua_specific_peaks_300bp_used.png", width=800, height=300)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_Decidua_specific_peaks_300bp_used.pdf", width=8, height=3)
geneGO_enrich_All_results_subset_plot

dev.off()



png(file="GO_hypo_methyl_Decidua_specific_peaks_300bp.png", width=800, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_Decidua_specific_peaks_300bp.pdf", width=8, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()


####************************************************************************####
####************************************************************************####

#### inner EPC only -- GO ####
#### promoter ####


bins_merge <-  region1$V5[-which(region1$V5 %in% region2$V5)]

ref_region_filter <- ref_region[which(ref_region$V4 %in% bins_merge),]

write.table(ref_region_filter[,1:3], ".../data/peak_PD/Placenta_specific_peaks_300bp.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


file_path <- ".../data/peak_PD/Placenta_specific_peaks_300bp.txt"

peakAnno <- annotatePeak(file_path, tssRegion=c(-2000, 500), TxDb=txdb, 
                         annoDb = 'org.Mm.eg.db', addFlankGeneInfo=TRUE, flankDistance=5000)

gene_list <- peakAnno@anno$geneId[which( peakAnno@anno$annotation == "Promoter (<=1kb)" & abs(peakAnno@anno$distanceToTSS) == 0 )]


#### GO ####

setwd(".../data/peak_PD")


geneGO_enrich_All_results <- enrichGO(gene = gene_list, minGSSize = 5, 
                                      OrgDb = org.Mm.eg.db, ont = "ALL", pAdjustMethod = "BH",
                                      keyType = 'ENTREZID', pvalueCutoff = 0.01, qvalueCutoff = 0.2, 
                                      readable = TRUE)

View(geneGO_enrich_All_results@result)


geneGO_enrich_All_results_sorted <- geneGO_enrich_All_results@result[order(geneGO_enrich_All_results@result$pvalue),]

write.xlsx(geneGO_enrich_All_results_sorted, "GO_hypo_methyl_Placenta_specific_peaks_300bp.xlsx")


geneGO_enrich_All_results_sorted$Description[1:20]



setwd(".../data/peak_PD")

geneGO_enrich_All_results_sorted <- read.xlsx("GO_hypo_methyl_Placenta_specific_peaks_300bp.xlsx",sheet=1)

Description <- c("GTPase regulator activity",  "positive regulation of GTPase activity",
                 "transmembrane receptor protein serine/threonine kinase signaling pathway","connective tissue development",
                 "cell-substrate junction", "ATP hydrolysis activity")

geneGO_enrich_All_results_subset <- geneGO_enrich_All_results_sorted[
  which(geneGO_enrich_All_results_subset$Description %in% Description),]

geneGO_enrich_All_results_subset$Description[
  which(geneGO_enrich_All_results_subset$Description == 
          "transmembrane receptor protein serine/threonine kinase signaling pathway")] <- "TRP serine/threonine kinase signaling pathway"

geneGO_enrich_All_results_subset$terms <- factor(geneGO_enrich_All_results_subset$Description,
                                                 levels = rev(geneGO_enrich_All_results_subset$Description))

mytheme <-  theme(legend.title = element_text(size=14, family="sans", face="plain"))+
  theme(legend.text = element_text(size=14, family="sans", face="plain"))+
  theme(plot.title = element_text(size = 14,family="sans", face="plain", hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14, family="sans", face="plain")) +
  theme(axis.title.y = element_text(size = 14, family="sans", face="plain"))+
  theme(axis.text.x = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(axis.text.y = element_text(size = 14, color="black", face="plain",
                                   family="sans", colour="black")) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, 
                                    linetype = "solid")) 

geneGO_enrich_All_results_subset_plot <- ggplot(data = geneGO_enrich_All_results_subset,
                                                aes(x = Count, y = terms, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "OrRd",direction = 1) +
  geom_bar(stat = "identity", width = 0.75) +
  theme_bw() +
  labs(x = "Number of Gene",
       y = "Gene terms",
       title = NULL) + mytheme

geneGO_enrich_All_results_subset_plot




png(file="GO_hypo_methyl_Placenta_specific_peaks_300bp_used.png", width=800, height=300)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_Placenta_specific_peaks_300bp_used.pdf", width=8, height=3)
geneGO_enrich_All_results_subset_plot

dev.off()



png(file="GO_hypo_methyl_Placenta_specific_peaks_300bp.png", width=1000, height=600)

geneGO_enrich_All_results_subset_plot

dev.off()


pdf(file="GO_hypo_methyl_Placenta_specific_peaks_300bp.pdf", width=10, height=6)
geneGO_enrich_All_results_subset_plot

dev.off()



