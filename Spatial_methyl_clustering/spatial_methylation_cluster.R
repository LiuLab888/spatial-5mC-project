####*************************************************************************####
####*************************************************************************####
#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)

#### mE55spM_A_7_split_merge only embryo ####

rm(list=ls())

setwd(".../data_analysis/mE55spM_A_7_split_merge")

#### raw data ####

input_name1 <- ".../data/mE55spM_A_7_split_merge/coverage_matrix_total.txt"

input_name2 <- ".../data/mE55spM_A_7_split_merge/methylation_matrix_total.txt"

refer_path <- ".../data/mE55spM_A_7_split_merge/mE55spM_A_7_split_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 5000 & 
                                               merge_barcode1$methylation >= 0.32),]


#### border filtering ####


index_filtered_EM <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE55spM_A_7_split_merge/cell_filtered_EM.txt", header=T)

merge_barcode_filter1 <- merge_barcode_filter[which(merge_barcode_filter$coord %in% index_filtered_EM$coord ),]



dim2 <-  merge_barcode_filter1$V2

dim1 <-  merge_barcode_filter1$V3


value1 <- scale(sqrt((dim1-1)^2 + (dim2-20)^2 ))

value2 <- scale(sqrt((dim1-15)^2 + (dim2-20)^2 ))

value3 <- scale(sqrt((dim1-1)^2 + (dim2-1)^2 ))

value4 <- scale(sqrt((dim1-15)^2 + (dim2-1)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=methylation)

PCA_selection <- data_merge


#### Kmeans ####

set.seed(100)

km_results <- kmeans(PCA_selection, 1, nstart = 25)

#fviz_cluster(km_results, PCA_selection) #显示聚类分布情况

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))


write.table(data_cluster,"data_cluster_EM_part.txt",
            row.names = FALSE,
            col.names = TRUE)


#### umap ####

config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.8 

config.params$n_neighbors=10

config.params$n_components=2

#config.params$set_op_mix_ratio=0.5

#config.params$local_connectivity=2

umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]

data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all_EM_part.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]


#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")     
                   
plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter

plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-6, 6.2),
                     breaks = c( -3, 0, 3, 6),
                     labels = c("-3", "0", "3", "6")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-6, 6.1),
                     breaks = c(-3, 0, 3, 6),
                     labels = c("-3", "0", "3", "6")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw_EM_part.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw_EM_part.pdf", width=8, height=8)

plot_scatter

dev.off()


#### plot location #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")   
                   

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=3,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,20,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,15,5)) +
  expand_limits(y=c(1,20),x=c(2,12))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location_raw_EM_part_new.png",width=150/1.6,height=250/1.6)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw_EM_part_new.pdf",width=1.5/1.35,height=2.5/1.35)

spatialPlot 

dev.off()


####*************************************************************************####
####*************************************************************************####
#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)


#### mE65spM_4_6_split_merge only embryo ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE65spM_4_6_split_merge")

#### raw data ####

input_name1 <- ".../data/mE65spM_4_6_split_merge/coverage_matrix_total.txt"

input_name2 <- ".../data/mE65spM_4_6_split_merge/methylation_matrix_total.txt"

refer_path <- ".../data/mE65spM_4_6_split_merge/mE65spM_4_6_split_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 1000 & 
                                               merge_barcode1$methylation >= 0.2),]

index_filtered_EM <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE65spM_4_6_split_merge/cell_filtered_EM.txt", header=T)

merge_barcode_filter <- merge_barcode_filter[which(merge_barcode_filter$coord %in% index_filtered_EM$coord ),]



#### border filtering ####

if(FALSE){
  
  source("/media/tangym/ETD/Spatial_Methylation/data_analysis/0_rm_outliner.shaorui.R")
  
  coord <- merge_barcode_filter$coord
  
  dat <- coord_to_dat(coord)
  
  value <- merge_barcode_filter$methylation
  
  dat$value <- value
  
  tissue_selector(dat)
  

  index_filtered <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE65spM_4_6_split_merge/cell_filtered.txt", header=T, sep = "\t")
  
  merge_barcode_filter1 <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered$coord),]
  
}


merge_barcode_filter1 <- merge_barcode_filter


dim1 <-  merge_barcode_filter1$V2

dim2 <-  merge_barcode_filter1$V3


value1 <- scale(sqrt((dim1-42)^2 + (dim2-1)^2 ))

value2 <- scale(sqrt((dim1-42)^2 + (dim2-20)^2 ))

value3 <- scale(sqrt((dim1-42)^2 + (dim2-39)^2 ))


value4 <- scale(sqrt((dim1-42)^2 + (dim2-58)^2 ))

value5 <- scale(sqrt((dim1-42)^2 + (dim2-77)^2 ))

value6 <- scale(sqrt((dim1-42)^2 + (dim2-96)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=methylation)

PCA_selection <- data_merge


#### Kmeans ####

set.seed(100)

km_results <- kmeans(PCA_selection, 6, nstart = 25)

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))



write.table(data_cluster,"data_cluster_EM.txt",
            row.names = FALSE,
            col.names = TRUE)


#### umap ####

config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.5 

config.params$n_neighbors=10

config.params$n_components=2

umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]

data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all_EM.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]


#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")
                   

plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-11, 10),
                     breaks = c( -8, 0, 8),
                     labels = c("-8", "0", "8")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-9, 13),
                     breaks = c(-8, 0, 8),
                     labels = c("-8", "0", "8")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw_EM.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw_EM.pdf", width=8, height=8)

plot_scatter

dev.off()



#### plot location #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")     
                   

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=3,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(breaks = seq(0,95,5)) +
  expand_limits(y=c(36,49),x=c(7,57))+ labs(color="Cluster") + theme_void() + coord_flip()

spatialPlot



png(file="Cluster_location_raw_EM_part_new.png",width=150/1.4,height=650/1.4)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw_EM_part_new.pdf",width=1.5/1.15,height=6.5/1.15)

spatialPlot 

dev.off()



####*************************************************************************####
####*************************************************************************####
#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)


#### mE75spM_6_6_split_merge only embryo ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE75spM_6_6_split_merge")

#### raw data ####

input_name1 <- ".../data/mE75spM_6_6_split_merge/coverage_matrix_total.txt"

input_name2 <- ".../data/mE75spM_6_6_split_merge/methylation_matrix_total.txt"

refer_path <- ".../data/mE75spM_6_6_split_merge/mE75spM_6_6_split_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 5000 & 
                                               merge_barcode1$methylation >= 0.32),]

index_filtered_EM <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE75spM_6_6_split_merge/mE75_6_6_cell_filter_4th.txt", header=T)

merge_barcode_filter <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered_EM$coord ),]



#### border filtering ####

if(FALSE){
  
  source("/media/tangym/ETD/Spatial_Methylation/data_analysis/0_rm_outliner.shaorui.new.R")
  
  coord <- merge_barcode_filter$coord
  
  dat <- coord_to_dat(coord)
  
  value <- merge_barcode_filter$methylation
  
  dat$value <- value
  
  tissue_selector(dat, 1000, 300)
  
  
}


index_filtered <- read.table(".../data/mE75spM_6_6_split_merge/cell_filtered_EM_and_TPs_and_cavity.txt", 
                             header=T, sep = "\t")

merge_barcode_filter1 <- merge_barcode_filter[which(merge_barcode_filter$coord %in% index_filtered$coord),]


dim2 <-  merge_barcode_filter1$V2

dim1 <-  merge_barcode_filter1$V3


value1 <- scale(sqrt((dim1-14)^2 + (dim2-1)^2 ))

value2 <- scale(sqrt((dim1-14)^2 + (dim2-20)^2 ))

value3 <- scale(sqrt((dim1-14)^2 + (dim2-39)^2 ))


value4 <- scale(sqrt((dim1-14)^2 + (dim2-58)^2 ))

value5 <- scale(sqrt((dim1-14)^2 + (dim2-77)^2 ))

value6 <- scale(sqrt((dim1-14)^2 + (dim2-96)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=methylation)

PCA_selection <- data_merge


#### Kmeans ####

set.seed(100)

km_results <- kmeans(PCA_selection, 6, nstart = 25)

#fviz_cluster(km_results, PCA_selection) #显示聚类分布情况

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))


write.table(data_cluster,"data_cluster_EM_part.txt",
            row.names = FALSE,
            col.names = TRUE)



#### umap ####

config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.8 

config.params$n_neighbors=10

config.params$n_components=2

umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]

data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all_EM_part.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]


#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")     
                   
plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-12, 13),
                     breaks = c( -10, 0, 10),
                     labels = c("-10", "0", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-14, 13),
                     breaks = c(-10, 0, 10),
                     labels = c("-10", "0", "10")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw_EM_part.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw_EM_part.pdf", width=8, height=8)

plot_scatter

dev.off()


#### plot location #### 

data_merge_all <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE75spM_6_6_split_merge/data_merge_all_EM_part.txt",
                             header=T)

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")   
                   

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=3,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,24,5)) +
  expand_limits(y=c(7,96),x=c(1,25))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location_raw_EM_part.png",width=250/1.25,height=850/1.25)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw_EM_part.pdf",width=2.5,height=8.5)

spatialPlot 

dev.off()


####*************************************************************************####
####*************************************************************************####
#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)


#### mE85spM_1_4_split_merge -- embryo and extraembryo ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE85spM_1_4_split_merge")

#### raw data ####

input_name1 <- ".../data/mE85spM_1_4_split_merge/coverage_matrix_total.txt"

input_name2 <- ".../data/mE85spM_1_4_split_merge/methylation_matrix_total.txt"

refer_path <- ".../data/mE85spM_1_4_split_merge/mE85spM_1_4_split_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 2000 & 
                                               merge_barcode1$methylation >= 0.2),]

index_filtered_EM <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE85spM_1_4_split_merge/cell_filtered_EM.txt", header=T)

merge_barcode_filter <- merge_barcode_filter[which(merge_barcode_filter$coord %in% index_filtered_EM$coord ),]



#### border filtering ####

if(FALSE){

source("/media/tangym/ETD/Spatial_Methylation/data_analysis/0_rm_outliner.shaorui.R")

coord <- merge_barcode_filter$coord

dat <- coord_to_dat(coord)

value <- merge_barcode_filter$methylation

dat$value <- value

tissue_selector(dat)


index_filtered <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE85spM_1_4_split_merge/cell_filtered.txt", header=T, sep = "\t")

merge_barcode_filter1 <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered$coord),]

}


merge_barcode_filter1 <- merge_barcode_filter


dim1 <-  merge_barcode_filter1$V2

dim2 <-  merge_barcode_filter1$V3


value1 <- scale(sqrt((dim1-48)^2 + (dim2-1)^2 ))

value2 <- scale(sqrt((dim1-48)^2 + (dim2-20)^2 ))

value3 <- scale(sqrt((dim1-48)^2 + (dim2-39)^2 ))


value4 <- scale(sqrt((dim1-48)^2 + (dim2-58)^2 ))

value5 <- scale(sqrt((dim1-48)^2 + (dim2-77)^2 ))

value6 <- scale(sqrt((dim1-48)^2 + (dim2-96)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=methylation)

PCA_selection <- data_merge


#### Kmeans ####

set.seed(100)

km_results <- kmeans(PCA_selection, 6, nstart = 25)

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))

# 1   2   3   4   5   6 
# 538 454 641 538 612 370 

write.table(data_cluster,"data_cluster_EM.txt",
            row.names = FALSE,
            col.names = TRUE)



#### umap ####

config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.8 

config.params$n_neighbors=10

config.params$n_components=2

#config.params$set_op_mix_ratio=0.5

#config.params$local_connectivity=2

umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]


data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all_EM.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]


#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")
                   

plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-18, 15),
                     breaks = c( -10, 0, 10),
                     labels = c("-10", "0", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-16, 14),
                     breaks = c(-10, 0, 10),
                     labels = c("-10", "0", "10")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw_EM.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw_EM.pdf", width=8, height=8)

plot_scatter

dev.off()



#### plot location #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")     
                   

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=3,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  expand_limits(y=c(0,96),x=c(0,96))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location_raw_EM.png",width=700,height=700)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw_EM.pdf",width=8/0.85,height=8/0.85)

spatialPlot 

dev.off()




####*************************************************************************####
####*************************************************************************####
#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)


#### mE65spM_4_2_split_merge -- Decidua ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE65spM_4_2_split_merge")

#### raw data ####

input_name1 <- ".../data/mE65spM_4_2_split_merge/coverage_matrix_total.txt"

input_name2 <- ".../data/mE65spM_4_2_split_merge/methylation_matrix_total.txt"

refer_path <- ".../data/mE65spM_4_2_split_merge/mE65spM_4_2_split_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 1000 & 
                                               merge_barcode1$methylation >= 0.32),]

#### border filtering ####

if(FALSE){
  
  source("/media/tangym/ETD/Spatial_Methylation/data_analysis/0_rm_outliner.shaorui.R")
  
  coord <- merge_barcode_filter$coord
  
  dat <- coord_to_dat(coord)
  
  value <- merge_barcode_filter$methylation
  
  dat$value <- value
  
  tissue_selector(dat)
  
}


index_filtered <- read.table("cell_filtered_EM.txt", header=T, sep = "\t")

merge_barcode_filter1 <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered$coord),]


dim1 <-  merge_barcode_filter1$V2

dim2 <-  merge_barcode_filter1$V3


value1 <- scale(sqrt((dim1-37)^2 + (dim2-7)^2 ))

value2 <- scale(sqrt((dim1-37)^2 + (dim2-23)^2 ))

value3 <- scale(sqrt((dim1-37)^2 + (dim2-38)^2 ))


value4 <- scale(sqrt((dim1-37)^2 + (dim2-54)^2 ))

value5 <- scale(sqrt((dim1-37)^2 + (dim2-69)^2 ))

value6 <- scale(sqrt((dim1-37)^2 + (dim2-85)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=methylation)

PCA_selection <- data_merge


#### Kmeans ####

set.seed(100)

km_results <- kmeans(PCA_selection, 5, nstart = 25)


PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))



write.table(data_cluster,"data_cluster.txt",
            row.names = FALSE,
            col.names = TRUE)



#### umap ####


config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.9 

config.params$n_neighbors=12

config.params$n_components=2

umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]


data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]


#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")


plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-13, 12),
                     breaks = c( -10, 0, 10),
                     labels = c("-10", "0", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-15, 17),
                     breaks = c(-10, 0, 10),
                     labels = c("-10", "0", "10")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw.pdf", width=8, height=8)

plot_scatter

dev.off()



#### plot location #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")    


spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=2,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  expand_limits(y=c(0,96),x=c(0,96))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location_raw.png",width=500,height=500)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw.pdf",width=8,height=8)

spatialPlot 

dev.off()



####*************************************************************************####
####*************************************************************************####
#### R packages ####

library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)


#### mE75spM_3_4_split -- Decidua ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE75spM_3_4_split")

#### raw data ####

input_name1 <- ".../data/mE75spM_3_4_split/coverage_matrix_total.txt"

input_name2 <- ".../data/mE75spM_3_4_split/methylation_matrix_total.txt"

refer_path <- ".../data/combine_barcode.round2round1_index1_index2.v3_big_methylation.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode2 <- read.table(input_name2,stringsAsFactors = F)


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 1000 & 
                                               merge_barcode1$methylation >= 0.32),]

index_filtered_EM <- read.table("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE75spM_3_4_split/cell_filtered_EM.txt", header=T)

merge_barcode_filter <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered_EM$coord ),]


#### border filtering ####

source("/media/tangym/ETD/Spatial_Methylation/data_analysis/0_rm_outliner.shaorui.R")

coord <- merge_barcode_filter$coord

dat <- coord_to_dat(coord)

value <- merge_barcode_filter$methylation

dat$value <- value

tissue_selector(dat)


index_filtered <- read.table("cell_filtered.txt", header=T, sep = "\t")

merge_barcode_filter1 <- merge_barcode_filter[-which(merge_barcode_filter$coord %in% index_filtered$coord),]


dim1 <-  merge_barcode_filter1$V2

dim2 <-  merge_barcode_filter1$V3



value1 <- scale(sqrt((dim1-44)^2 + (dim2-29)^2 ))

value2 <- scale(sqrt((dim1-44)^2 + (dim2-37)^2 ))

value3 <- scale(sqrt((dim1-44)^2 + (dim2-45)^2 ))


value4 <- scale(sqrt((dim1-44)^2 + (dim2-53)^2 ))

value5 <- scale(sqrt((dim1-44)^2 + (dim2-61)^2 ))

value6 <- scale(sqrt((dim1-44)^2 + (dim2-69)^2 ))


value7 <- scale(sqrt((dim1-44)^2 + (dim2-77)^2 ))

value8 <- scale(sqrt((dim1-44)^2 + (dim2-85)^2 ))

value9 <- scale(sqrt((dim1-44)^2 + (dim2-96)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/4)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=value7,
                         PCA8=value8,
                         PCA9=value9,
                         PCA10=methylation)

PCA_selection <- data_merge


#### Kmeans ####

set.seed(100)

km_results <- kmeans(PCA_selection, 8, nstart = 25)

#fviz_cluster(km_results, PCA_selection) #显示聚类分布情况

PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))


write.table(data_cluster,"data_cluster.txt",
            row.names = FALSE,
            col.names = TRUE)



#### umap ####


config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.9 

config.params$n_neighbors=10

config.params$n_components=2


umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]


data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]



#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                   "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")

                   
plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-12, 15),
                     breaks = c( -10, 0, 10),
                     labels = c("-10", "0", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-16, 18),
                     breaks = c(-10, 0, 10),
                     labels = c("-10", "0", "10")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw_new.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw_new.pdf", width=8, height=8)

plot_scatter

dev.off()



#### plot location #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
                    "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")
                    

spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=2,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  expand_limits(y=c(0,96),x=c(0,96))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location_raw_new.png",width=500,height=500)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw_new.pdf",width=8,height=8)

spatialPlot 

dev.off()

            

####*************************************************************************####
####*************************************************************************####
#### R packages ####
library(ggsci)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(factoextra)
library(umap)
library(devtools)
library(ggbiplot)
library(cluster)
library(factoextra)


#### mE85spM_1_4_split_merge -- Decidua ####

rm(list=ls())

setwd("/media/tangym/ETD/Spatial_Methylation/data_analysis/mE85spM_1_4_split_merge")

#### raw data ####

input_name1 <- ".../data/mE85spM_1_4_split_merge/coverage_matrix_total.txt"

input_name2 <- ".../data/mE85spM_1_4_split_merge/methylation_matrix_total.txt"

refer_path <- ".../data/mE85spM_1_4_split_merge/mE85spM_1_4_split_ref_barcode.txt"


# 1
ref_barcode <- read.table(refer_path,stringsAsFactors = F)

our_barcode1 <- read.table(input_name1,stringsAsFactors = F)

our_barcode1 <- our_barcode1[which(duplicated(our_barcode1$V1)==FALSE), ]


our_barcode2 <- read.table(input_name2,stringsAsFactors = F)

our_barcode2 <- our_barcode2[which(duplicated(our_barcode2$V1)==FALSE), ]


# 2
merge_barcode0 <- merge(our_barcode1, our_barcode2, by="V1", all=TRUE)

merge_barcode1 <- merge(ref_barcode, merge_barcode0, by="V1", all=TRUE)

colnames(merge_barcode1) <- c("V1","V2","V3","coverage","methylation")

merge_barcode1$coord <- paste0(merge_barcode1$V2,"x", merge_barcode1$V3)


#### filtering ####

merge_barcode_filter <- merge_barcode1[which(merge_barcode1$coverage >= 1000 & 
                                               merge_barcode1$methylation >= 0.32),]

#### border filtering ####

source("/media/tangym/ETD/Spatial_Methylation/data_analysis/0_rm_outliner.shaorui.R")

coord <- merge_barcode_filter$coord

dat <- coord_to_dat(coord)

value <- merge_barcode_filter$methylation

dat$value <- value

tissue_selector(dat)


index_filtered <- read.table("cell_filtered_EM.txt", header=T, sep = "\t")

index_filtered_other <- read.table("cell_filtered_other.txt", header=T, sep = "\t")


merge_barcode_filter1 <- merge_barcode_filter[-which(merge_barcode_filter$coord %in%
                                                       c(index_filtered$coord, index_filtered_other$coord)),]


dim1 <-  merge_barcode_filter1$V2

dim2 <-  merge_barcode_filter1$V3


value1 <- scale(sqrt((dim1-48)^2 + (dim2-1)^2 ))

value2 <- scale(sqrt((dim1-48)^2 + (dim2-20)^2 ))

value3 <- scale(sqrt((dim1-48)^2 + (dim2-39)^2 ))


value4 <- scale(sqrt((dim1-48)^2 + (dim2-58)^2 ))

value5 <- scale(sqrt((dim1-48)^2 + (dim2-77)^2 ))

value6 <- scale(sqrt((dim1-48)^2 + (dim2-96)^2 ))



methylation <- (merge_barcode_filter1$methylation-mean(merge_barcode_filter1$methylation))/
  (sd(merge_barcode_filter1$methylation)/2)


data_merge <- data.frame(PCA1=value1,
                         PCA2=value2,
                         PCA3=value3,
                         PCA4=value4,
                         PCA5=value5,
                         PCA6=value6,
                         PCA7=methylation)

PCA_selection <- data_merge


#### Kmeans####

set.seed(100)

km_results <- kmeans(PCA_selection, 6, nstart = 25)


PCA_selection_cluster <- as.data.frame(km_results$cluster)

colnames(PCA_selection_cluster)[1] <- "cluster"

data_cluster <- data.frame(V1 = merge_barcode_filter1$coord,
                           V2 = PCA_selection_cluster$cluster)

summary(factor(data_cluster$V2))


write.table(data_cluster,"data_cluster.txt",
            row.names = FALSE,
            col.names = TRUE)



#### umap ####

config.params <- umap.defaults

config.params$random_state=100

config.params$min_dist=0.9 

config.params$n_neighbors=12

config.params$n_components=2


umap_result <- umap(PCA_selection, config = config.params)

umap1 <- umap_result$layout[,1]

umap2 <- umap_result$layout[,2]

data_merge_all <- cbind(merge_barcode_filter1, cluster=data_cluster$V2, umap1, umap2)


write.table(data_merge_all,"data_merge_all.txt",
            row.names = FALSE,
            col.names = TRUE)



#### t-SNE ####

##install.packages("Rtsne")

library(Rtsne)

set.seed(100) 

tsne_result <- Rtsne(PCA_selection,
                     dims = 2, pca = F, max_iter = 1000,
                     theta = 0.4, perplexity = 10,check_duplicates = F,
                     verbose = F) 

tsne1 <- tsne_result$Y[,1]

tsne2 <- tsne_result$Y[,2]



#### plot UMAP #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

umap1 <- data_merge_all$umap1

umap2 <- data_merge_all$umap2

data_plot <- data.frame(umap1, umap2, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", "#0076B9", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")


plot_scatter <- ggplot(data=data_plot, aes(x=umap1, y=umap2, 
                                           color = data_group)) + 
  geom_point(size=2, alpha=1,shape=19 ) +
  
  scale_color_manual(values = cols) 

# stat_ellipse(aes(fill=data_group2), geom='polygon', type="norm",
#              level=0.68, alpha=0.5, show.legend = F) 

plot_scatter


plot_scatter <- plot_scatter + theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(axis.ticks = element_line(size = 0.75, color ="black"),
        axis.ticks.length=unit(1.25, 'mm'))


plot_scatter <-plot_scatter + theme(legend.title = element_text(size=18, family = "sans", face = "plain"))+
  theme(legend.text = element_text(size=16, family = "sans", face = "plain")) +
  theme(plot.title = element_text(size = 18, family = "sans", face = "plain", hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 18, family = "sans", face = "plain")) +
  theme(axis.title.y = element_text(size = 18, family = "sans", face = "plain"))+
  theme(axis.text.x = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +
  theme(axis.text.y = element_text(size = 18, color = "black", face = "plain",
                                   family = "sans", colour = "black")) +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.25, 
                                    linetype = "solid")) + #调节图形的四周边框
  labs(title = "",x = "UMAP1", y = "UMAP2") +
  theme(plot.margin = unit(c(0.5,0.5,0.2,0.2),"cm"))+
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-15, 13),
                     breaks = c( -10, 0, 10),
                     labels = c("-10", "0", "10")) +  
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-20, 19),
                     breaks = c(-15, 0, 15),
                     labels = c("-15", "0", "15")) +
  
  # theme(legend.position=c(0.97, 0.97),
  #       legend.justification = c(0.97, 0.97)) +
  # 
  # theme(legend.background = element_rect(fill="white",colour="black")) +
  
  guides(color=guide_legend(title = "Cluster"))

plot_scatter


png(file="Cluster_UMAP_raw.png", width=600, height=600)

plot_scatter

dev.off()


pdf(file="Cluster_UMAP_raw.pdf", width=8, height=8)

plot_scatter

dev.off()



#### plot location #### 

data_group <- factor(data_merge_all$cluster)

summary(data_group)

nA <- data_merge_all$V2

nB <- data_merge_all$V3

data_plot <- data.frame(nA, nB, data_group)

cols <- c("#699ECA", "#FF8C00", "#F898CB", "#4DAF4A", "#D65190", "#731A73", "#FFCB5B", "#e66f51", 
          "#3D505A", "#0098B2","#FBEA2E", "#F8B072", "#8582BD", "#4F99C9", "#A8D3A0", "#A6D0E6", "#EC3E31")    


spatialPlot <- ggplot(data=data_plot, aes(x=nB,y=nA,color=data_group))+
  geom_point(shape=16,size=2,show.legend=F)+
  scale_color_manual(values =cols) +
  scale_y_continuous(breaks = seq(0,95,5)) +
  scale_x_continuous(trans = "reverse",breaks = seq(0,95,5)) +
  expand_limits(y=c(0,96),x=c(0,96))+ labs(color="Cluster") + theme_void() 

spatialPlot


png(file="Cluster_location_raw.png",width=500,height=500)

spatialPlot 

dev.off()


pdf(file="Cluster_location_raw.pdf",width=8,height=8)

spatialPlot 

dev.off()



