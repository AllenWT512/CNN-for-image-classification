Sys.setenv(LANGUAGE="en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
library(Seurat)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(dplyr) 

A_J_S_C_W <- readRDS("/PROJ01/share/zhouwutong/AJSCW/A_J_S_C_W.rds")

A_J_S_C_W <- FindClusters(A_J_S_C_W, resolution = 3)
A_J_S_C_W <- RunUMAP(A_J_S_C_W, dims = 1:20, reduction = "integrated.Harmony")


A_J_S_C_W_markers = FindAllMarkers(A_J_S_C_W, only.pos = TRUE)
write.csv(A_J_S_C_W_markers,"/PROJ01/share/zhouwutong/AJSCW/A_J_S_C_W_markers.csv")
DimPlot(A_J_S_C_W, reduction = "umap",split.by = "orig.ident",raster=FALSE,label = T)
DimPlot(A_J_S_C_W, reduction = "umap",split.by = "seurat_clusters",raster=FALSE,label = T)
DimPlot(A_J_S_C_W, reduction = "umap",raster=FALSE,label = T)
FeaturePlot(A_J_S_C_W, features = "TLR7",raster=FALSE)

######Oligo#######
FeaturePlot(A_J_S_C_W_removed_3, features = "ST18",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "TMEM144",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "TF",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "MBP",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "MOBP",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3_removed_2, features = "PLP1",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "CNP",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "MAG",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "MYRF",raster=FALSE)


###OPC marker############################################################
FeaturePlot(A_J_S_C_W, features = "TNR",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "LHFPL3",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "MEGF11",raster=FALSE)


###Neuro progenitor marker##########cluster 28##################################################
FeaturePlot(A_J_S_C_W, features = "SPAG17",raster=FALSE)

####Tcell marker####
FeaturePlot(A_J_S_C_W, features = "SKAP1",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "TC2N",raster=FALSE)

####B cell marker####
FeaturePlot(A_J_S_C_W, features = "FCRL5",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "IGKC",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "IGHG1",raster=FALSE)

####Astrocyte marker####
FeaturePlot(A_J_S_C_W, features = "AQP4",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "ADCY2",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_2, features = "GFAP",raster=FALSE)

FeaturePlot(A_J_S_C_W, features = "ERBB2IP",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "BAI3",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "KIAA1598",raster=FALSE)

FeaturePlot(A_J_S_C_W, features = "AC025946.1",raster=FALSE)

####SC marker####
FeaturePlot(A_J_S_C_W, features = "LAMA2",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "ABCA9",raster=FALSE)
FeaturePlot(A_J_S_C_W, features = "CEMIP",raster=FALSE)


FeaturePlot(A_J_S_C_W, features = "LRMDA",raster=FALSE)

####
A_J_S_C_W <- FindClusters(A_J_S_C_W, resolution = 5)
A_J_S_C_W <- RunUMAP(A_J_S_C_W, dims = 1:20, reduction = "integrated.Harmony")
saveRDS(A_J_S_C_W,"A_J_S_C_W.rds")

A_J_S_C_W$cluster_color <- ifelse(A_J_S_C_W$seurat_clusters == "84", "gray", "red")

# 绘制 UMAP 图，使用新的颜色映射
DimPlot(A_J_S_C_W, group.by = "cluster_color")

####patient remove####
Idents(A_J_S_C_W) <- "seurat_clusters"
A_J_S_C_W_removed <- A_J_S_C_W[, !(colnames(A_J_S_C_W) %in% WhichCells(A_J_S_C_W, idents = c("41","84","102","105","106","107","108")))]
A_J_S_C_W_removed <- FindNeighbors(A_J_S_C_W_removed, reduction = "integrated.Harmony", dims = 1:20)
A_J_S_C_W_removed <- FindClusters(A_J_S_C_W_removed, resolution = 5)
A_J_S_C_W_removed <- RunUMAP(A_J_S_C_W_removed, dims = 1:20, reduction = "integrated.Harmony")
DimPlot(A_J_S_C_W_removed, reduction = "umap",raster=FALSE,label = T)
saveRDS(A_J_S_C_W_removed,"A_J_S_C_W_removed.rds")

######
A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "87", "gray", "red")
# A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "85", "gray", "red")
# A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "86", "gray", "red")
A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "89", "gray", "red")
A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "92", "gray", "red")
A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "94", "gray", "red")
# A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "95", "gray", "red")
# A_J_S_C_W_removed$cluster_color <- ifelse(A_J_S_C_W_removed$seurat_clusters == "101", "gray", "red")
# 绘制 UMAP 图，使用新的颜色映射
DimPlot(A_J_S_C_W_removed, group.by = "cluster_color",raster=FALSE)


######remove doublets#####
Idents(A_J_S_C_W_removed) <- "seurat_clusters"
A_J_S_C_W_removed_2 <- A_J_S_C_W_removed[, !(colnames(A_J_S_C_W_removed) %in% WhichCells(A_J_S_C_W_removed, idents = c("87","89","92","94")))]
A_J_S_C_W_removed_2 <- FindNeighbors(A_J_S_C_W_removed_2, reduction = "integrated.Harmony", dims = 1:20)
A_J_S_C_W_removed_2 <- FindClusters(A_J_S_C_W_removed_2, resolution = 1)
A_J_S_C_W_removed_2 <- RunUMAP(A_J_S_C_W_removed_2, dims = 1:20, reduction = "integrated.Harmony")
DimPlot(A_J_S_C_W_removed_2, reduction = "umap",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_removed_2, reduction = "umap",group.by = "cell_types_origin",raster=FALSE,label = T)
saveRDS(A_J_S_C_W_removed_2,"A_J_S_C_W_removed_2.rds")
FeaturePlot(A_J_S_C_W_removed_3, features = "HAP1",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "CACNA1G",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "C10orf11",raster=FALSE)
FeaturePlot(A_J_S_C_W_removed_3, features = "FYB",raster=FALSE)
######remove doublets#####
Idents(A_J_S_C_W_removed_2) <- "seurat_clusters"
A_J_S_C_W_removed_3 <- A_J_S_C_W_removed_2[, !(colnames(A_J_S_C_W_removed_2) %in% WhichCells(A_J_S_C_W_removed_2, idents = "24"))]
A_J_S_C_W_removed_3 <- FindNeighbors(A_J_S_C_W_removed_3, reduction = "integrated.Harmony", dims = 1:20)
A_J_S_C_W_removed_3 <- FindClusters(A_J_S_C_W_removed_3, resolution = 1)
A_J_S_C_W_removed_3 <- RunUMAP(A_J_S_C_W_removed_3, dims = 1:20, reduction = "integrated.Harmony")
DimPlot(A_J_S_C_W_removed_3, reduction = "umap",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_removed_3, split.by = "orig.ident", reduction = "umap",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_removed_3, reduction = "umap",group.by = "cell_types_origin",raster=FALSE,label = T)

######annotation#####
A_J_S_C_W_removed_3$cell_types <- "Ng"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(0,1,2,4,5,12,14,18)] <- "Oligo"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(19,15)] <- "In neuron"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(3,9,13,23,29,32)] <- "Ex neuron"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(6,33,34)] <- "OPC"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(17,22,27)] <- "Vascular"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(7,8,20,21,26,30,31)] <- "Astrocyte"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(10,11,16,25)] <- "Microglia"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(24)] <- "T cell"
A_J_S_C_W_removed_3$cell_types[A_J_S_C_W_removed_3@meta.data$seurat_clusters%in%c(28)] <- "B cell"
DimPlot(A_J_S_C_W_removed_3, reduction = "umap",group.by = "cell_types",raster=FALSE,label = T)
saveRDS(A_J_S_C_W_removed_3,"A_J_S_C_W_removed_3.rds")
A_J_S_C_W_removed_3$Condition[A_J_S_C_W_removed_3$Donor_ID %in% c("14_043","11_069","12_002")] <- "CTRL"
A_J_S_C_W_removed_3@meta.data$MS_type[A_J_S_C_W_removed_3$Donor_ID %in% c("14_043","11_069","12_002")] <- "ctrl"


pathology_sample_counts <- A_J_S_C_harmony@meta.data %>%
  group_by(Pathology) %>%                  # 按 pathology 分组
  summarise(n_samples = n_distinct(Sample)) # 计算每组中不同 sample 的数量

# 查看结果
print(pathology_sample_counts,n = 50)


pathology_per_dataset_str <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(Pathology) %>%
  summarise(
    Pathologies = paste(unique(orig.ident), collapse = ", "),
    .groups = "drop"
  )

# 打印结果
print(pathology_per_dataset_str)

#####################################
Donor_ID_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(Condition,Sex) %>%                  # 按 pathology 分组
  summarise(n_samples = n_distinct(Donor_ID)) # 计算每组中不同 sample 的数量

# 查看结果
print(Donor_ID_counts,n = 50)

#####################################
Donor_ID_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(Pathology) %>%                  # 按 pathology 分组
  summarise(n_samples = n()) # 计算每组中不同 sample 的数量

# 查看结果
print(Donor_ID_counts,n = 50)


######Bar plot########
# 计算每种病理类型的病人数量
donor_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(Pathology,Sex) %>%
  summarise(Num_Patients = n_distinct(Donor_ID))

# 绘制柱状图
ggplot(donor_counts, aes(x = Pathology, y = Num_Patients, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Pathology") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )

# 计算每种性别的病人数量
donor_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(Sex) %>%
  summarise(Num_Patients = n_distinct(Donor_ID))

# 绘制柱状图
ggplot(donor_counts, aes(x = Sex, y = Num_Patients)) +
  geom_bar(stat = "identity", fill = "black") +
  xlab("Sex") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 12),
        panel.grid.major = element_blank(),                             # 移除主网格线
        axis.line = element_line(color = "black")
  )

########## 计算MS年龄的病人数量
MS_A_J_S_C_W_removed_3 <- subset(A_J_S_C_W_removed_3, subset = Condition == "MS")
donor_counts <- MS_A_J_S_C_W_removed_3@meta.data %>%
  group_by(Age,Sex) %>%
  summarise(Num_Patients = n_distinct(Donor_ID))

# 绘制柱状图
donor_counts$Age <- as.numeric(as.character(donor_counts$Age))

# 创建年龄区间
# 创建年龄区间，并为每个区间指定标签
donor_counts$Age_Group <- cut(donor_counts$Age, 
                              breaks = seq(min(donor_counts$Age), max(donor_counts$Age), by = 5), 
                              include.lowest = TRUE, 
                              right = FALSE, 
                              labels = paste(seq(min(donor_counts$Age), max(donor_counts$Age) - 5, by = 5), 
                                             seq(min(donor_counts$Age) + 5, max(donor_counts$Age), by = 5), sep = "-"))

# 按年龄区间分组并计算患者数量
donor_counts_grouped <- donor_counts %>%
  group_by(Age_Group,Sex) %>%
  summarise(Num_Patients = sum(Num_Patients))

# 绘制柱状图
ggplot(donor_counts_grouped, aes(x = Age_Group, y = Num_Patients, fill = Sex)) +
  geom_bar(stat = "identity",position = "stack" ) +
  xlab("Age Group") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )

########## 计算CTRL年龄的病人数量
CTRL_A_J_S_C_W_removed_3 <- subset(A_J_S_C_W_removed_3, subset = Condition == "CTRL")
donor_counts <- CTRL_A_J_S_C_W_removed_3@meta.data %>%
  group_by(Age,Sex) %>%
  summarise(Num_Patients = n_distinct(Donor_ID))

# 绘制柱状图
donor_counts$Age <- as.numeric(as.character(donor_counts$Age))

# 创建年龄区间
# 创建年龄区间，并为每个区间指定标签
donor_counts$Age_Group <- cut(donor_counts$Age, 
                              breaks = seq(min(donor_counts$Age), max(donor_counts$Age), by = 5), 
                              include.lowest = TRUE, 
                              right = FALSE, 
                              labels = paste(seq(min(donor_counts$Age), max(donor_counts$Age) - 5, by = 5), 
                                             seq(min(donor_counts$Age) + 5, max(donor_counts$Age), by = 5), sep = "-"))

# 按年龄区间分组并计算患者数量
donor_counts_grouped <- donor_counts %>%
  group_by(Age_Group,Sex) %>%
  summarise(Num_Patients = sum(Num_Patients))

# 绘制柱状图
ggplot(donor_counts_grouped, aes(x = Age_Group, y = Num_Patients, fill = Sex)) +
  geom_bar(stat = "identity",position = "stack" ) +
  xlab("Age Group") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )



# 计算每种数据来源的病人数量
donor_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(orig.ident, Condition) %>%
  summarise(Num_Patients = n_distinct(Donor_ID), .groups = "drop")

# 绘制柱状图
ggplot(donor_counts, aes(x = orig.ident, y = Num_Patients, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("orig.ident") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )

# 计算每种Condition的病人数量
donor_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(Condition,Sex) %>%
  summarise(Num_Patients = n_distinct(Donor_ID))

# 绘制柱状图
ggplot(donor_counts, aes(x = Condition, y = Num_Patients, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Condition") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )



######
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "Oligo", "Oligo","other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "OPC", "OPC", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "Microglia", "Microglia", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "Astrocyte", "Astrocyte", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "Ex neuron", "Ex neuron", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "In neuron", "In neuron", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "T cell", "T cell", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "B cell", "B cell", "other")
A_J_S_C_W_removed_3$group_color <- ifelse(A_J_S_C_W_removed_3$cell_types == "Vascular", "Vascular", "other")

DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("Oligo" = "lightblue", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("OPC" = "#6495ED", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("Microglia" = "#32CD32", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("Astrocyte" = "lightcoral", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("Ex neuron" = "darkgreen", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("In neuron" = "lightgreen", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("T cell" = "purple", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("B cell" = "brown", "other" = "lightgray"), raster = FALSE)
DimPlot(A_J_S_C_W_removed_3, group.by = "group_color", cols = c("Vascular" = "#FF1493", "other" = "lightgray"), raster = FALSE)

''' MARKER: 
Oligo："PLP1","ST18","MOBP"
OPC:"PCDH15","TNR","LHFPL3"
Microglia:"LRMDA","DOCK8","TBXAS1"
Astrocyte:SLC1A2,GLIS3,GFAP
Ex neuron:SV2B,KCNIP4,LRRTM4
In neuron:GAD1,VIP
T cell:PTPRC,SKAP1,PARP8
B cell:IGHG1,IGKC,IGHG3
Vascular:"FLT1", "EBF1", "CLDN5" '''
FeaturePlot(A_J_S_C_W_removed_3, features = "KCNIP4",raster=FALSE)
DotPlot(A_J_S_C_W_removed_3, features = c("PLP1","ST18","MOBP",
                                          "LRMDA","DOCK8","TBXAS1",
                                          "SLC1A2","GLIS3","GFAP",
                                          "SLC17A7","SV2B","CBLN2",
                                          "PCDH15","TNR","LHFPL3",
                                          "FLT1", "EBF1", "CLDN5",
                                          "PTPRC","SKAP1","PARP8",
                                          "GAD1","GAD2","GALNTL6",
                                          "IGHG1","IGKC","IGHG3")) +
  scale_color_gradientn(colors = c("#FFFFE0", "#FFD700", "#FF4500", "#800000"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 计算不同cell_typed的病理类型的病人数量
donor_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(cell_types,Pathology) %>%
  summarise(Num_Patients = n_distinct(Donor_ID))

# 绘制柱状图
ggplot(donor_counts, aes(x = cell_types, y = Num_Patients, fill = Pathology)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Pathology") +
  ylab("Number of Donors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )







cell_counts <- A_J_S_C_W_removed_3@meta.data %>%
  group_by(cell_types, Pathology) %>%
  summarise(Num_Cells = n(), .groups = "drop")

# 计算总细胞数量
total_cell_counts <- cell_counts %>%
  group_by(cell_types) %>%
  summarise(Total_Cells = sum(Num_Cells), .groups = "drop")

# 计算百分比
cell_percentages <- cell_counts %>%
  left_join(total_cell_counts, by = "cell_types") %>%
  mutate(Percentage = (Num_Cells / Total_Cells) * 100)

# 绘制柱状图显示百分比
ggplot(cell_percentages, aes(x = cell_types, y = Percentage, fill = Pathology)) +
  geom_bar(stat = "identity", position = "fill") +  # 使用position = "fill"来归一化堆叠高度
  xlab("Cell Types") +
  ylab("Percentage of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
  )