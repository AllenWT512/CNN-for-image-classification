Sys.setenv(LANGUAGE="en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
library(Seurat)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(dplyr) 

A_v5_new_patient_removed <- readRDS("/PROJ01/share/zhouwutong/AJSCW/A/A_v5_new_patient_removed.rds")
J_v5_new_patient_removed <- readRDS("/PROJ01/share/zhouwutong/AJSCW/J/J_v5_new_patient_removed.rds")
S_v5_new <- readRDS("/PROJ01/share/zhouwutong/AJSCW/S/S_v5_new.rds")
C_integration_16_harmony <- readRDS("/PROJ01/share/zhouwutong/AJSCW/Nat_neuro_MS_dataset_2024/C_integration_16_harmony.rds")
W_revise <- readRDS("/PROJ01/share/zhouwutong/AJSCW/W/W_revise.rds")

DefaultAssay(A_v5_new_patient_removed) <- "RNA"
DefaultAssay(J_v5_new_patient_removed) <- "RNA"
DefaultAssay(S_v5_new) <- "RNA"
DefaultAssay(C_integration_16_harmony) <- "RNA"
DefaultAssay(W_revise) <- "RNA"

# dataset integration
A_J_S_C_W <- merge(C_integration_16_harmony, y = list(A_v5_new_patient_removed,J_v5_new_patient_removed, S_v5_new,W_revise))
A_J_S_C_W@assays[["RNA"]] = as(object = A_J_S_C_W@assays[["RNA"]], Class = "Assay5")

# pipeline
A_J_S_C_W <- NormalizeData(A_J_S_C_W)
A_J_S_C_W <- FindVariableFeatures(A_J_S_C_W, assay = "RNA", selection.method = "vst", nfeatures = 2000)
A_J_S_C_W <- ScaleData(A_J_S_C_W)
A_J_S_C_W <- RunPCA(A_J_S_C_W, verbose=F)
ElbowPlot(A_J_S_C_W)
# harmony去批次
A_J_S_C_W_harmony <- IntegrateLayers(object = A_J_S_C_W, 
                                     method = HarmonyIntegration, 
                                     orig.reduction = "pca", 
                                     new.reduction = "integrated.Harmony",
                                     max.iter.cluster = 100,
                                     verbose = FALSE)
A_J_S_C_W_harmony[["RNA"]] <- JoinLayers(A_J_S_C_W_harmony[["RNA"]])

# 细胞聚类
set.seed(42)
A_J_S_C_W_harmony <- FindNeighbors(A_J_S_C_W_harmony, reduction = "integrated.Harmony", dims = 1:20)
A_J_S_C_W_harmony <- FindClusters(A_J_S_C_W_harmony, resolution = 0.5)
A_J_S_C_W_harmony <- RunUMAP(A_J_S_C_W_harmony, dims = 1:20, reduction = "integrated.Harmony")
DimPlot(A_J_S_C_W_harmony, reduction = "umap", split.by = "orig.ident",raster=FALSE)
saveRDS(object = A_J_S_C_W_harmony,file = "/PROJ01/share/zhouwutong/AJSCW/A_J_S_C_W_harmony.rds",compress = FALSE)

#####change resolution to 1
A_J_S_C_W_harmony <- FindClusters(A_J_S_C_W_harmony, resolution = 1)
A_J_S_C_W_harmony <- RunUMAP(A_J_S_C_W_harmony, dims = 1:20, reduction = "integrated.Harmony")
DimPlot(A_J_S_C_W_harmony, reduction = "umap", raster=FALSE,label = T)

#####change resolution to 1.5
A_J_S_C_W_harmony <- FindClusters(A_J_S_C_W_harmony, resolution = 1)
A_J_S_C_W_harmony <- RunUMAP(A_J_S_C_W_harmony, dims = 1:20, reduction = "integrated.Harmony")
DimPlot(A_J_S_C_W_harmony, reduction = "umap", raster=FALSE,label = T)

A_J_S_C_W_harmony@meta.data$cell_type_new = NA
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("0","1","2","4","5","8","12","15")] <- "Oligo"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("24","36","23","35","9","18","26","10","28")] <- "Astrocyte"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("19","22")] <- "Vascular"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("25")] <- "T cell"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("30")] <- "B cell"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("7","37","38")] <- "OPCs+COPs"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("27","20","13","6","29")] <- "Microglia"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("14","21")] <- "Inhibitory neurons"
A_J_S_C_W_harmony@meta.data$cell_type_new[A_J_S_C_W_harmony$seurat_clusters %in% c("33","34","16","11","3","32","31","17")] <- "Excitatory neurons"

A_J_S_C_W_harmony@meta.data$disease[A_J_S_C_W_harmony$disease %in% c("PPMS","RRMS","SPMS")] <- "ms"
A_J_S_C_W_harmony@meta.data$disease[A_J_S_C_W_harmony$disease %in% c("CTR")] <- "ctrl"

A_J_S_C_W_harmony@meta.data$pathology[A_J_S_C_W_harmony$pathology %in% c("ctrl_GM")] <- "ctrl"
A_J_S_C_W_harmony@meta.data$pathology[A_J_S_C_W_harmony$pathology %in% c("ctrl_WM")] <- "ctrl"

A_J_S_C_W_harmony@meta.data$region[A_J_S_C_W_harmony$region %in% c("white_matter")] <- "WM"
A_J_S_C_W_harmony@meta.data$region[A_J_S_C_W_harmony$region %in% c("white&grey_matter")] <- "WM&GM"

DimPlot(A_J_S_C_W_harmony, reduction = "umap", group.by = "cell_type_new",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_harmony, reduction = "umap", split.by = "orig.ident",group.by = "cell_type_new",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_harmony, reduction = "umap", split.by = "sex",group.by = "cell_type_new",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_harmony, reduction = "umap", split.by = "region",group.by = "cell_type_new",raster=FALSE,label = T)
DimPlot(A_J_S_C_W_harmony, reduction = "umap", split.by = "disease",group.by = "cell_type_new",raster=FALSE,label = T)

AllMarkers = FindAllMarkers(A_J_S_C_W_harmony)






















# 绘图
DimPlot(A_J_S_C_W_harmony, reduction = "umap", split.by = "orig.ident",raster=FALSE)

DimPlot(A_J_S_C_W_harmony, reduction = "umap", raster=FALSE, label = T)

A_J_S_C_harmony_all_markers <- FindAllMarkers(A_J_S_C_harmony, only.pos = TRUE)

###Oligo_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "ST18",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "TMEM144",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "TF",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "MBP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "MOBP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "PLP1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "CNP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "MAG",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "MYRF",raster=FALSE)

###Astrocyte_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "GLIS3",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "PITPNC1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "NRG3",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "RFX4",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "TPD52L1",raster=FALSE)

###IMM_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "ARHGAP24",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "DOCK8",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "TBXAS1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "APBB1IP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "ARHGAP15",raster=FALSE)

###vascular_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "EPAS1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "PITPNC1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "CLDN5",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "COBLL1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "EBF1",raster=FALSE)

###Neuro_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "NELL2",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "SYT1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "RBFOX1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GABRB2",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GRIN2B",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GRIN2A",raster=FALSE)

###EX Neuro_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "SYT1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "RBFOX1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GABRB2",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GRIN2B",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GRIN2A",raster=FALSE)

###In Neuro_marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "DAB1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GAD1",raster=FALSE)

###OPC marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "PDGFRA",raster=FALSE)

###Neuro progenitor marker############################################################
FeaturePlot(A_J_S_C_harmony, features = "BMPR1B",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "ANGPT1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "LAMA1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "SPAG17",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GFAP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "ANGPT1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "LAMA1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "SPAG17",raster=FALSE)

FeaturePlot(A_J_S_C_harmony, features = "PTPRC",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "ARHGAP15",raster=FALSE)

####Tcell marker####
FeaturePlot(A_J_S_C_harmony_removed_annotation, features = "SKAP1",raster=FALSE)

A_J_S_C_harmony_annotation = RenameIdents(A_J_S_C_harmony,
                                "0" = "Oligodendrocyte",
                                "1" = "Oligodendrocyte",
                                "2" = "Astrocyte",
                                "3" = "Neuro",
                                "4" = "IMM",
                                "5" = "OPC",
                                "6" = "Neuro",
                                "7" = ,
                                "8" = "Vascular",
                                "9" = "Neuro",
                                "10" = "T cell",
                                "11" = "Neuro progenitor")

saveRDS(A_v5_new_patient_removed,"A_v5_new_patient_removed.rds")
saveRDS(J_v5_new_patient_removed,"J_v5_new_patient_removed.rds")
write.csv(A_J_S_C_harmony_all_markers, "A_J_S_C_harmony_all_markers.csv")

FeaturePlot(A_J_S_C_harmony, features = "ABCA9",raster=FALSE) ###SC##
FeaturePlot(A_J_S_C_harmony, features = "SKAP1",raster=FALSE) ##T cell###
FeaturePlot(A_J_S_C_harmony, features = "FcγRIII",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "MOBP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "GFAP",raster=FALSE) ##reactive AS##
FeaturePlot(A_J_S_C_harmony, features = "AQP4",raster=FALSE) ##AS general##
FeaturePlot(A_J_S_C_harmony, features = "IRF7",raster=FALSE) ##AS general##
FeaturePlot(A_J_S_C_harmony, features = "ARHGAP24",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "DOCK8",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "CENP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "P2RY12",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "SOX6",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "VCAN",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "PCDH15",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "CD19",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "CD38",raster=FALSE)
####70 ：neuro&microglia 83：OPC & P2RY12###  remove 70,83

FeaturePlot(A_J_S_C_harmony, features = "IQGAP2",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "ADAMTSL1",raster=FALSE)
### remove 82

FeaturePlot(A_J_S_C_harmony, features = "HIF3A",raster=FALSE)


######
Idents(A_J_S_C_harmony) <- "RNA_snn_res.4"
A_J_S_C_harmony_removed <- subset(A_J_S_C_harmony, idents = c("70", "80", "82","83"),invert = TRUE)
A_J_S_C_harmony_removed <- FindNeighbors(A_J_S_C_harmony_removed, reduction = "integrated.Harmony", dims = 1:20)
A_J_S_C_harmony_removed <- FindClusters(A_J_S_C_harmony_removed, resolution = 0.7)
A_J_S_C_harmony_removed <- RunUMAP(A_J_S_C_harmony_removed, dims = 1:20, reduction = "integrated.Harmony")
A_J_S_C_harmony_removed_all_markers <- FindAllMarkers(A_J_S_C_harmony_removed, only.pos = TRUE)


# 绘图
DimPlot(A_J_S_C_harmony, reduction = "umap", split.by = "orig.ident",raster=FALSE)
DimPlot(A_J_S_C_harmony, reduction = "umap", group.by = "orig.ident",raster=FALSE)
DimPlot(A_J_S_C_harmony_removed, reduction = "umap", raster=FALSE, label = T)

####Astrocyte marker####
FeaturePlot(A_J_S_C_harmony_removed_annotation, features = "SKAP1",raster=FALSE)

FeaturePlot(A_J_S_C_harmony_removed, features = "CEMIP",raster=FALSE) ##SC cell##
FeaturePlot(A_J_S_C_harmony_removed, features = "MOBP",raster=FALSE)

FeaturePlot(A_J_S_C_harmony_removed, features = "EPAS1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "PITPNC1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "CLDN5",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "COBLL1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony, features = "EBF1",raster=FALSE)

FeaturePlot(A_J_S_C_harmony_removed, features = "PDGFRA",raster=FALSE)


FeaturePlot(A_J_S_C_harmony_removed, features = "DAB1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "GAD1",raster=FALSE)

FeaturePlot(A_J_S_C_harmony_removed, features = "GLIS3",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "PITPNC1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "NRG3",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "RFX4",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "TPD52L1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "GFAP",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "AQP4",raster=FALSE)

FeaturePlot(A_J_S_C_harmony_removed, features = "SYT1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "RBFOX1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "GABRB2",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "GRIN2B",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "GRIN2A",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "CAMK2A",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "NRGN",raster=FALSE)


###Neuro progenitor marker############################################################
FeaturePlot(A_J_S_C_harmony_removed, features = "BMPR1B",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "ANGPT1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "LAMA1",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "SPAG17",raster=FALSE)
Idents(A_J_S_C_harmony) <- "RNA_snn_res.0.7"
A_J_S_C_harmony_removed_annotation = RenameIdents(A_J_S_C_harmony_removed,
                                          "0" = "Oligodendrocyte",
                                          "1" = "Oligodendrocyte",
                                          "2" = "Oligodendrocyte",
                                          "3" = "Oligodendrocyte",
                                          "4" = "Oligodendrocyte",
                                          "5" = "Astrocyte",
                                          "6" = "Microglia",
                                          "7" = "OPC",
                                          "8" = "Astrocyte",
                                          "9" = "Ex neuro",
                                          "10" = "Ex neuro",
                                          "11" = "Oligodendrocyte",
                                          "12" = "Microglia",
                                          "13" = "Astrocyte",
                                          "14" = "In Neuro",
                                          "15" = "In Neuro",
                                          "16" = "Ex neuro",
                                          "17" = "Ex neuro",
                                          "18" = "Vascular",
                                          "19" = "Ex neuro",
                                          "20" = "T cell",
                                          "21" = "Astrocyte",
                                          "22" = "Neuro progenitor",
                                          "23" = "Stromal cell",
                                          "24" = "Ex neuro")
A_J_S_C_harmony_removed_annotation$cell_type <- Idents(A_J_S_C_harmony_removed_annotation)
DimPlot(A_J_S_C_harmony_removed_annotation, reduction = "umap", label = T, raster=FALSE)
DimPlot(A_J_S_C_harmony_removed_annotation, group.by = "orig.ident", reduction = "umap", label = T, raster=FALSE)
DimPlot(A_J_S_C_harmony_removed_annotation, split.by = "orig.ident", reduction = "umap", label = T, raster=FALSE)
DimPlot(A_J_S_C_harmony_removed_annotation, split.by = "pathology", reduction = "umap", label = T, raster=FALSE)
DimPlot(A_J_S_C_harmony_removed_annotation, split.by = "sex", reduction = "umap", label = T, raster=FALSE)
DimPlot(A_J_S_C_harmony_removed_annotation, split.by = "region", reduction = "umap", label = T, raster=FALSE)
saveRDS(A_J_S_C_harmony_removed_annotation,"A_J_S_C_harmony_removed_annotation.rds")

write.csv(A_J_S_C_harmony_removed_annotation@meta.data,"A_J_S_C_harmony_removed_annotation@meta.data.csv",  quote = F)
write.csv(C_integration_16_harmony@meta.data,"C_integration_16_harmony@meta.data.csv",  quote = F)
FeaturePlot(A_J_S_C_harmony_removed, features = "S1PR1",raster=FALSE) 
FeaturePlot(A_J_S_C_harmony_removed, features = "NFKBIA",raster=FALSE)
FeaturePlot(A_J_S_C_harmony_removed, features = "CRTC1",raster=FALSE)


