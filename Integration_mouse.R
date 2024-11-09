library(dplyr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(clustree)
dir1 = c('DRG.matrix/')
dir2 = c("SJXW.matrix/")
names(dir1) = c('DRG')  
names(dir2) = c( 'Tumor')  
counts1 <- Read10X(data.dir =dir1)
counts2 <- Read10X(data.dir =dir2)
sc_drg = CreateSeuratObject(counts1,min.cells = 3, min.features = 50,project = "DRG")
sc_tumor = CreateSeuratObject(counts2,min.cells = 3, min.features = 50,project = "Tumor")
sc_drg[["percent.mt"]] <- PercentageFeatureSet(sc_drg, pattern = "^mt-")
sc_tumor[["percent.mt"]] <- PercentageFeatureSet(sc_tumor, pattern = "^mt-")
#计算红细胞比例
HB.genes <- c("Hba-a1","Hba-a2", "Hbb-bs", "Hbb-bt")
HB_m <- match(HB.genes, rownames(sc_drg@assays$RNA)) 
HB.genes <- rownames(sc_drg@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
sc_drg[["percent.HB"]]<-PercentageFeatureSet(sc_drg, features=HB.genes) 
HB.genes <- c("Hba-a1","Hba-a2", "Hbb-bs", "Hbb-bt")
HB_m <- match(HB.genes, rownames(sc_tumor@assays$RNA)) 
HB.genes <- rownames(sc_tumor@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
sc_tumor[["percent.HB"]]<-PercentageFeatureSet(sc_tumor, features=HB.genes)

VlnPlot(sc_drg,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 4)
VlnPlot(sc_tumor,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 4)
sc_drg <- subset(sc_drg,subset = nFeature_RNA > 200& nFeature_RNA < 6500 & percent.mt < 15 & percent.HB < 3 & nCount_RNA < 50000)
sc_tumor <- subset(sc_tumor,subset = nFeature_RNA > 200& nFeature_RNA < 6500 & percent.mt < 15 & percent.HB < 3 & nCount_RNA < 50000)
rm(counts1,counts2)


mu_list <- list(sc_drg,sc_tumor)
mu_list <- lapply(X = mu_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = mu_list,nfeatures = 2000)
mu_list <- lapply(X = mu_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
mu.anchors <- FindIntegrationAnchors(object.list = mu_list, anchor.features = features, reduction = "rpca",k.anchor = 5)

mu.combined <- IntegrateData(anchorset = mu.anchors)
rm(mu.anchors)
DefaultAssay(mu.combined) <- "integrated"

mu.combined <- ScaleData(mu.combined, verbose = FALSE)
mu.combined <- RunPCA(mu.combined, npcs = 30, verbose = FALSE)
ElbowPlot(mu.combined, ndims=30)
mu.combined <- FindNeighbors(mu.combined, reduction = "pca", dims = 1:20)
mu.combined <- FindClusters(mu.combined, resolution = seq(0.1,1.2,by=0.1))
clustree(mu.combined)

mu.combined <- RunUMAP(mu.combined, reduction = "pca", dims = 1:20)
Idents(mu.combined) = mu.combined@meta.data$integrated_snn_res.0.8

plot1 =DimPlot(mu.combined, reduction = "umap",label = T) 
plot2 = DimPlot(mu.combined, reduction = "umap", group.by='orig.ident') 
plot1+plot2

FeaturePlot(mu.combined,features = c("Fosb","Gfra3","Aldh1a1","Vegfa"))

my_color <- c("#ffa726","#ef9a9a","#673ab7","#00b0ff","#ab47bc","#b71c1c",
                       "#a1887f","#5c6bc0","#9ccc65","#4caf50",
                       "#C653A5","#e7298a")
DimPlot(mu.combined,reduction = "umap",group.by = "integrated_snn_res.0.2",label = T,cols = my_color)
DimPlot(mu.combined, reduction = "umap", group.by='orig.ident',cols = c("#ffa726","#00b0ff"))

mu_all_mark <- FindAllMarkers(mu.combined,only.pos = T,logfc.threshold = 0.3,min.pct = 0.2)
mu_top10 <- mu_all_mark %>% group_by(cluster) %>% top_n(10,wt=avg_log2FC)

mu.combined = RenameIdents(mu.combined,"0"="Immature_SCs","15"="Mye_SCs","14"="Cycling_SCs","9"="Alah1a1_SCs",
                           "2"="iCAF","3"="myCAF","10"="Dpt_CAF",
                           "12"="Macrophages","1"="Macrophages","8"="Macrophages",
                           "4"="Neutrophils","6"="Neutrophils","7"="Neutrophils","11"="T_cells",
                           "5"="B_cells","13"="Endo_cells")
#mu.combined = RenameIdents(mu.combined,"0"="Immature_SCs","10"="Mye_SCs","11"="Cycling_SCs","1"="iCAF","5"="myCAF",
#                                        "2"="Macrophages","4"="Macrophages","3"="Neutrophils","7"="Neutrophils","8"="T_cells",
 #                          "6"="B_cells","9"="Endo_cells")
mu.combined@meta.data$celltype = mu.combined@active.ident
DimPlot(mu.combined,reduction = "umap",group.by = "celltype",cols = my_color)


###########celltype markers heatmap#########
mu_type_mark <- FindAllMarkers(mu.combined,only.pos = T,logfc.threshold = 0.3,min.pct = 0.2)
mk_sig <- mu_type_mark %>% 
  filter(pct.1 > 0.3 &
           (pct.1-pct.2) > 0.1) %>%
  filter(cluster != "other") %>%
  arrange(cluster, -(pct.1-pct.2)) # 按power排序，这是热图呈好看的对角线排列的关键
table(mk_sig$cluster)
degs_top30<- mk_sig %>% 
  group_by(cluster) %>% 
  top_n(30, (pct.1-pct.2)) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, -(pct.1-pct.2))
write.csv(degs_top30,file = "celltype_deg30.csv")
avgData <- mu.combined@assays$integrated@data[degs_top30$gene,] %>% 
  apply(1, function(x){
    tapply(x, mu.combined@active.ident, mean) # ExpMean
  }) %>% t

phData <- MinMax(scale(avgData), -2, 2) # z-score
rownames(phData) <- 1:nrow(phData)
library(pheatmap)
ann_colors = list(
    cluster = c(Immature_SCs="#ffa726",Mye_SCs="#ef9a9a",Cycling_SCs="#673ab7",
                iCAF="#00b0ff",myCAF="#ab47bc",Macrophages="#b71c1c",
                Neutrophils= "#a1887f",T_cells="#5c6bc0",B_cells="#9ccc65",Endo_cells="#4caf50")
)
phres <- pheatmap(
  phData, 
  color = colorRampPalette(c("#b3e5fc", "white", "red3"))(99), #配色
  scale = "row",
  cluster_rows = F, #不按行聚类
  cluster_cols = F, #按列聚类
  clustering_method = "complete",
  show_rownames = F, #显示cluster名
  annotation_row = data.frame(cluster = degs_top30$cluster) ,
  annotation_colors = ann_colors,
  annotation_legend = F
)

########2 samples column proportion#####
sc_df <- table(mu.combined@meta.data$celltype,mu.combined@meta.data$orig.ident) %>% melt()
colnames(sc_df) <- c("Cluster","Sample","Number")
sc_df$Sample = c(rep("DRG",10),rep("Dermal",10))
#sc_df <- sc_df[sc_df$Number>0,]
pB5 <- ggplot(data = sc_df, aes(x =Sample, y =Number , fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = my_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black",angle = 45,hjust = 1))
pB5
pB6 <- ggplot(data = sc_df, aes(x =Number, y =Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = my_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black",angle = 45,hjust = 1))
pB6

###save new S4##
save(mu.combined,file = "mu.combined230107.RData")
load(file = "mu.combined230107.RData")
###SCs####
t_sc <- subset(mu.combined,idents = c("Immature_SCs","Mye_SCs","Alah1a1_SCs"))
t_sc <- subset(mu.combined,idents = c("0","5","15"))
t_sc = ScaleData(t_sc)
FeaturePlot(t_sc,features = c("Vegfa","Igfbp5","Aldh1a1"))
FeaturePlot(t_sc,features = c("Amotl2","Ccn1","Ccn2"))
t_sc <- SCTransform(t_sc)
t_sc <- FindNeighbors(t_sc, reduction = "pca", dims = 1:20)
t_sc <- FindClusters(t_sc, resolution = seq(0.1,1,by=0.1))
clustree(t_sc)
t_sc <- RunUMAP(t_sc,dims = 1:20)
DimPlot(t_sc)
Idents(t_sc) = "integrated_snn_res.0.4"
DimPlot(t_sc,group.by = "integrated_snn_res.0.4")
t_sc_1 <- subset(t_sc,idents = 0:5)
Idents(t_sc_1) = "celltype"
t_sc_1 <- RunUMAP(t_sc_1,dims = 1:20)
DimPlot(t_sc_1,label = F,cols = c("#ffa726","#ef9a9a","#00b0ff"))
DimPlot(t_sc_1, reduction = "umap", group.by='orig.ident',cols = c("#ffa726","#00b0ff"))
FeaturePlot(t_sc_1,features = c("Sox10","S100b","Mbp","Mpz","Ncam1","Lamc1","Gpr37l1","Cadm4"),ncol = 4)
FeaturePlot(t_sc_1,features = c("Amotl2","Ccn1","Ccn2","Ccnd1"))
VlnPlot(t_sc_1,features = c("Amotl2","Ccn1","Ccn2","Ccnd1"),ncol=4)
t_sc_mark <- FindAllMarkers(t_sc,only.pos = T,logfc.threshold = 0.3,min.pct = 0.2)
FeaturePlot(sc_tumor,features = c("Vegfa","Bgn","Rpl10","Rps25"))
####hippo score###
hippo_gene <- list(c("Amotl2","Ccn1","Ccn2", "Ccnd1"))
t_sc_1 <- AddModuleScore(
  object = t_sc_1,
  features = hippo_gene,
  ctrl = 50,
  name = "Hippo"
)
VlnPlot(t_sc_1,features = "Hippo1",sort = T)
FeaturePlot(t_sc_1,features = "Hippo1",min.cutoff = -0.5,max.cutoff = 1)

##DE##
sw_de <- FindMarkers(t_sc_1,ident.1 ="Mye_SCs")
sw_de_sig <- sw_de %>% 
  filter(pct.1 > 0.3 &
           abs((pct.1-pct.2)) > 0.2) %>%
  arrange(-abs((pct.1-pct.2))) # 按power排序，这是热图呈好看的对角线排列的关键

de_top20<- sw_de_sig %>% 
  top_n(20, abs((pct.1-pct.2))) %>%
  top_n(20, avg_log2FC) %>%
  arrange(-abs((pct.1-pct.2)))
DotPlot(t_sc_1,features=rownames(de_top20)) + coord_flip()
DotPlot(t_sc_1,features=c("Pllp", "Prx", "Dusp15", "Mpz", "Drp2",
                          "Ccn1", "Amotl2", "Vwa1", "Lgals3","Itga6")) + coord_flip() + scale_color_gradient(low = "grey",high = "red3")+
                          theme(axis.text.x = element_text(size=12, colour = "black",angle = 45,hjust = 1))
####
save(t_sc_1,file = "t_sc_120230110.Rdata")
load("t_sc_120230110.Rdata")

##
