
load("hm.comb0207.Rdata")
Idents(hm.combined) = hm.combined@meta.data$integrated_snn_res.0.4
hm.combined = RenameIdents(hm.combined,"1"="NF2_like","4"="NF2_like","16"="Myelinating SCs","12"="NF1_like","6"="Trans_SC","7"="MES_like","2"="CAF","8"="CAF","18"="CAF","13"="Fibroblast",
                           "0"="Macrophage","10"="Macrophage","17"="Macrophage","19"="Macrophage","11"="B_cells",
                           "3"="Neutrophil","5"="NK/T cells",
                           "9"="Endothelial/VSMC","15"="Endothelial/VSMC","14"="Cycling")
hm.combined@meta.data$celltype_b <- hm.combined@active.ident
my_color <- c("#ef9a9a","#b71c1c","#00b0ff","#ab47bc","#673ab7",
                       "#a1887f","#8d6e63","#ffa726","#4db6ac",
                       "#9ccc65","#4caf50","#e7298a","#C653A5")
DimPlot(hm.combined,cols = my_color,label = T)+NoLegend()
Idents(hm.combined) = "newident"
pdim_1 = DimPlot(subset(hm.combined, downsample = 4000),reduction = "umap",split.by = "newident",group.by = "celltype_b",cols = my_color)
ggsave("hm_sample_dim.pdf",pdim_1,width = 12,height = 4)
VlnPlot(hm.combined,features = "CD79A")
pB2_df <- table(hm.combined@meta.data$celltype_b,hm.combined@meta.data$newident) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB4 <- ggplot(data = pB2_df, aes(x =Sample, y =Number , fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = my_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black",angle = 45,hjust = 1))
pB4
ggsave("hm_sample_prop.pdf",pB4,width = 3.6,height = 5)
table(hm.combined@meta.data$orig.ident)

Idents(hm.combined) = "celltype_b"
sc.combined <- subset(hm.combined,idents = c("NF2_like","Myelinating SCs","NF1_like","Trans_SC","MES_like"))

sc.combined <- FindVariableFeatures(sc.combined,nfeatures = 2000)
sc.combined <- ScaleData(sc.combined)

sc.combined <- RunPCA(sc.combined,verbose=FALSE)
ElbowPlot(sc.combined,reduction = "pca",ndims = 15)

sc.combined <- RunUMAP(sc.combined, reduction = "pca", dims = 1:8)
hm_sc_dim = DimPlot(sc.combined, reduction = "umap", group.by='celltype_b',cols = my_color,label = T)+NoLegend() 
ggsave("202302 finaloutput/F7J.hm_sc_dim.pdf",hm_sc_dim,width = 5,height = 5)
FeaturePlot(sc.combined,features= c("MPZ"),min.cutoff = -1)

Idents(sc.combined) = "newident" 
hm_sc_dim_sample = DimPlot(subset(sc.combined, downsample = 500),reduction = "umap",split.by = "newident",group.by = "celltype",cols = my_color)

table(sc.combined@meta.data$newident)

sc_df <- table(sc.combined@meta.data$celltype_b,sc.combined@meta.data$newident) %>% melt()
colnames(sc_df) <- c("Cluster","Sample","Number")
sc_df <- sc_df[sc_df$Number>0,]
pB5 <- ggplot(data = sc_df, aes(x =Sample, y =Number , fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = my_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black",angle = 45,hjust = 1))
pB5
ggsave("202302 finaloutput/SF8D.hm_scprop.pdf",pB5,width = 3.6,height = 5)

Idents(sc.combined) <- sc.combined@meta.data$newident
table(sc.combined@meta.data$newident)
sc.combined <- RenameIdents(sc.combined,"NF1"="humanNF1","mpn"="MPNST")
sc.combined@meta.data$newident = sc.combined@active.ident
sc.combined@meta.data$newident <- factor(sc.combined@meta.data$newident,levels = c("NWW","humanNF2","humanNF1","MPNST"))
DimPlot(sc.combined,reduction = "umap",split.by = "newident",group.by = "celltype",cols = my_color)

FeaturePlot(sc.combined,features= c("MPZ"),min.cutoff = -1)

Idents(sc.combined) = "newident" 
p7k = DimPlot(subset(sc.combined, downsample = 500),reduction = "umap",split.by = "newident",group.by = "celltype_b",cols = my_color)
ggsave("230308 finaloutput/F7K.pdf",p7k,width = 12,height = 4)

sc_df <- table(sc.combined@meta.data$celltype_b,sc.combined@meta.data$newident) %>% melt()
colnames(sc_df) <- c("Cluster","Sample","Number")
sc_df <- sc_df[sc_df$Number>0,]
pB5 <- ggplot(data = sc_df, aes(x =Sample, y =Number , fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values = my_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black",angle = 45,hjust = 1))
pB5
