#library(devtools)
#install_github("navinlabcode/copykat")
library(copykat)
#load("mouse2.Rdata")
exp.rawdata <- as.matrix(mu_sc@assays$RNA@counts)
#install.packages("quantreg")
copykat.test <- copykat(rawmat=exp.rawdata, 
                        id.type="S", 
                        ngene.chr=5, 
                        win.size=25, 
                        KS.cut=0.1, 
                        sam.name="mouse", 
                        distance="euclidean", 
                        norm.cell.names="",
                        output.seg="FLASE", 
                        plot.genes="TRUE", 
                        genome="mm10",
                        n.cores=1)
save(copykat.test,file = "mu_sc_copykat_test230107.Rdata")

pred.test <- data.frame(copykat.test$prediction)
table(pred.test$copykat.pred)
pred.test$copykat.pred[pred.test$copykat.pred == "not.defined"] = "diploid"

CNA.test <- data.frame(copykat.test$CNAmat)
my_palette <- colorRampPalette(c("#00b0ff", "white", "red3"))(999)

chr <- as.numeric(CNA.test$chrom) %% 2+1
rbPal1 <- colorRampPalette(c('black','grey'))
CHR <- rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
com.preN <- pred.test$copykat.pred
pred <- rbPal5(2)[as.numeric(factor(com.preN))]

cells <- rbind(pred,pred)
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

heatmap.3(t(CNA.test[,8:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
library(dplyr)
#######tumor cells
tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
tumor.mat <- CNA.test[,which(colnames(CNA.test)[8:ncol(CNA.test)] %in% tumor.cells)]
hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
hc.umap <- cutree(hcc,2)

rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
cells <- rbind(subpop,subpop)

heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =4, method = "euclidean"), 
          hclustfun = function(x) hclust(x, method="ward.D2"),
          ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          keysize=1, density.info="none", trace="none",
          cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))

legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
###
table(mu.combined@meta.data$integrated_snn_res.0.2)
mu_sc@meta.data$copykat <- pred.test$copykat.pred
Idents(mu_sc)= "copykat"
mu_sc@meta.data$copykat[mu_sc@meta.data$copykat  == "not.defined"] = "diploid"
p1 = DimPlot(mu_sc,group.by = "copykat",cols = c("red3","#a1887f"))
p1 
ggsave("202302 finaloutput/S7Ccopykat.pdf",p1,width =5,height = 5)
