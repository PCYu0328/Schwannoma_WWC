library(monocle)
data <- as(as.matrix(mu_sc@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = mu_sc@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
#print(head(fData(monocle_cds)))

HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.2 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)


HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')


HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "celltype3")+  
  scale_colour_manual(
    values = col4
   )

plot_cell_trajectory(HSMM, color_by = "integrated_snn_res.0.4")
?plot_cell_trajectory

plot_cell_trajectory(HSMM, color_by = "State")

plot_cell_trajectory(HSMM, color_by = "Pseudotime")


plot_cell_trajectory(HSMM, color_by = "State") +
  facet_wrap(~State, nrow = 1)



blast_genes <- row.names(subset(fData(HSMM),
                                gene_short_name %in% c("Gapdh", "Sox17")))
plot_genes_jitter(HSMM[blast_genes,],
                  grouping = "integrated_snn_res.0.4",
                  min_expr = 0.1)


#筛选影响轨迹的基因
HSMM_expressed_genes <-  row.names(subset(fData(HSMM),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("Ccn1", "Gapdh", "Plp1")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "integrated_snn_res.0.4")



plot_genes_in_pseudotime(cds_subset, color_by =  "Pseudotime")

diff_test_res <- differentialGeneTest(HSMM_filtered,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene <- subset(diff_test_res, qval < 1e-2)
sig_gene_names <- rownames(subset(sig_gene, use_for_ordering == "TRUE"))
p1 = plot_pseudotime_heatmap(HSMM_filtered[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap = T)
p1$tree_row
table(diff_test_res$use_for_ordering)

#lung <- load_lung()
#plot_cell_trajectory(lung, color_by = "Time")
#画分向的heatmap
BEAM_res <- BEAM(HSMM_filtered, branch_point = 1, cores = 6)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
table(BEAM_res$qval < 1e-120)
branch_genes1 = subset(BEAM_res, qval < 1e-10)
branch_genes2 = subset(BEAM_res, pval < 1e-110)
branch_genes = rbind(branch_genes1,branch_genes2)
plot_genes_branched_heatmap(HSMM_filtered[row.names(branch_genes1),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F)

#提取monocle2 heatmap中的基因名,branch的和单向的是不一样的
plot_genes_branched_heatmap(HSMM_filtered[row.names(branch_genes1),],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = F)
pm = plot_genes_branched_heatmap(HSMM_filtered[row.names(branch_genes1),],
                                 branch_point = 1,
                                 num_clusters = 3,
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T) #导PDF会重叠2个图
gene_group=pm$annotation_row
gene_group$gene=rownames(gene_group)
write.csv(gene_group,"momnocle_heatmap_gene.csv")

#分向点图
point_genes <- row.names(subset(fData(HSMM),
                               gene_short_name %in% c("Egr1", "Fos", "Rpl26")))
plot_genes_branched_pseudotime(HSMM[point_genes,],
                               branch_point = 1,
                               color_by = "celltype3",
                               ncol = 1)+  
  scale_colour_manual(
    values = col4
  )
