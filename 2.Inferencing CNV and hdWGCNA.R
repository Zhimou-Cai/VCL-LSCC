##====================================================================================##
##===================Step01：Analysis of CNV in Epithelial cells======================##
##====================================================================================##
if(F){
  ###### infercnv
  library(Seurat)
  library(infercnv)
  library(scPRIT)
  options(future.globals.maxSize = 30 * 1024^3)
  setwd("~/project/VCL_LSCC/Step2_maEC/")

  scRNA <- readRDS("~/project/VCL_LSCC/Step1_QC&Define/sco.hny_anno.rds")

  scRNAVCL <- subset(scRNA, Group %in% c("VCL"))
  scRNAVCL <- subset(scRNAVCL, CellType %in% c("Epithelial cell",'T/NK cell','B cell','Mast cell','Endothelial cell'))
  
  ppdata <- infercnv.Prepare(scRNAVCL, celltype.key = "CellType", samleid.key = "orig.ident", 
                             tumorOrigin = "Epithelial cell", nRefcell = 10000)
  
  gene.pos <- read.table("~/database/infercnv/hg38_gencode_v27.txt", header = F, row.names = 1)
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = ppdata[["counts"]],
                                      annotations_file = ppdata[["annodata"]],
                                      gene_order_file = gene.pos,
                                      ref_group_names = ppdata[["refnames"]]) 
  saveRDS(infercnv_obj, "infercnv_pp.rds")
  
  infercnv_obj <- readRDS("infercnv_pp.rds")
  cnv.res <- infercnv::run(infercnv_obj,
                           cutoff=0.1,
                           out_dir="infercnv_out", 
                           cluster_by_groups=TRUE, 
                           leiden_resolution = 0.0001,
                           denoise=TRUE,
                           num_threads = 10,
                           HMM=FALSE)
  
  ### Analysis of CNV Results
  library(infercnv)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  
  infercnv_obj = readRDS("infercnv_VCL/run.final.infercnv_obj")
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- c(normal_loc$`T/NK cell`,normal_loc$`Endothelial cell`,normal_loc$`B cell`,normal_loc$`Mast cell`)
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc <- c(test_loc$VCL1,test_loc$VCL2,test_loc$VCL3,test_loc$VCL4,test_loc$VCL5,
                test_loc$VCL6,test_loc$VCL7,test_loc$VCL8,test_loc$VCL9,test_loc$VCL10)
  
  anno.df=data.frame(
    CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
    class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
  )
  
  geneFile <- read.table("~/database/infercnv/hg38_gencode_v27.txt",header = F,sep = "\t",stringsAsFactors = F)
  rownames(geneFile)=geneFile$V1
  usedgene.sort=intersect(geneFile$V1,rownames(expr))
  sma11_geneFile <- geneFile[usedgene.sort,]
  expr=expr[usedgene.sort,]
  
  #kmeans clustering
  set.seed(20250607)
  kmeans.result <- kmeans(t(expr), 6)
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB")
  kmeans_df_s=arrange(kmeans_df,kmeans_class)
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class)
  head(kmeans_df_s)
  
  #Define heatmap annotations and color matching
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=RColorBrewer::brewer.pal(6, "Set1")[1:6]
  names(color_v)=as.character(1:6)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))
  
  pdf("try1.pdf",width = 15,height = 10)
  ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
               col = colorRamp2(c(0.6,1,1.4), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sma11_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
               column_gap = unit(2, "mm"),
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",
                                           title_position = "leftcenter-rot",at=c(0.6,1,1.4),legend_height = unit(3, "cm")),
               top_annotation = top_anno,left_annotation = left_anno, #添加注释
               row_title = NULL,column_title = NULL)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  write.table(kmeans_df_s, file = "kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
  #Combined with graphs and tables, tumor cells were identified
  table(kmeans_df_s$kmeans_class,kmeans_df_s$class)
  epi_anno_new=kmeans_df_s
  epi_anno_new$type <- ifelse(
    #epi_anno_new$kmeans_class == "1",
    epi_anno_new$kmeans_class == "2" | epi_anno_new$kmeans_class == "4",
    "non-mailgnant",
    "mailgnant"
  )
  table(epi_anno_new$kmeans_class,epi_anno_new$type)
  write.table(epi_anno_new, file = "epi_anno_new.txt",quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
  epi_anno_new$kmeans_class=NULL
  epi_anno_new$type=factor(epi_anno_new$type,levels = c("non-mailgnant","mailgnant"))
  epi_anno_new$class=factor(epi_anno_new$class,levels = c("normal","test"))
  #epi_anno_new=epi_anno_new%>%arrange(type,class)
  epi_anno_new=epi_anno_new%>%arrange(type)
  
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  left_anno <- rowAnnotation(df=epi_anno_new,col=list(class=c("normal"="blue","test"="red"),
                                                      type=c("non-mailgnant"="#ddddda","mailgnant"="#D30404")))
  
  pdf("new_heatmap.pdf",width = 15,height = 8)
  ht = Heatmap(t(expr)[rownames(epi_anno_new),],
               col = colorRamp2(c(0.9,1,1.1), c("#377EB8", "#F0F0F0", "#E41A1C")), 
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sma11_geneFile$V2, paste("chr",1:22,sep = "")), 
               column_gap = unit(0, "mm"),
               row_split = c(rep("2",5609),rep("4",10956),rep("3",4224),rep("1",3111),rep("5",924),rep("6",1740)),
               row_gap = unit(0,"mm"),
               border_gp = gpar(col="black",lty=1,lwd=2),
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.9,1,1.1),legend_height = unit(3, "cm")),
               top_annotation = top_anno,left_annotation = left_anno,
               row_title = NULL,column_title = NULL)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  ####### Violin plot of CNV scores
  expr2=expr-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  kmeans_df_s$CB=rownames(kmeans_df_s)
  CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")
  
  CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
    scale_y_continuous(limits = c(0, 0.01),breaks = seq(0, 0.02, 0.005))+
    scale_fill_manual(values = color_v)+
    theme_bw()
  ggsave("CNV level.pdf",width = 12,height = 8,units = "cm")
  
  expr2=expr-1
  expr2=expr2 ^ 2
  CNV_score=as.data.frame(colMeans(expr2))
  colnames(CNV_score)="CNV_score"
  CNV_score$CB=rownames(CNV_score)
  epi_anno_new$CB=rownames(epi_anno_new)
  CNV_score=CNV_score%>%inner_join(epi_anno_new,by="CB")
  write.table(CNV_score, file = "CNV_score.txt",quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
  CNV_score%>%ggplot(aes(type,CNV_score))+geom_violin(aes(fill=type),color="NA")+
    scale_y_continuous(limits = c(0, 0.01),breaks = seq(0, 0.02, 0.005))+
    scale_fill_manual(values = c("#ddddda","#D30404"))+
    theme_bw()
  ggsave("CNV level_type.pdf",width = 12,height = 8,units = "cm")
  
  CNV_score <- read.table("./project/HPC_CT/inferCNV/CNV_score.txt")
  rownames(CNV_score) <- CNV_score$CB
  sco.CNV <- scRNAVCL[,CNV_score$CB]
  sco.CNV<-AddMetaData(sco.CNV,metadata = CNV_score)
  #sco.CNV$CNV_score <- CNV_score$CNV_score
  #sco.CNV$type <- CNV_score$type
  colnames(sco.CNV@meta.data)
  
  
  saveRDS(sco.CNV, "sco.scRNAVCL_CNV.rds")
  scRNA<-readRDS('infercnv_VCL/sco.scRNAVCL_CNV.rds')
  
}

##====================================================================================##
##===========================Step02：hdWGCNA==========================================##
##====================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(enrichR)
  library(GeneOverlap)
  
  seurat_obj <- readRDS("./hdWGCNA_object.rds")
  
  ModuleCorrelogram(seurat_obj)
  
  # get hMEs from seurat object
  MEs <- GetMEs(seurat_obj, harmonized=TRUE)
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module)
  mods <- mods[mods !='grey']
  
  # add hMEs to Seurat meta-data:
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  
  # plot with Seurat's DotPlot function
  p <- DotPlot(seurat_obj, features=mods, group.by ='label')
  
  # flip the x/y axes, rotate the axis labels, and change color scheme:
  p <- p +
    RotatedAxis() +coord_flip()+
    scale_color_gradientn(colours = c("#268987", "#F9FAFC","#BB0722"))
  
  # plot output
  p
  
  seurat_obj@misc[["Malig.Tumor"]][["wgcna_net"]][["TOMFiles"]] <- '/home/rstudio/project/VCL_LSCC/Step2_maEC/TOM/Epithelial cell_TOM.rda'
  ModuleNetworkPlot(
    seurat_obj,
    outdir='/home/rstudio/project/VCL_LSCC/Step2_maEC/TOM',# new folder name
    n_inner =10,# number of genes in inner ring
    n_outer =15,# number of genes in outer ring
    n_conns =Inf,# show all of the connections
    plot_size=c(10,10),# larger plotting area
    vertex.label.cex=2.5# font size
  )
  
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs =10,# number of hub genes to include for the UMAP embedding
    n_neighbors=15,# neighbors parameter for UMAP
    min_dist=0.1# min distance between points in UMAP space
  )
  
  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)
  
  # plot with ggplot
  ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
    geom_point(
      color=umap_df$color,# color each point by WGCNA module
      size=umap_df$kME*2# size of each point based on intramodular connectivity
    ) +
    umap_theme()
  
  options(future.globals.maxSize = 3 * 1024^3)
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1,# proportion of edges to sample (20% here)
    label_hubs=2,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
  )
  
  
  ######################
  #dir.create("14-hdWGCNA")
  #setwd("14-hdWGCNA")
  
  # enrichr databases
  dbs <- c('MSigDB_Hallmark_2020','GO_Biological_Process_2025')
  
  # 富集分析
  seurat_obj <- RunEnrichr(
    seurat_obj,
    dbs=dbs,
    max_genes =100# use max_genes = Inf to choose all genes
  )
  
  # 检索输出表
  enrich_df <- GetEnrichrTable(seurat_obj)
  
  # 查看结果
  head(enrich_df)
  
  # make GO term plots:
  EnrichrBarPlot(
    seurat_obj,
    outdir ="enrichr_plots",# name of output directory
    n_terms =10,# number of enriched terms to show (sometimes more are shown if there are ties)
    plot_size = c(5,7),# width, height of the output .pdfs
    logscale=TRUE# do you want to show the enrichment as a log scale?
  )
  
  # enrichr dotplot
  EnrichrDotPlot(
    seurat_obj,
    mods ="all",# use all modules (default)
    database ="MSigDB_Hallmark_2020",# this must match one of the dbs used previously
    n_terms=2,# number of terms per module
    term_size=8,# font size for the terms
    p_adj =FALSE# show the p-val or adjusted p-val?
  ) + scale_color_stepsn(colors=rev(viridis::magma(256)))
}