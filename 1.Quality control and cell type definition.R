##====================================================================================##
##===============================Step01：QC===========================================##
##====================================================================================##
if(F){
  #Attention, please! The scPRIT package relies on specific R and Python environments and database resources, 
  #which the developers integrate into the docker image analysis environment. 
  #Specific access method please see the following link: https://mp.weixin.qq.com/s/sR12NVP179xZK5VEgRSxmA
  
  library(Seurat)
  library(data.table)
  library(parallel)
  library(tidyverse)
  library(patchwork)
  library(scPRIT)
  rm(list = ls())  
  
  setwd("~/project/VCL_LSCC/")
  if(!file.exists("QC")) dir.create("QC")
  if(!file.exists("Data")) dir.create("Data")
  
  mat.list <- list(
    VCP1 = Read10X("SourceData/VCP1/"),
    VCP2 = Read10X("SourceData/VCP2/"),
    VCP3 = Read10X("SourceData/VCP3/"),
    VCP4 = Read10X("SourceData/VCP4/"),
    VCP5 = Read10X("SourceData/VCP5/"),
    VCP6 = Read10X("SourceData/VCP6/"),
    VCP7 = Read10X("SourceData/VCP7/"),
    VCP8 = Read10X("SourceData/VCP8//"),
    VCP9 = Read10X("SourceData/VCP9/"),
    VCP10 = Read10X("SourceData/VCP10/"),
    VCL1 = Read10X("SourceData/VCL1//"),
    VCL2 = Read10X("SourceData/VCL2/"),
    VCL3 = Read10X("SourceData/VCL3/"),
    VCL4 = Read10X("SourceData/VCL4/"),
    VCL5 = Read10X("SourceData/VCL5/"),
    VCL6 = Read10X("SourceData/VCL6/"),
    VCL7 = Read10X("SourceData/VCL7/"),
    VCL8 = Read10X("SourceData/VCL8//"),
    VCL9 = Read10X("SourceData/VCL9/"),
    VCL10 = Read10X("SourceData/VCL10/"),
    LSCC1 = Read10X("SourceData/LSCC1/"),
    LSCC2 = Read10X("SourceData/LSCC2/"),
    LSCC3 = Read10X("SourceData/LSCC3/"),
    LSCC4 = Read10X("SourceData/LSCC4/"),
    LSCC5 = Read10X("SourceData/LSCC5/"),
    LSCC6 = Read10X("SourceData/LSCC6/"),
    LSCC7 = Read10X("SourceData/LSCC7/"),
    LSCC8 = Read10X("SourceData/LSCC8//"),
    LSCC9 = Read10X("SourceData/LSCC9/"),
    LSCC10 = Read10X("SourceData/LSCC10/"),
    LSCCP1 = Read10X("SourceData/LSCCP1/"),
    LSCCP2 = Read10X("SourceData/LSCCP2/"),
    LSCCP3 = Read10X("SourceData/LSCCP3/"),
    LSCCP4 = Read10X("SourceData/LSCCP4/"),
    LSCCP5 = Read10X("SourceData/LSCCP5/"),
    LSCCP6 = Read10X("SourceData/LSCCP6/"),
    LSCCP7 = Read10X("SourceData/LSCCP7/"),
    LSCCP8 = Read10X("SourceData/LSCCP8//"),
    LSCCP9 = Read10X("SourceData/LSCCP9/"),
    LSCCP10 = Read10X("SourceData/LSCCP10/"))
  
  sn = names(mat.list)
  scRNAlist <- mclapply(sn, function(i){
    sco <- CreateSeuratObject(mat.list[[i]], project = i, min.cells=3, min.features = 200)
    sco <- RenameCells(sco, add.cell.id = i)
    
    sco[["percent.mt"]] <- PercentageFeatureSet(sco, pattern = "^MT-")
    
    sco[["percent.rb"]] <- PercentageFeatureSet(sco, pattern = "^RP[LS]")
    
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(sco))
    sco[["percent.HB"]]<-PercentageFeatureSet(sco, features=HB.genes) 
    return(sco)}, mc.cores = length(sn))
  names(scRNAlist) <- sn
  rm(mat.list, sn)
  saveRDS(scRNAlist, file = "Data/scRNAlist0.rds")
  
  
  ### The theoretical proportion of doublets was calculated from the Poisson distribution
  scRNAlist<-readRDS("Data/scRNAlist0.rds")
  expected_doublet_rate = sapply(scRNAlist, calcDBRate)
  save(expected_doublet_rate, file = "QC/expected_doublet_rate.rda")
  
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
  scRNA <- JoinLayers(scRNA)
  featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
  group = "orig.ident"
  
  # Data before quality control
  plots = list()
  for(i in seq_along(featrures)){
    p = VlnPlot(scRNA, group.by=group, pt.size = 0, features = featrures[i], layer = "counts")
    plots[[i]] = p + NoLegend() + theme(axis.title.x=element_blank())
  }
  violin <- wrap_plots(plots = plots, nrow=2)    
  ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 25, height = 8) 
  
  df1 <- as.data.frame(table(scRNA$orig.ident))
  colnames(df1) <- c("Sample", "nCells")
  #df1$Sample <- paste0(df1$Sample, "_before")
  df1$Group <- "before"
  
  # Quality control indicators
  quantile(scRNA$nFeature_RNA, seq(0.01, 0.2, 0.01))
  quantile(scRNA$nFeature_RNA, seq(0.90, 1, 0.002))
  quantile(scRNA$nCount_RNA, seq(0.90, 1, 0.002))
  quantile(scRNA$percent.mt, seq(0.9, 1, 0.01))
  quantile(scRNA$percent.HB, seq(0.9, 1, 0.01))
  plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 5000)
  plots[[2]] + geom_hline(yintercept = 25000)
  plots[[3]] + geom_hline(yintercept = 15)
  plots[[5]] + geom_hline(yintercept = 1)
  
  minGene=500
  maxGene=5000
  maxUMI=25000
  pctMT=15
  pctHB=1
  
  # Quality control
  scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                    nCount_RNA < maxUMI & percent.mt < pctMT & percent.HB < pctHB)
  
  ### Cell cycle score
  g2m_genes <- cc.genes$g2m.genes
  g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(scRNA))
  s_genes <- cc.genes$s.genes    
  s_genes <- CaseMatch(search=s_genes, match=rownames(scRNA))
  sco.tmp <- NormalizeData(scRNA)
  sco.tmp <- CellCycleScoring(sco.tmp, g2m.features=g2m_genes, s.features=s_genes)
  scRNA <- AddMetaData(scRNA, metadata = sco.tmp@meta.data[,c("S.Score","G2M.Score","Phase")])
  
  ### Predicted cell types with doublets
  ## SingleR
  refdata <- readRDS("~/database/SingleR/Human_all.rds")
  sco.list <- SplitObject(scRNA, split.by = "orig.ident")
  for(i in seq_along(sco.list)) {
    sco.list[[i]] <- NormalizeData(sco.list[[i]]) %>% FindVariableFeatures() %>% ScaleData()
    sco.list[[i]] <- RunPCA(sco.list[[i]], verbose = FALSE)
    sco.list[[i]] <- RunUMAP(sco.list[[i]], dims = 1:15)
    sco.list[[i]] <- FindNeighbors(sco.list[[i]], dims = 1:15)
    sco.list[[i]] <- FindClusters(sco.list[[i]], resolution = 2)
    sco.list[[i]] <- runSingleR(sco.list[[i]], refdata = refdata)
  }
  
  ## DoubletFinder
  load("Step1_QC&Define/QC/expected_doublet_rate.rda")
  for(i in names(sco.list)){
    sco <- sco.list[[i]]
    db.rate <- as.numeric(expected_doublet_rate[i])
    sco.list[[i]] <- runDoubletFinder(sco, celltype.key = "SingleR", doublet_rate = db.rate)
  }
  
  scRNA <- merge(sco.list[[1]], sco.list[2:length(sco.list)])
  saveRDS(scRNA, file = "Data/sco.qc.rds")
  
  ##Removing doubles
  scRNA1<-readRDS("Data/sco.qc.rds")
  scRNA1 <- subset(scRNA1, DF.doublets != "Doublet")
  
  # Data after quality control
  df2 <- as.data.frame(table(scRNA1$orig.ident))
  colnames(df2) <- c("Sample", "nCells")
  #df2$Sample <- paste0(df2$Sample, "_after")
  df2$Group <- "after"
  
  df <- rbind(df1, df2)
  df$Group <- relevel(factor(df$Group), ref = "before")
  p <- ggplot(df, aes(x=Sample, y=nCells, fill=Group)) + 
    geom_bar(stat= "identity", position = "dodge", width = 0.6) + 
    geom_text(aes(label=nCells, y=nCells+200), size=4) + 
    ggsci::scale_fill_jco() +
    theme_bw()
  p <- p + RotatedAxis()
  ggsave("QC/nCells.pdf", plot = p, width = 18, height = 8)
  
  featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
  group = "orig.ident"
  plots = list()
  for(i in seq_along(featrures)){
    p = VlnPlot(scRNA1, group.by=group, pt.size = 0, features = featrures[i], layer = "counts")
    plots[[i]] = p + NoLegend() + theme(axis.title.x=element_blank())
  }
  violin <- wrap_plots(plots = plots, nrow=2)     
  ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 25, height = 8)
  
  saveRDS(scRNA1, file = "Data/sco.qc.rds")
}


##====================================================================================##
##=======step02: Data integration and dimensionality reduction clustering=============##
##====================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(parallel)
  library(future)
  library(scPRIT)
  rm(list = ls())
  setwd("~/project/VCL_LSCC/Step1_QC&Define/")
  options(future.globals.maxSize = 60 * 1024^3)
  
  ### Data preprocessing
  scRNA <- readRDS("~/project/VCL_LSCC/Data/sco.qc.rds")
  scRNA <- NormalizeData(scRNA)
  scRNA <- FindVariableFeatures(scRNA)
  scRNA <- ScaleData(scRNA, vars.to.regress = c("S.Score","G2M.Score"))
  scRNA <- RunPCA(scRNA, verbose = F)
  ElbowPlot(scRNA, ndims = 50)
  saveRDS(scRNA, file = "Step1_QC&Define/sco.log_pp.rds")

  ## Harmony整合 
  rm(list = ls())
  scRNA <- readRDS("sco.log_pp.rds")
  scRNA <- IntegrateLayers(
    scRNA, method = HarmonyIntegration, lambda = 1,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE)
  # 整合结果降维聚类
  scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:20)
  scRNA <- FindClusters(scRNA, resolution = 1)
  scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20)
  scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:20)
  
  p1 <- FeaturePlot(scRNA, reduction = "umap", features = "DF.scores")
  p2 <- DimPlot(scRNA, reduction = "umap", group.by = "DF.doublets") + 
    scale_color_manual(values = c("red","gray"))
  p3 <- FeaturePlot(scRNA, reduction = "umap", features = "scr_score")
  p4 <- DimPlot(scRNA, reduction = "umap", group.by = "scr_pred") + 
    scale_color_manual(values = c("gray","red"))
  p <- (p1|p2)/(p3|p4)
  ggsave("hny_Doublets.pdf", p, width = 12, height = 9)
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
  p2 <- DimPlot(scRNA, reduction = "umap", group.by = "SingleR", label = T) + NoLegend()
  p <- p1 + p2
  ggsave("hny_SingleR_Cluster.pdf", p, width = 15, height = 7)
  
  p3 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", label = T) + NoLegend()
  p4 <- DimPlot(scRNA, reduction = "tsne", group.by = "SingleR", label = T) + NoLegend()
  p <- p3 + p4
  ggsave("hny_SingleR_Cluster_tsne.pdf", p, width = 15, height = 7)
  
  p1 <- FeaturePlot(scRNA, reduction = "tsne", features = "DF.scores")
  p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "DF.doublets") + 
    scale_color_manual(values = c("red","gray"))
  p3 <- FeaturePlot(scRNA, reduction = "tsne", features = "scr_score")
  p4 <- DimPlot(scRNA, reduction = "tsne", group.by = "scr_pred") + 
    scale_color_manual(values = c("gray","red"))
  p <- (p1|p2)/(p3|p4)
  ggsave("hny_Doublets_tsne.pdf", p, width = 12, height = 9)
  
  saveRDS(scRNA, file = "sco.hny_dr.rds")
  
}


##====================================================================================##
##============================step03: Cell Type Annotation============================##
##====================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(parallel)
  library(scPRIT)
  setwd("~/project/VCL_LSCC/Step1_QC&Define/")
  options(future.globals.maxSize = 30 * 1024^3)
  
  scRNA <- readRDS("~/project/VCL_LSCC/Step1_QC&Define/sco.hny_dr.rds")
  scRNA <- JoinLayers(scRNA)
  
  ### Cluster differentially expressed genes
  dir.create("ClusterDEG")
  clusterDEG <- FindAllMarkers(scRNA, only.pos = T, logfc.threshold = 0.5, slot = "data")
  write.csv(clusterDEG, file = "ClusterDEG/clusterDEG.csv", row.names = F)
  saveRDS(clusterDEG, file = "ClusterDEG/clusterDEG.rds")
  # Screening top markers
  topDEGs <- filterAllMarkers(clusterDEG, topN = 50, log2FC.T = 0.5, 
                              p_val_adj.T = 0.05, rm.MT = T, rm.RB = T)
  write.csv(topDEGs, file = "ClusterDEG/clusterDEG_top50.csv", row.names = F)
  saveRDS(topDEGs, file = "ClusterDEG/clusterDEG_top50.rds")
  head(topDEGs)
  
  ## SingleR annotation
  refdata <- readRDS("~/database/SingleR/Human_all.rds")
  scRNA <- runSingleR(scRNA, refdata = refdata)
  DimPlot(scRNA, group.by = "SingleR", label = T) + NoLegend()

  saveRDS(scRNA, file = "sco.hny_dr.rds")
  
  ## The prediction results are plotted
  p1 <- DimPlot(scRNA, group.by = "seurat_clusters", label = T) + NoLegend()
  p2 <- DimPlot(scRNA, group.by = "SingleR", label = T) + NoLegend()
  p <- p1|p2
  ggsave("celltype_reference.pdf", p, width = 18, height = 6)
  

  ### Manually defined cells
  scRNA$CellType <- recode(scRNA$seurat_clusters,
                           "0" = "T/NK cell",
                           "1" = "Fibroblast",
                           "2" = "T/NK cell",
                           "3" = "T/NK cell",
                           "4" = "Fibroblast",
                           "5" = "T/NK cell",
                           "6" = "Myeloid cell",
                           "7" = "Endothelial cell",
                           "8" = "B cell",
                           "9" = "Myeloid cell",
                           "10" = "Epithelial cell",
                           "11" = "Myeloid cell",
                           "12" = "Fibroblast",
                           "13" = "Epithelial cell",
                           "14" = "Fibroblast",
                           "15" = "Plasma cell",
                           "16" = "Myeloid cell",
                           "17" = "Epithelial cell",
                           "18" = "Epithelial cell",
                           "19" = "Fibroblast", 
                           "20" = "Mast cell",
                           "21" = "Epithelial cell", 
                           "22" = "Fibroblast",
                           "23" = "Endothelial cell",
                           "24" = "Epithelial cell",
                           "25" = "Fibroblast",
                           "26" = "Myeloid cell")
  
  saveRDS(scRNA, file = "~/project/VCL_LSCC/Step1_QC&Define/sco.hny_anno.rds")  
  
  scRNA$CellType <-  factor(scRNA@meta.data$CellType,levels=c("Epithelial cell","Myeloid cell","Fibroblast","Endothelial cell",
                                                              "B cell","Plasma cell","T/NK cell","Mast cell"))
  scRNA$Group <-  factor(scRNA@meta.data$Group,levels=c("LSCC","VCL","VCP","LSCCP"))
  
  mycol <- c("#E41A1C", "#FF7F00","#984EA3","#FFD121","#F781BF","#377EB8","#4DAF4A","#A65628")
  groupcol <- c("#1F78B4",'#BC80BD',"#D20A13","#ebb2d2")
  
  p1 <- DimPlot(scRNA, group.by = "Group", label = F,raster = F,cols = groupcol)+
    labs(x = "UMAP1", y = "UMAP2",title = "Treatment") +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p2 <- DimPlot(scRNA, group.by = "CellType", label = F,raster = F,cols = mycol)+
    labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p <- p1|p2
  ggsave("Celltype_UMAP.pdf", p, width = 16, height = 6.5)
  
  ###Marker gene heatmap & dot map
  Idents(scRNA)<-"CellType"
  clusterDEG <- FindAllMarkers(scRNA, only.pos = T, logfc.threshold = 0.5, slot = "data")
  # top markers
  topDEGs <- filterAllMarkers(clusterDEG, topN = 15, log2FC.T = 0.5, 
                              p_val_adj.T = 0.05, rm.MT = T, rm.RB = T)
  
  gene <- c('KRT5','KRT14','KRT13',   # Epithelial
            'CD14', 'CCL3','ITGAX', # Myeloid
            'DCN','COL1A2', 'COL3A1',  # Fibroblast
            'VWF','PLVAP', 'CLDN5',  # Endothelial
            'BANK1','MS4A1','CD79A',# B cell
            'MZB1', 'XBP1','JCHAIN',  # Plasma
            'CD2','CD3D', 'NKG7', # T cell
            'TPSAB1','CPA3', 'KIT')  # Mast cell 
  
  pdf(file="DotPlot.pdf", w=10, h=4.5, family = "sans")
  DotPlot(scRNA, 
          features = gene,
          group.by = "CellType") +  # 使用排序后的细胞类型
    scale_color_gradientn(colours = c('#336699',"#87CEEB","#FFA07A","#D41F26")) +
    labs(x = "Markers", y = "Cell type") +
    theme(legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          axis.title.x =element_text(size=15, face="bold"), 
          axis.title.y=element_text(size=15, face="bold"),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, color = "black", size = 13.5, hjust = 1, vjust=1),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_line(colour = "grey95", linewidth = 0.2), 
          panel.grid.minor = element_line(colour = "grey95", linewidth = 0.1)) +  
    geom_point(aes(size = pct.exp, color = avg.exp.scaled), 
               shape = 21, colour = "black", stroke = 0.5) 
  dev.off()
  
  saveRDS(scRNA, file = "~/project/VCL_LSCC/Step1_QC&Define/sco.hny_anno.rds")  
  
  ###Proportion statistics of cell types
  library(ggpubr)
  library(RColorBrewer)
  library(tidyr)
  library(Seurat)
  library(patchwork)
  scRNA <- readRDS("sco.hny_anno.rds")

  levels(factor(scRNA$CellType))
  mycol <- c("#E41A1C", "#FF7F00","#984EA3","#FFD121","#F781BF","#377EB8","#4DAF4A","#A65628")
  
  bar.df=scRNA@meta.data
  text.df1=as.data.frame(table(bar.df$CellType))
  p1=bar.df%>%ggplot(aes(x=CellType))+geom_bar(aes(fill=Group),position = "fill")+
    scale_x_discrete("")+
    scale_y_continuous("cell ratio",expand = c(0.02,0),labels = scales::label_percent())+
    scale_fill_manual(values = c('#93CBAE','#F3AE63','#73558B','#AC1B1F'))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      legend.title = element_blank())
  ggsave("cell_group.bar1.pdf",width = 12,height = 18,units = "cm")
  
  p2=bar.df%>%ggplot(aes(x=Group))+geom_bar(aes(fill=CellType),position = "fill")+
    scale_x_discrete("")+
    scale_y_continuous("cell ratio",expand = c(0.02,0),labels = scales::label_percent())+
    scale_fill_manual(values = mycol)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      legend.title = element_blank())
  ggsave("cell_group.bar2.pdf",width = 12,height = 18,units = "cm")
  
  p1+p2+plot_layout(ncol = 2)
  ggsave("PatientCelltype_bar.pdf",height = 9,width = 12)
  
}
