##===============================================================================================================##
##=========================Step01：Data preprocessing for interaction analysis===================================##
##===============================================================================================================##
if(F){
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(patchwork)
  library(reshape2)
  rm(list = ls())
  setwd("~/project/VCL_LSCC/Step5_CellChat/")
  
  ### Organize the data
  scRNA <- readRDS("~/project/VCL_LSCC/Step1_QC&Define/sco.hny_anno.rds")
  scRNA <- subset(scRNA, CellType != "T/NK cell")
  scRNA <- subset(scRNA, CellType != "Myeloid cell")
  counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
  colnames(scRNA@meta.data)
  meta <- scRNA@meta.data[,c(4:10,15)]
  scRNA <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco.MM <- readRDS("~/project/VCL_LSCC/Step3.2_Myeloid/sco.scRNAMM_hny_anno.rds")
  sco.MM$CellType <- recode(sco.MM$celltype,"Monocyte" = "Monocyte","ISG15_macrophage" = "Macrophage",
                            "C1QC_macrophage" = "Macrophage","LYVE1_macrophage" = "Macrophage","SPP1_macrophage" = "Macrophage")
  colnames(sco.MM@meta.data)
  counts <- GetAssayData(sco.MM, assay = "RNA", layer = "counts")
  meta <- sco.MM@meta.data[,c(4:10,22)]
  sco.MM <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco.DC <- readRDS("~/project/VCL_LSCC/Step3.2_Myeloid/sco.scRNADC_hny_anno.rds")
  sco.DC$CellType <- recode(sco.DC$celltype,"cDC1" = "DC","cDC2" = "DC","pDC" = "DC","LAMP3_DC" = "DC")
  colnames(sco.DC@meta.data)
  counts <- GetAssayData(sco.DC, assay = "RNA", layer = "counts")
  meta <- sco.DC@meta.data[,c(4:10,21)]
  sco.DC <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco.tcell <- readRDS("~/project/VCL_LSCC/Step4.1_Tcell/sco.TNK_hny_anno.rds")
  sco.tcell$CellType <- recode(sco.tcell$celltype,"CD4_Tnaive" = "CD4 T cell","CD4_Tfh" = "CD4 T cell","CD4_Th17" = "CD4 T cell","CD4_Treg" = "CD4 T cell",
                               "CD8_Teff" = "CD8 T cell","CD8_Tm" = "CD8 T cell","CD8_Tex" = "CD8 T cell","CD8_ISG+ T" = "CD8 T cell",
                               "NKT cell" = "NKT cell","NK cell" = "NK cell")
  colnames(sco.tcell@meta.data)
  counts <- GetAssayData(sco.tcell, assay = "RNA", layer = "counts")
  meta <- sco.tcell@meta.data[,c(4:10,27)]
  sco.tcell <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco <- merge(scRNA,list(sco.MM,sco.DC,sco.tcell))
  sco <- JoinLayers(sco)
  sco$CellType <- as.character(sco$CellType)
  sco <- NormalizeData(sco)
  unique(sco$CellType)
  saveRDS(sco, "scRNA_Cellchat.rds")
  
  scRNA<-readRDS("scRNA_Cellchat.rds")
}


##===============================================================================================================##
##=========================Step02：Strength of EC/EpC interaction between each cell type=========================##
##===============================================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(CellChat)
  library(scPRIT)
  library(reshape2)
  library(reshape)
  setwd("~/project/VCL_LSCC/Step5_CellChat/")
  scRNA <- readRDS("~/project/VCL_LSCC/Step5_CellChat/scRNA_Cellchat.rds")
  scRNA <- subset(scRNA, Group %in% "VCL")
  scRNA <- subset(scRNA, CellType != "Epithelial cell")
  counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
  colnames(scRNA@meta.data)
  meta <- scRNA@meta.data[,c(4:11)]
  scRNA <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco.EpC <- readRDS("~/project/VCL_LSCC/Step2_maEC/sco.maEC_hny_dr.rds")
  sco.EpC <- subset(sco.EpC, Group %in% "VCL")
  colnames(sco.EpC@meta.data)
  counts <- GetAssayData(sco.EpC, assay = "RNA", layer = "counts")
  meta <- sco.EpC@meta.data[,c(4:10,15)]
  sco.EpC <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco <- merge(scRNA,list(sco.EpC))
  sco <- JoinLayers(sco)
  sco$CellType <- as.character(sco$CellType)
  sco <- NormalizeData(sco)
  unique(sco$CellType)
  saveRDS(sco, "VCL_Cellchat.rds")
  
  ###cellchat
  sco<-readRDS("VCL/VCL_Cellchat.rds")
  runCellChat(sco, celltype.key = "CellType", output = "VCL", species = "human", min.cells = 5)
  
  setwd("~/project/VCL_LSCC/Step5_CellChat/VCL/")
  scRNA_cellchat <- readRDS('VCL/cellchat_obj.rds')
  
  NetWeightAll <- scRNA_cellchat@net$weight
  NetWeightAll <- as.data.frame(NetWeightAll)
  
  NetWeightAll <- bind_rows(NetWeightAll)
  NetWeightAll_melt <- melt(NetWeightAll)
  NetWeightAll_melt$cluster <- rep(colnames(NetWeightAll)[1:13],times = 1)
  NetWeightAll_melt$cellchat <- paste(NetWeightAll_melt$cluster,'to',NetWeightAll_melt$variable,sep = ' ')
  
  unique(NetWeightAll_melt$cellchat)
  NetWeightAll_melt$target <- rep("EC",times = 1)
  
  NetWeightAll_meltF <- NetWeightAll_melt %>% filter(variable == 'Endothelial cell')#挑选作用于Epithelial cell的分数
  NetWeightAll_meltF$cellchat2 <- paste(NetWeightAll_meltF$cellchat,NetWeightAll_meltF$target)
  NetWeightAll_meltF$cellchat2 <- stringr::str_remove(NetWeightAll_meltF$cellchat2,"Endothelial cell")
  #NetWeightAll_meltF$cellchat2 <- paste(NetWeightAll_meltF$cellchat2,NetWeightAll_meltF$variable)
  
  #调整顺序
  NetWeightAll_meltF <- NetWeightAll_meltF %>% arrange(desc(value))
  
  
  mycol <- c("#B15928","#1F78B4","#68228B","#223D6C","#FFD121","#D20A13",
             "#088247","#FF7F00","#FB9A99","#7A142C","#FDBF6F","#EE82EE","grey")
  
  
  P <- ggplot(NetWeightAll_meltF,aes(x=value,y=cellchat2))+
    ylab('')+
    xlab('weight')+
    #ggtitle(sub)+
    geom_segment(aes(yend=cellchat2),xend=0,colour='grey50')+  ###绘制以数据点为端点的线段
    geom_point(size=3,aes(colour=cellchat))+   ###此处我们将以正负相关(postive  negative)映射其颜色
    scale_color_manual(values = mycol)+ ###颜色加深  
    # scale_y_discrete(limits= rev(c('NT','Pre','E','A','R')))+
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_y_discrete(limits= rev(NetWeightAll_meltF$cellchat2))+
    NoLegend()###删除网格线
  dev.off()
  
  ggsave("CelltoEC_VCL_Lolliplot.pdf",width = 10,height = 8.5,units = "cm")
}

##===============================================================================================================##
##=========================Step03：Strength of interaction between Fib and EC/EpC================================##
##===============================================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(CellChat)
  library(scPRIT)
  library(reshape2)
  library(reshape)
  library(dplyr)
  
  setwd("~/project/VCL_LSCC/Step5_CellChat/VCL/Fib_EC/")
  sco.Fib <- readRDS("~/project/VCL_LSCC/Step3.1_Fibro/sco.Fibro_hny_anno.rds")
  sco.Fib <- subset(sco.Fib, Group %in% "LSCC")
  colnames(sco.Fib@meta.data)
  counts <- GetAssayData(sco.Fib, assay = "RNA", layer = "counts")
  meta <- sco.Fib@meta.data[,c(4:10,17)]
  sco.Fib <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco.EC <- readRDS("~/project/VCL_LSCC/Step3.3_Endo/sco.EC_hny_anno.rds")
  sco.EC <- subset(sco.EC, Group %in% "LSCC")
  colnames(sco.EC@meta.data)
  counts <- GetAssayData(sco.EC, assay = "RNA", layer = "counts")
  meta <- sco.EC@meta.data[,c(4:10,17)]
  sco.EC <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco <- merge(sco.Fib,list(sco.EC))
  sco <- JoinLayers(sco)
  sco$celltype <- as.character(sco$celltype)
  sco <- NormalizeData(sco)
  unique(sco$celltype)
  saveRDS(sco, "Fib_EC_Cellchat.rds")
  
  ######Cellchat
  sco<-readRDS("Fib_EC_Cellchat.rds")
  runCellChat(sco, celltype.key = "celltype", output = "./", species = "human", min.cells = 5)
  
  scRNA_EC_cellchat <- readRDS('cellchat_obj.rds')
  
  NetWeightAll <- scRNA_EC_cellchat@net$weight
  NetWeightAll <- as.data.frame(NetWeightAll)
  
  NetWeightAll <- bind_rows(NetWeightAll)
  NetWeightAll_melt <- melt(NetWeightAll)
  NetWeightAll_melt$cluster <- rep(colnames(NetWeightAll)[1:10],times = 1)
  NetWeightAll_melt$cellchat <- paste(NetWeightAll_melt$cluster,'to',NetWeightAll_melt$variable,sep = ' ')
  
  unique(NetWeightAll_melt$cellchat)
  
  NetWeightAll_meltF <- NetWeightAll_melt %>% filter(variable %in% c('Venous EC','Arterial EC','Proliferative EC',"Lymphatic EC"))
  NetWeightAll_meltF <- NetWeightAll_meltF %>% filter(cluster %in% c('CFD_fibroblast','POSTN_fibroblast','APCDD1_fibroblast',
                                                                     'CD74_fibroblast','ACTA2_myofibroblast','Proliferating fibroblast'))
  
  NetWeightAll_meltF <- NetWeightAll_meltF %>% arrange(desc(value))
  variable_order <- c('Proliferative EC','Venous EC',"Lymphatic EC",'Arterial EC')
  NetWeightAll_meltF$variable <- factor(NetWeightAll_meltF$variable, levels = variable_order)
  NetWeightAll_meltF <- NetWeightAll_meltF %>%arrange(variable)
  
  
  mycol <- c("#D20A13","#1F78B4",'#e184ca',"#FF7F00","#68228B","#223D6C","#FFD121","#B15928",
             "#088247","#FB9A99","#7A142C","#FDBF6F","#EE82EE","grey")
  
  
  P <- ggplot(NetWeightAll_meltF,aes(x=value,y=cellchat))+
    ylab('')+
    xlab('weight')+
    geom_segment(aes(yend=cellchat),xend=0,colour='grey50')+  
    geom_point(size=3,aes(colour=variable))+   
    scale_color_manual(values = mycol)+ 
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_y_discrete(limits= rev(NetWeightAll_meltF$cellchat))+
    NoLegend()
  dev.off()
  
  ggsave("FibtoEC_VCL_Lolliplot.pdf",width = 12,height = 9.5,units = "cm")
  
  ######CellphoneDB
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(scPRIT)
  
  setwd("~/project/VCL_LSCC/Step5_CellChat/LSCC/Fib_EC/")
  sco<-readRDS("Fib_EC_Cellchat.rds")
  runCellPhoneDB(sco, group.by = "celltype", output = "./", ncore = 60)
  
  setwd("~/project/VCL_LSCC/Step5_CellChat/VCL/Fib_EC/cellphonedb/")
  pvals <- read.delim("statistical_analysis_pvalues.txt", check.names = FALSE)
  means <- read.delim("statistical_analysis_means.txt", check.names = FALSE)
  decon = read.delim("statistical_analysis_deconvoluted.txt", check.names = FALSE)
  
  ### Statistics of the number of cell significant interactions
  res <- cpdb.CountHeatmap(cpdb_pvals = pvals, display_numbers = T, return_data = T)
  ggsave("cci_count_heatmap.pdf", res$plot, width = 8, height = 7.5)
  
  ### Cell interaction weight statistics
  res <- cpdb.WeightHeatmap(cpdb_means = means, cpdb_pvals = pvals, display_numbers = T,
                            color = c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                      "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A",
                                      "#FD705F", "#E12234","#840729"),
                            significant_only = T, return_data = T)
  ggsave("cci_weight_heatmap.pdf", res$plot, width = 8, height = 7.5)

  pdf("cci_weight_graph_POSTN_fibroblast.pdf", width = 8, height = 8)
  cpdb.WeightGraph(cpdb_means = means, cpdb_pvals = pvals, significant_only = T, 
                   source.use = c('POSTN_fibroblast'), edge.label = T)
  dev.off()
  
}

##===============================================================================================================##
##=========================Step04：Interaction RL of each Fib on EpC/Venous/Lymphatic EC=========================##
##===============================================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(CellChat)
  library(scPRIT)
  library(reshape2)
  library(reshape)
  library(jjAnno)
  
  setwd("~/project/VCL_LSCC/Step5_CellChat/Fib_EC_detail/")
  sco.Fib <- readRDS("~/project/VCL_LSCC/Step3.1_Fibro/sco.Fibro_hny_anno.rds")
  sco.Fib <- subset(sco.Fib, celltype %in% c('POSTN_fibroblast','APCDD1_fibroblast'))
  colnames(sco.Fib@meta.data)
  counts <- GetAssayData(sco.Fib, assay = "RNA", layer = "counts")
  meta <- sco.Fib@meta.data[,c(4:10,17)]
  sco.Fib <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco.EC <- readRDS("~/project/VCL_LSCC/Step3.3_Endo/sco.EC_hny_anno.rds")
  sco.EC <- subset(sco.EC, celltype %in% c('Venous EC','Lymphatic EC'))
  colnames(sco.EC@meta.data)
  counts <- GetAssayData(sco.EC, assay = "RNA", layer = "counts")
  meta <- sco.EC@meta.data[,c(4:10,17)]
  sco.EC <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  sco <- merge(sco.Fib,list(sco.EC))
  sco <- JoinLayers(sco)
  sco$celltype <- as.character(sco$celltype)
  sco <- NormalizeData(sco)
  unique(sco$celltype)
  saveRDS(sco, "Fib_EC_detail_Cellchat.rds")
  
  ###cellchat
  sco<-readRDS("Fib_EC_detail_Cellchat.rds")
  sco.list <- SplitObject(sco, split.by = "Group")
  sco.list <- lapply(sco.list, FUN = NormalizeData)
  sco <- sco.list[['LSCCP']]
  runCellChat(sco, celltype.key = "celltype", output = "LSCCP", species = "human", min.cells = 5)
  sco <- sco.list[['VCP']]
  runCellChat(sco, celltype.key = "celltype", output = "VCP", species = "human", min.cells = 5)
  sco <- sco.list[['VCL']]
  runCellChat(sco, celltype.key = "celltype", output = "VCL", species = "human", min.cells = 5)
  sco <- sco.list[['LSCC']]
  runCellChat(sco, celltype.key = "celltype", output = "LSCC", species = "human", min.cells = 5)
  
  LSCCP <- readRDS("~/project/VCL_LSCC/Step5_CellChat/Fib_EC_detail/LSCCP/cellchat_obj.rds")
  LCP <- readRDS("~/project/VCL_LSCC/Step5_CellChat/Fib_EC_detail/VCP/cellchat_obj.rds")
  LCV <- readRDS("~/project/VCL_LSCC/Step5_CellChat/Fib_EC_detail/VCL/cellchat_obj.rds")
  LSCC <- readRDS("~/project/VCL_LSCC/Step5_CellChat/Fib_EC_detail/LSCC/cellchat_obj.rds")
  
  cellchat_list <- list('LSCCP' = LSCCP, 'LCP' = LCP, 'LCV' = LCV, 'LSCC' = LSCC)
  uterus_cellchat <- mergeCellChat(cellchat_list,add.names = names(cellchat_list))
  
  netVisual_bubble(uterus_cellchat, 
                   sources.use = c('POSTN_fibroblast','APCDD1_fibroblast'), 
                   targets.use = c('Venous EC','Lymphatic EC'), 
                   comparison = c(1, 2, 3,4))
  

  pairLR.use <- as.data.frame(c("TNC_ITGA9_ITGB1",
                                "THBS2_CD47",
                                "THBS2_CD36",
                                "LAMC3_ITGA9_ITGB1",
                                "FN1_ITGA5_ITGB1",
                                'LAMC1_ITGA6_ITGB1'))
  colnames(pairLR.use) <- 'interaction_name'
  
  netVisual_bubble(uterus_cellchat,
                   sources.use = c('POSTN_fibroblast','APCDD1_fibroblast'),
                   targets.use = c('Venous EC','Lymphatic EC'),
                   comparison = c(1, 2, 3,4),
                   pairLR.use = pairLR.use)
  
  #Differences in horizontal communication between ligand-receptor pairs: a selective pathway
  uterus_cellchat@netP[["LSCC"]][["pathways"]]
  pathways.show <- c("COLLAGEN","FN1","MIF","TENASCIN") 
  
  netVisual_bubble(uterus_cellchat, 
                   sources.use = c('POSTN_fibroblast','APCDD1_fibroblast'), 
                   targets.use = c('Venous EC','Lymphatic EC'), 
                   comparison = c(1, 2, 3,4), 
                   signaling = pathways.show)
  
  p <- netVisual_bubble(uterus_cellchat,
                        sources.use = c('POSTN_fibroblast','APCDD1_fibroblast'),
                        targets.use = c('Venous EC','Lymphatic EC'),
                        comparison = c(1, 2, 3,4),
                        signaling = pathways.show)+
    theme(legend.position = 'top',
          legend.key.width = unit(1,'cm'),
          legend.margin = margin(0.5,0.5,0,0,'cm'),
          plot.margin=unit(c(2.5, 2.5, 2.5, 2.5),'cm'))+
    coord_cartesian(clip = 'off')
  
  p2 <- annoSegment(object = p,
                    annoPos = 'top',
                    aesGroup = T,
                    aesGroName = 'source',
                    yPosition = 34, 
                    segWidth = 0.4, 
                    pCol=c("#EB746A",'#26A5DF'), 
                    addText=T, 
                    textSize = 8, 
                    textCol = c("black","black"))
  
  
}

##===============================================================================================================##
##=========================Step05：Nichenet analysis=============================================================##
##===============================================================================================================##
if(F){
  library(Seurat)
  library(cowplot)
  library(patchwork)
  library(ggplot2)
  library(nichenetr)
  library(tidyverse)
  library(dplyr)
  library(circlize)
  library(RColorBrewer)
  
  setwd("~/project/VCL_LSCC/Step5_CellChat/Fib_EpC_detail/NicheNet/")
  seurat_obj<-readRDS("Fib_EpC_detail_Cellchat.rds")
  seurat_obj <- subset(seurat_obj, celltype %in% c("Epithelial cell","POSTN_fibroblast"))
  seurat_obj@meta.data$celltype_Group = paste(seurat_obj@meta.data$celltype, seurat_obj@meta.data$Group,sep = "_")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = F)
  ElbowPlot(seurat_obj, ndims = 50)
  pcs <- 1:15
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = pcs)
  seurat_obj <- FindClusters(seurat_obj, resolution = 1)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = pcs)
  seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = pcs)
  
  seurat_obj$celltype_Group <- factor(seurat_obj$celltype_Group,levels = c("Epithelial cell_LSCCP","Epithelial cell_VCP","Epithelial cell_VCL","Epithelial cell_LSCC",
                                                                           "POSTN_fibroblast_LSCCP","POSTN_fibroblast_VCL","POSTN_fibroblast_VCP","POSTN_fibroblast_LSCC"))
  
  celltype_id = "celltype_Group"
  seurat_obj <- SetIdent(seurat_obj, value = seurat_obj@meta.data[[celltype_id]])
  
  ####Read in the NicheNet ligand-receptor network and ligand-target matrix
  # ##the following of human origin
  ligand_target_matrix <- readRDS("~/database/Resource/NicheNet/ligand_target_matrix.rds")
  lr_network <- readRDS("~/database/Resource/NicheNet/lr_network.rds")
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)
  
  organism = "human" # user adaptation required on own dataset
  
  
  
  ####nichnet
  #####1. Define the niches/microenvironments of interest
  #! Important: your receiver cell type should consist of 1 cluster!
  
  niches = list(
    "LSCCP" = list(
      "sender" = "POSTN_fibroblast_LSCCP", 
      "receiver" = "Epithelial cell_LSCCP"), 
    "VCP" = list(
      "sender" = "POSTN_fibroblast_VCP",
      "receiver" = "Epithelial cell_VCP"),
    "VCL" = list(
      "sender" = "POSTN_fibroblast_VCL",
      "receiver" = "Epithelial cell_VCL"),
    "LSCC" = list(
      "sender" = "POSTN_fibroblast_LSCC",
      "receiver" = "Epithelial cell_LSCC")
  ) 
  
  
  #####2. Calculate differential expression between the niches
  #calculate_niche_de Calculate differential expression of cell types in one niche versus all other niches of interest. 
  assay_oi = "RNA" # other possibilities: RNA,...
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
  
  DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  expression_pct = 0.10
  DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
  
  ###Combine sender-receiver DE based on L-R pairs:
  specificity_score_LR_pairs = "min_lfc" 
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
  
  
  #####3. Optional: Calculate differential expression between the different spatial regions
  include_spatial_info_sender = FALSE
  include_spatial_info_receiver = FALSE
  #spatial_info = tibble(celltype_region_oi = "CAF_High", celltype_other_region = "myofibroblast_High", niche =  "pEMT_High_niche", celltype_type = "sender") # user adaptation required on own dataset
  #specificity_score_spatial = "lfc"
  # this is how this should be defined if you don't have spatial info
  # mock spatial info
  if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 
  
  if(include_spatial_info_sender == TRUE){
    sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "sender"))
    sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, type = "sender", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
    sender_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
    
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
    
  } else {
    # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
    sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
    
  }
  ## [1] "Calculate Spatial DE between: CAF_High and myofibroblast_High"
  
  if(include_spatial_info_receiver == TRUE){
    receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), spatial_info = spatial_info %>% filter(celltype_type == "receiver"))
    receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, type = "receiver", lr_network = lr_network, expression_pct = expression_pct, specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
    receiver_spatial_DE_others = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
    
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
    
  } else {
    # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
    receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }
  
  #####4. ligand activities were calculated and active ligand-target links were inferred
  ##We recommend having between 20 and 1000 genes in the geneset of interest
  lfc_cutoff = 0.25 # recommended for 10x as min_lfc cutoff. 
  specificity_score_targets = "min_lfc"
  
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
  View(DE_receiver_processed_targets)
  
  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  geneset_niche3 = DE_receiver_processed_targets %>% filter(receiver == niches[[3]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  geneset_niche4 = DE_receiver_processed_targets %>% filter(receiver == niches[[4]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche3 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche4 %>% setdiff(rownames(ligand_target_matrix))
  
  print(length(geneset_niche1))
  print(length(geneset_niche2))
  print(length(geneset_niche3))
  print(length(geneset_niche4))
  
  
  top_n_target = 250
  
  
  niche_geneset_list = list(
    "LSCCP" = list(
      "receiver" = niches[["LSCCP"]]$receiver,
      "geneset" = geneset_niche1,
      "background" = background),
    "VCP" = list(
      "receiver" = niches[["VCP"]]$receiver,
      "geneset" = geneset_niche2,
      "background" = background),
    "VCL" = list(
      "receiver" = niches[["VCL"]]$receiver,
      "geneset" = geneset_niche3,
      "background" = background),
    "LSCC" = list(
      "receiver" = niches[["LSCC"]]$receiver,
      "geneset" = geneset_niche4,
      "background" = background)
  )
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
  
  #####5. Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% dplyr::rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% dplyr::rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
  exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% dplyr::rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% dplyr::rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
  dotplot
  
  exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  
  exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  #####6. Expression fraction and receptor
  exprs_sender_receiver = lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 
  View(ligand_scaled_receptor_expression_fraction_df)
  #####7. Prioritization of ligand-receptor and ligand-target links
  prioritizing_weights = c("scaled_ligand_score" = 5,
                           "scaled_ligand_expression_scaled" = 1,
                           "ligand_fraction" = 1,
                           "scaled_ligand_score_spatial" = 0, 
                           "scaled_receptor_score" = 0.5,
                           "scaled_receptor_expression_scaled" = 0.5,
                           "receptor_fraction" = 1, 
                           "ligand_scaled_receptor_expression_fraction" = 1,
                           "scaled_receptor_score_spatial" = 0,
                           "scaled_activity" = 0,
                           "scaled_activity_normalized" = 1,
                           "bona_fide" = 1)
  output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
  prioritization_tables$prioritization_tbl_ligand_target %>% filter(receiver == niches[[1]]$receiver) %>% head(10)
  
  ####8. Visualization of the Differential NicheNet output
  top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% dplyr::rename(top_niche = niche)
  top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% dplyr::rename(top_niche = niche)
  ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(20, prioritization_score) %>% ungroup() # get the top50 ligands per niche
  
  table(ligand_prioritized_tbl_oi$receiver)
  
  outpath <- "./"
  
  receiver_oi = "Epithelial cell_LSCC" 
  sender_oi="POSTN_fibroblast_LSCC"
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% filter(sender == sender_oi) %>% pull(ligand) %>% unique()
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
  lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, prioritized_tbl_oi, prioritization_tables$prioritization_tbl_ligand_receptor, plot_legend = FALSE, heights = NULL, widths = NULL)
  
  pdf(paste(outpath,"/1.Top20ligand_plot.pdf",sep = ""),width = 10,height = 10)
  lfc_plot
  dev.off()
  
  
  ##Ligand expression, activity and target genes
  exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = T, heights = NULL, widths = NULL)
  
  pdf(paste(outpath,"/2.LigandExpressionActivityTargetgenes_plotLegend.pdf",sep = ""),width = 35,height = 15)
  exprs_activity_target_plot$combined_plot
  dev.off()
  
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% filter(sender == sender_oi)%>% top_n(20, prioritization_score) %>% pull(ligand) %>% unique()
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
  exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, prioritized_tbl_oi,  prioritization_tables$prioritization_tbl_ligand_receptor,  prioritization_tables$prioritization_tbl_ligand_target, output$exprs_tbl_ligand,  output$exprs_tbl_target, lfc_cutoff, ligand_target_matrix, plot_legend = FALSE, heights = NULL, widths = NULL)
  
  pdf(paste(outpath,"/3.LigandExpressionActivityTargetgenes_Filterplot.pdf",sep = ""),width = 20,height = 10)
  exprs_activity_target_plot$combined_plot
  dev.off()
  
  ####### L-R pairs
  # barplot(c(1:5),col = c("lavender","#F7F7F7","#CCCCCC", "#969696","#636363"))
  receiver_oi <- c("Epithelial cell_LSCCP","Epithelial cell_VCP","Epithelial cell_VCL","Epithelial cell_LSCC")
  sender_oi <-  c("POSTN_fibroblast_LSCCP","POSTN_fibroblast_VCL","POSTN_fibroblast_VCP","POSTN_fibroblast_LSCC")
  
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(25, prioritization_score) %>% pull(ligand) %>% unique()
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
  
  prioritized_tbl_oi$sender <- factor(prioritized_tbl_oi$sender,levels = sender_oi)
  prioritized_tbl_oi$receiver <- factor(prioritized_tbl_oi$receiver,levels = receiver_oi )
  
  
  colors_sender = c('#9B5B33','#7623f1','#3C77AF','#e184ca') %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
  colors_receiver = c('#93CBAE','#F3AE63','#73558B','#AC1B1F')  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())
  circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
  
  pdf(paste(outpath,"/4.Circosplot.pdf",sep = ""),width = 7,height =6)
  circos_output$p_circos
  dev.off()
  pdf(paste(outpath,"/4.Circosplotlegend.pdf",sep = ""),width = 5,height = 5)
  circos_output$p_legend
  dev.off()
  
}
