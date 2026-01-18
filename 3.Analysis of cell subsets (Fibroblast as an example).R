##====================================================================================##
##=========================Step1：Annotation of fibroblast subsets====================##
##====================================================================================##
if(F){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(parallel)
  library(scPRIT)
  library(trqwe)
  library(ggplot2)
  library(ggsci)
  library(tidydr)
  setwd("~/project/VCL_LSCC/Step3.1_Fibro/")
  
  scRNA <- readRDS("~/project/VCL_LSCC/Step1_QC&Define/sco.hny_anno.rds")
  table(scRNA$CellType, scRNA$Group)
  scRNA <- subset(scRNA, CellType %in% c("Fibroblast"))
  colnames(scRNA@meta.data)
  counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
  meta <- scRNA@meta.data[,c(4:14)]
  scRNA <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)
  
  # harmony
  scRNA <- NormalizeData(scRNA)
  scRNA <- FindVariableFeatures(scRNA)
  scRNA <- ScaleData(scRNA)
  scRNA <- RunPCA(scRNA, verbose = F)
  ElbowPlot(scRNA, ndims = 50)
  scRNA <- harmony::RunHarmony(scRNA, group.by.vars = "orig.ident")
  pcs <- 1:20
  scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = pcs)
  scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = pcs) %>% FindClusters()
  saveRDS(scRNA, file = "sco.Fibro_hny.rds")
  
  DimPlot(scRNA, group.by = c("seurat_clusters","orig.ident","Group"), label = T)
  
  DotPlot(scRNA, features = c("COL1A1","COL1A2","COL3A1","COL14A1","DCN"))
  DotPlot(scRNA, features = c("APOD","CXCL12","CRABP1","WNT5A"))
  
  VlnPlot(scRNA, features = c("CCL17","CCL19","CCL22","CD274","PDCD1LG2"),pt.size= 0)
  
  ### cluster marker
  if(F){
    dir.create("ClusterDEG")
    # cluster DEG
    clusterDEG <- FindAllMarkers(scRNA, only.pos = T, logfc.threshold = 0.5, slot = "data")
    write.csv(clusterDEG, file = "ClusterDEG/clusterDEG.csv", row.names = F)
    saveRDS(clusterDEG, file = "ClusterDEG/clusterDEG.rds")
    # cluster top markers
    topDEGs <- filterAllMarkers(clusterDEG, topN = 30, log2FC.T = 0.5, 
                                p_val_adj.T = 0.05, rm.MT = T, rm.RB = T)
    write.csv(topDEGs, file = "ClusterDEG/clusterDEG_top30.csv", row.names = F)
    saveRDS(topDEGs, file = "ClusterDEG/clusterDEG_top50.rds")
    head(topDEGs)
    }
  
  scRNA$celltype <- recode(scRNA$seurat_clusters,
                           "0" = "CFD_fibroblast",
                           "1" = "CFD_fibroblast",
                           "2" = "CFD_fibroblast",
                           "3" = "POSTN_fibroblast",
                           "4" = "ACTA2_myofibroblast",
                           "5" = "CFD_fibroblast",
                           "6" = "CFD_fibroblast",
                           "7" = "APCDD1_fibroblast",
                           "8" = "ACTA2_myofibroblast",
                           "9" = "CFD_fibroblast",
                           "10" = "POSTN_fibroblast",
                           "11" = "CFD_fibroblast",
                           "12" = "POSTN_fibroblast",
                           "13" = "CFD_fibroblast",
                           "14" = "CD74_fibroblast",
                           "15" = "Proliferating fibroblast",
                           "16" = "CFD_fibroblast")
  scRNA$celltype <- factor(scRNA$celltype,levels = c('CFD_fibroblast','POSTN_fibroblast','APCDD1_fibroblast',
                                                     'CD74_fibroblast','ACTA2_myofibroblast','Proliferating fibroblast'))
  saveRDS(scRNA, file = "sco.Fibro_hny_anno.rds")
  
  colors1<-c('#df772b','#de554c','#306311','#e184ca','#b2dfe8','#5ec9b8','#e8e684')
  colors1<-c('#EE934E','#D1352B','#9B5B33','#7623f1','#3C77AF','#e184ca','#D2EBC8','#FCED82',
             '#7DBFA7','#AECDE1','#F5CFE4','#8FA4AE','#F5D2A8','#BBDD78','#314e9b')
  

  DimPlot(scRNA, reduction = "umap", group.by = "celltype",label = F,cols = colors1,pt.size = 1) + 
    theme_dr(xlength = 0.2, 
             ylength = 0.2,
             arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(face = 2,hjust = 0.03))

  
  DimPlot(scRNA, group.by = "celltype", reduction = "umap",label = F,raster = F,cols = colors1)+
    labs(x = "UMAP1", y = "UMAP2",title = "CellType") +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  ##Changes in the proportion of cell subsets
  scRNA$Group <- factor(scRNA$Group,levels = c('LSCCP','VCP','VCL','LSCC'))
  subF <- scRNA@meta.data
  SufibrosDisA <- table(subF[,c("Group","celltype")]) %>% 
    data.frame %>% set_colnames(c("Group","CellTypes","Number"))
  
  SufibrosDisA_Tissue <- lapply(split(SufibrosDisA,SufibrosDisA$Group),function(X){
    X%>% dplyr::mutate(Per=100*Number/sum(Number))
  }) %>% dplyr::bind_rows(.)
  
  SufibrosDisA_Tissue$Group <- factor(SufibrosDisA_Tissue$Group,levels = c('LSCCP','VCP','VCL','LSCC'))
  SufibrosDisA_Tissue$CellTypes <- factor(SufibrosDisA_Tissue$CellTypes,levels = names(table(scRNA$celltype)))
  SufibrosDisA_TissueP <- arrange(SufibrosDisA_Tissue,CellTypes) %>% dplyr::filter(Number != 0 )
  install.packages("openxlsx")
  library(openxlsx)
  openxlsx::write.xlsx(SufibrosDisA_TissueP,file="Fibro_Proportion.xlsx")
  
  library(ggplot2)
  library(dplyr)
  
  ggplot(SufibrosDisA_TissueP, aes(x = Group, y = Per, group = CellTypes, color = CellTypes)) +
    geom_line(linewidth = 1.2) +  
    geom_point(size = 3) +         
    labs(
      title = "Percentage of Cell Types Across Groups",
      x = "Group",
      y = "Percentage (%)",
      color = "Cell Types"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "right",
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) 
  

  selected_celltypes <- c("POSTN_fibroblast")
  df_filtered <- SufibrosDisA_TissueP %>% 
    filter(CellTypes %in% selected_celltypes)
  celltype_colors <- c(
    "CFD_fibroblast" = '#df772b',
    "POSTN_fibroblast" = '#D1352B')

  ggplot(df_filtered, aes(x = Group, y = Per, group = CellTypes, color = CellTypes)) +
    geom_line(size = 1.2, alpha = 0.8, linetype = "dashed") +
    geom_point(size = 6, shape = 19) +
    geom_text(aes(label = sprintf("%.1f%%", Per)), 
              vjust = -1.2, size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = celltype_colors) +
    labs(
      title = "Percentage of CFD_fibroblast Across Sample Types",
      x = "Group",
      y = "Percentage (%)",
      color = "Cell Types"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(face = "bold", size = 13),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.border = element_rect(color = "grey30", fill = NA, size = 0.8)
    ) +
    scale_y_continuous(limits = c(0, 30), expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(labels = c("LSCCP" = "LSCCP\n(n=16,200)", 
                                "VCP" = "VCP\n(n=13,875)", 
                                "VCL" = "VCL\n(n=29,900)", 
                                "LSCC" = "LSCC\n(n=5,900)"))
  
  
  ##########Dot plot of marker genes
  scRNA <- readRDS("sco.Fibro_hny_anno.rds")
  Idents(scRNA) <- "celltype"
  celltypeDEG <- FindAllMarkers(scRNA, only.pos = T,logfc.threshold = 0.5, slot = "data")
  topDEGs <- filterAllMarkers(celltypeDEG, topN = 10, log2FC.T = 0.5, 
                              p_val_adj.T = 0.05, rm.MT = T, rm.RB = T)
  
  gene <- c('CFD','MFAP5',  # CFD_fibroblast
            'POSTN', 'TNC', # POSTN_fibroblast
            'APCDD1', 'CACNA2D3', #APCDD1_fibroblast
            'CD74','HLA-DRA',# CD74_fibroblast
            'ACTA2','RGS5.1',# ACTA2_myofibroblast
            'MKI67','TOP2A' # Proliferating fibroblast
  ) 
  
  pdf(file="DotPlot.pdf", w=8, h=4, family = "sans")
  DotPlot(scRNA, 
          features = gene,
          group.by = "celltype") + 
    scale_color_gradientn(colours = c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                      "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A",
                                      "#FD705F", "#E12234","#840729")) +
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
  
  ###GSVA分析
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(ggpubr)
  setwd("~/project/VCL_LSCC/Step3.1_Fibro/")
  dir.create("GSVA")
  sco <- readRDS("sco.Fibro_hny_anno.rds")
  
  genesets <- read.gmt("~/database/msigdb/h.all.v2024.1.Hs.symbols.gmt")
  genesets <- split(genesets$gene, genesets$term)
  
  avg.es <- runGSVA(sco, genesets, method = "aucell", useAverage = T, group.by = "celltype")
  rownames(avg.es) <- gsub("HALLMARK_", "", rownames(avg.es))
  write.csv(avg.es, "GSVA/gsva_hallmark.csv", row.names = T)

  group1 <- data.frame(Sample = factor(colnames(avg.es)), row.names = colnames(avg.es))
  ann_cols = c('#EE934E','#D1352B','#9B5B33','#7623f1','#3C77AF')
  
  gsvaHeatmap(avg.es, group = group1, ann_cols = ann_cols, fontsize = 8,
              show_colnames = T, cluster_rows = T, cluster_cols = T, scale = "row", 
              width = 8, height = 8.5, filename = "GSVA/CellType_GSVA_Heatmap.pdf")
  

  scRNA <- readRDS("sco.Fibro_hny_anno.rds")
  sco <- subset(scRNA, celltype %in% c("POSTN_fibroblast",'CFD_fibroblast'))
  genesets <- read.gmt("~/database/msigdb/h.all.v2024.1.Hs.symbols.gmt")
  genesets <- split(genesets$gene, genesets$term)
  sc.es <- runGSVA(sco, genesets, method = "aucell", ncores = 8)
  rownames(sc.es) <- gsub("HALLMARK_", "", rownames(sc.es))
  
  sc.diff <- gsvaDiff(gsva.result = sc.es, ident.1 = "POSTN_fibroblast", ident.2 = 'CFD_fibroblast',
                      sc.level = T, sc.group = sco$celltype)
  
  write.csv(sc.diff, "GSVA/Fibro_gsva_diff_hallmark.csv",, row.names = F)

  p <- gsvaDivBarplot(gsva.diff = "GSVA/Fibro_gsva_diffsig_hallmark.csv", value.key = "t", ggcharts = F,cols = c('#D1352B','#EE934E'),)
  ggsave("GSVA/gsva_diff_sig.pdf", p, width = 10, height = 5)
  
  ###Proportion statistics of cell types
  library(ggpubr)
  library(RColorBrewer)
  library(tidyr)
  library(Seurat)
  library(patchwork)
  scRNA <- readRDS("sco.Fibro_hny_anno.rds")
  levels(factor(scRNA$celltype))
  pal <- c('#EE934E','#D1352B','#9B5B33','#7623f1','#3C77AF','#e184ca')
  
  bar.df=scRNA@meta.data
  text.df1=as.data.frame(table(bar.df$celltype))
  p1=bar.df%>%ggplot(aes(x=celltype))+geom_bar(aes(fill=Group),position = "fill")+
    scale_x_discrete("")+
    scale_y_continuous("cell ratio",expand = c(0.02,0),labels = scales::label_percent())+
    scale_fill_manual(values = c('#93CBAE','#F3AE63','#73558B','#AC1B1F'))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      legend.title = element_blank())
  ggsave("cell_group.bar1.pdf",width = 12,height = 18,units = "cm")
  
  p2=bar.df%>%ggplot(aes(x=Group))+geom_bar(aes(fill=celltype),position = "fill")+
    scale_x_discrete("")+
    scale_y_continuous("cell ratio",expand = c(0.02,0),labels = scales::label_percent())+
    scale_fill_manual(values = pal)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=45,vjust=1,hjust=1),
      legend.title = element_blank())
  ggsave("cell_group.bar2.pdf",width = 12,height = 18,units = "cm")

  p1+p2+plot_layout(ncol = 2)
  ggsave("PatientCelltype_bar.pdf",height = 9,width = 12)
  
  ###GSEA analysis
  library(Seurat)
  library(fgsea)
  library(msigdbr)
  library(tidyverse)
  library(pheatmap)
  library(scPRIT)
  setwd("~/project/VCL_LSCC/Step3.1_Fibro/GSEA/")
  scRNA <- readRDS("sco.Fibro_hny_anno.rds")
  
  # Download the Hallmark gene set
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol) %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
  # Download the KEGG gene set
  mdb <- msigdbr(species = "Homo sapiens", category = "C2")
  mdb_kegg = mdb [grep("^KEGG",mdb $gs_name),]
  fgsea_sets<- mdb_kegg %>% split(x = .$gene_symbol, f = .$gs_name)
  
  # Download the GOBP gene set
  mdb <- msigdbr(species = "Homo sapiens", category = "C5")
  mdb_GOBP = mdb [grep("^GOBP",mdb $gs_name),]
  fgsea_sets<- mdb_GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
  
  Idents(scRNA) <- "celltype"
  
  markers <- FindAllMarkers(scRNA,test.use = "wilcox",only.pos = FALSE,
                            logfc.threshold = 0.1,min.pct = 0.1,verbose = FALSE)
  
  gsea_results <- list()
  
  for (celltype in unique(markers$cluster)) {
    celltype_markers <- markers %>% 
      filter(cluster == celltype) %>%
      arrange(desc(avg_log2FC)) 
    
    ranked_genes <- setNames(celltype_markers$avg_log2FC, celltype_markers$gene)
    
    # The GSEA analysis was run
    set.seed(123) 
    gsea_res <- fgsea(
      pathways = fgsea_sets,
      stats = ranked_genes,
      minSize = 15,     
      maxSize = 500,    
      eps = 0.0,        
      nPermSimple = 10000 
    )
    
    gsea_res$celltype <- celltype
    gsea_results[[celltype]] <- gsea_res
  }
  
  combined_gsea <- bind_rows(gsea_results)
  
  # Creating the NES Matrix
  nes_matrix <- combined_gsea %>%
    dplyr::select(pathway, celltype, NES) %>%
    pivot_wider(
      names_from = celltype,
      values_from = NES,
      values_fill = 0 
    ) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  

  pheatmap(
    mat = nes_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",          
    cluster_rows = TRUE,    
    cluster_cols = TRUE,     
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = 'white',       
    main = "Hallmark Pathway NES Scores by Cell Type",
    angle_col = 45  
  )
  
  ###Draw the picture after selecting GSEA
  library(openxlsx)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)

  gsea_data <- read.xlsx("gsea_CAF.xlsx", sheet = 1)

  colnames(gsea_data) <- c("Pathway_Category", "Pathway", "CFD_fibroblast", "POSTN_fibroblast", 
                           "APCDD1_fibroblast", "CD74_fibroblast", "ACTA2_myofibroblast", 
                           "Proliferating_fibroblast")

  gsea_data <- gsea_data[!is.na(gsea_data$Pathway), ]
  
  nes_cols <- c("CFD_fibroblast", "POSTN_fibroblast", "APCDD1_fibroblast", 
                "CD74_fibroblast", "ACTA2_myofibroblast", "Proliferating_fibroblast")
  
  gsea_data[nes_cols] <- lapply(gsea_data[nes_cols], function(x) {
    as.numeric(as.character(x))
  })
  
  row_annotation <- data.frame(
    Pathway_Category = gsea_data$Pathway_Category,
    row.names = gsea_data$Pathway
  )
  
  nes_matrix <- gsea_data %>%
    select(-Pathway_Category) %>%
    column_to_rownames("Pathway") %>%
    as.matrix()
  
  category_colors <- list(
    Pathway_Category = c(
      GOBP = "#ff7f0e",
      KEGG = "#1f77b4"
    )
  )
  
  gobp_count <- sum(gsea_data$Pathway_Category == "GOBP")
  kegg_count <- sum(gsea_data$Pathway_Category == "KEGG")
  
  pheatmap(
    mat = nes_matrix,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_row = row_annotation,
    annotation_colors = category_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 12,
    fontsize_col = 12,
    fontsize = 15,
    angle_col = 45,
    border_color = 'white',
    cellwidth = 16,cellheight = 16,  
    gaps_row = gobp_count,  
    na_col = "grey90")

}

