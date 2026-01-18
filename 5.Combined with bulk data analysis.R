library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
rm(list = ls())
setwd("~/project/VCL_LSCC/Step7_correlation/")

sco.Fibro <- readRDS("~/project/VCL_LSCC/Step3.1_Fibro/sco.Fibro_hny_anno.rds")
colnames(sco.Fibro@meta.data)
counts <- GetAssayData(sco.Fibro, assay = "RNA", layer = "counts")
meta <- sco.Fibro@meta.data[,c(4:10,17)]
sco.Fibro <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)

sco.MM <- readRDS("~/project/VCL_LSCC/Step3.2_Myeloid/sco.scRNAMM_hny_anno.rds")
colnames(sco.MM@meta.data)
counts <- GetAssayData(sco.MM, assay = "RNA", layer = "counts")
meta <- sco.MM@meta.data[,c(4:10,19)]
sco.MM <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)

sco.DC <- readRDS("~/project/VCL_LSCC/Step3.2_Myeloid/sco.scRNADC_hny_anno.rds")
colnames(sco.DC@meta.data)
counts <- GetAssayData(sco.DC, assay = "RNA", layer = "counts")
meta <- sco.DC@meta.data[,c(4:10,19)]
sco.DC <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)

sco.tcell <- readRDS("~/project/VCL_LSCC/Step4.1_Tcell/sco.TNK_hny_anno.rds")
colnames(sco.tcell@meta.data)
counts <- GetAssayData(sco.tcell, assay = "RNA", layer = "counts")
meta <- sco.tcell@meta.data[,c(4:10,15)]
sco.tcell <- CreateSeuratObject(counts, min.cells = 3, meta.data = meta)

sco <- merge(sco.Fibro,list(sco.MM,sco.DC,sco.tcell))
sco <- JoinLayers(sco)
sco <- NormalizeData(sco)
unique(sco$celltype)
saveRDS(sco, "scRNA_Correlation.rds")

###Marker genes for all subpopulations were obtained using FindAllMarkers
scRNA<-readRDS("scRNA_Correlation.rds")
scRNA <- ScaleData(scRNA)
Idents(scRNA) <- "celltype"

all_markers <- FindAllMarkers(object = scRNA,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,           
                              test.use = "wilcox",verbose = FALSE)          

# Extract the top 150 marker genes for each subgroup (ranked by adjusted p-value and logFC)
top_genes_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 150) %>%
  ungroup()

# Create a list of gene sets (one gene set per subgroup)
gene_sets <- lapply(unique(top_genes_per_cluster$cluster), function(cluster_id) {
  top_genes_per_cluster %>%
    filter(cluster == cluster_id) %>%
    pull(gene)
})

# Set the gene set name
names(gene_sets) <- unique(top_genes_per_cluster$cluster)

saveRDS(gene_sets, "150cell_subset_gene_sets.rds")

gene<-readRDS("150cell_subset_gene_sets.rds")

###Analysis of correlation
library(limma)
library(sva)
library(GSVA)
library(ggplot2)
library(ggpubr)
setwd("~/project/VCL_LSCC/Step7_correlation/")

#geo gene expression files were read and the data were processed
rt=read.table("GSE27020/GSE27020.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)
#geo=log2(geo+1)

qx=as.numeric(quantile(geo, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  geo[geo<0]=0
  geo=log2(geo+1)}
geo=normalizeBetweenArrays(geo)
boxplot(geo, main="After standardization", las=2)

gene_sets <- readRDS("150cell_subset_gene_sets.rds")
gene_sets <- lapply(gene_sets, function(x) intersect(x, rownames(geo)))

# Create the ssgseaParam object and run GSVA
param <- ssgseaParam(
  exprData = geo,
  geneSets = gene_sets)

es <- gsva(param, verbose = TRUE)
write.csv(es, "ssgesa.csv", row.names = T)

data <- data.frame(X_expr = as.numeric(es["LAMP3+ DC",]),
                   Y_expr = as.numeric(es["CD8+ Tex",]))
cor.test(data$X_expr,data$Y_expr,method = 'pearson')
ggplot(data,aes(x = X_expr, y = Y_expr)) + 
  #xlim(-20,15) + ylim(-15,10) +
  labs(x = "LAMP3+ DC", y = "CD8+ Tex",title = "TCGA_LSCC") +
  geom_point(colour = "#368ad9",size = 4) +
  geom_smooth(method ='lm',color="red") +
  geom_rug(colour = '#93CBAE')+
  stat_cor(method = "pearson",size = 8) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text=element_text(size=15))

library(ggstatsplot)
ggscatterstats(data = data,
               x = X_expr,
               y = Y_expr,
               xlab = "LAMP3+ DC",
               ylab = "SPP1_macrophage",
               bf.message = F)+ 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16),
        legend.title = element_blank(),
        legend.text=element_text(size=15))

###survival analysis
library(readxl)      
library(tidyverse)   
library(survival)    
library(survminer)  
library(gridExtra)   

# 1.Clinical data were read and preprocessed
clinical_data <- read_excel("GSE27020/clinical.xlsx", sheet = "Sheet1") %>%
  mutate(
    DFS_time = as.numeric(str_extract(DFS_time, "\\d+")), 
    DFS_event = as.numeric(str_extract(DFS_event, "\\d")) 
  ) %>%
  select(Sample_ID, DFS_time, DFS_event, Age, Grade)

# 2.The ssGSEA results were processed
ssgsea_data <- t(es) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample_ID") %>%
  mutate(across(-Sample_ID, as.numeric))

##The length of Sample_ID in TCGA was treated differently
ssgsea_data$Sample_ID <- sapply(ssgsea_data$Sample_ID, function(x) {
  parts <- unlist(strsplit(x, "-"))
  paste(parts[1:3], collapse = "-")
})

# 3.Clinical data and ssGSEA results were combined
combined_data <- clinical_data %>%inner_join(ssgsea_data, by = "Sample_ID")
combined_data$futime <- combined_data$futime / 30

res.cut <- surv_cutpoint(combined_data, 
                         time = "DFS_time", 
                         event = "DFS_event", 
                         minprop = 0.3,     
                         variables = c("POSTN_fibroblast")) 

res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(DFS_time, DFS_event) ~POSTN_fibroblast, data = res.cat)
# Customized survival curves
ggsurvplot(fit, data = res.cat,
           surv.median.line = "hv", 
           legend.title = "POSTN_fibroblast",
           legend.labs = c("High", "Low"), 
           xlab = "Time (months)", ylab = "DFS Probability",
           # Add p-value and tervals
           pval = TRUE, 
           pval.method = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           tables.height = 0.2,
           tables.theme = theme_cleantable(), 
           palette = c('#D1352B','#3C77AF'), 
           ggtheme = theme_bw())

