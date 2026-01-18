Single-cell transcriptomic analysis of the progression from normal laryngeal mucosa to precancerous lesions (vocal cord leukoplakia) to laryngeal squamous cell carcinoma (LSCC), including the assessment of changes in cell type proportions, 
functional characteristics of cell types, tumor cell module features, as well as intercellular interactions and communication during LSCC development. 
The analytical framework integrates several key methodologies, including quality control, batch effect correction, unsupervised clustering, copy number variation (CNV) inference, hdWGCNA, pathway enrichment analysis, 
and ligand-receptor interaction mapping.   Together, these approaches provide a comprehensive view of the LSCC tumor microenvironment.

The analysis is conducted using standard scRNA-seq workflows implemented in the Seurat package (v5.1.0) for R (v4.4.1) and other supporting R packages.
Requirements

R v4.4.1

Seurat v5.1.0

DoubletFinder v2.0.4

Harmony v1.2.1

InferCNV v1.20.0

hdWGCNA v0.4.08

fgsea v1.30.0

GSVA v1.52.3

Monocle v2.30.1

CellChat v2.1.2

NicheNet v2.2.0

 #Attention, please! The scPRIT package relies on specific R and Python environments and database resources, 
  which the developers integrate into the docker image analysis environment. 
  Specific access method please see the following link: https://mp.weixin.qq.com/s/sR12NVP179xZK5VEgRSxmA
