##HLH-T0-QC
##HLH-T0-QC
##HLH-T0-QC
##HLH-T0-QC
##HLH-T0-QC
##HLH-T0-QC
****************************************
1-cluster-pre-workflow-new-way
****************************************
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(Matrix)
  library(proxy)
  library(gplots)
  library(Rtsne)
  library(densityClust)
  library(irlba)
  library(monocle)
  library(plyr)
  library(DOSE)
  library(clusterProfiler)
  library(topGO)
  library(pathview)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(cowplot)
  library(ggplot2)
})
XY_Read10X <- function(data.dir = NULL){
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (! dir.exists(run)){
      stop("Directory provided does not exist")
    }
    if(!grepl("\\/$", run)){
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv.gz")
    gene.loc <- paste0(run, "features.tsv.gz")
    matrix.loc <- paste0(run, "matrix.mtx.gz")
    if (!file.exists(barcode.loc)){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc)){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc)){
      stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(
        x = as.character(
          x = sapply(
            X = cell.names,
            FUN = ExtractField, field = 1,
            delim = "-"
          )
        )
      )
    }
    rownames(x = data) <- make.unique(
      names = as.character(
        x = sapply(
          X = gene.names,
          FUN = ExtractField,
          field = 2,
          delim = "\\t"
        )
      )
    )
    if (is.null(x = names(x = data.dir))) {
      if(i < 2){
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}
environment(XY_Read10X) <- asNamespace('Seurat')
Human_UH.data <- XY_Read10X(data.dir = "/mnt/data/user_data/xuelan/project/3_singlecell/6_liupengpeng/matrix/humanT1_2_all/outs/filtered_feature_bc_matrix")
colnames(Human_UH.data) <- paste("PP1",colnames(Human_UH.data),sep="_")

aaa_test <- Human_UH.data["EBV",]
setwd(".")
Human_UH <- CreateSeuratObject(counts = Human_UH.data, min.cells = 3, min.features = 200,project = "Patient_state_1")
mito.features <- grep(pattern = "^MT-", x = rownames(x = Human_UH), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = Human_UH, slot = "counts")[mito.features,])/Matrix::colSums(x = GetAssayData(object = Human_UH, slot = "counts"))
Human_UH[["percent.mito"]] <- percent.mito
VlnPlot(object = Human_UH, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size=.1)
FeatureScatter(object = Human_UH, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = Human_UH, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Human_UH <- subset(x = Human_UH, subset = percent.mito < 0.1 & nFeature_RNA <3000)
Human_UH <- NormalizeData(object = Human_UH)
Human_UH <- FindVariableFeatures(object = Human_UH)
length(x = VariableFeatures(object = Human_UH))

Human_UH <- ScaleData(object = Human_UH, features = rownames(x = Human_UH), vars.to.regress = c("nCount_RNA","percent.mito"))
Human_UH <- RunPCA(object = Human_UH, features = VariableFeatures(object = Human_UH))
Human_UH <- ProjectDim(object = Human_UH)
ElbowPlot(object = Human_UH)
Human_UH <- FindNeighbors(object = Human_UH, dims = 1:15)
Human_UH <- FindClusters(object = Human_UH, resolution = 0.4)
Human_UH <- RunTSNE(object = Human_UH, dims = 1:15)
Human_UH <- RunUMAP(object = Human_UH, dims = 1:15)
p1 <- DimPlot(object = Human_UH, reduction = "tsne",label=T)
p2 <- DimPlot(object = Human_UH, reduction = "umap",label=T)
p3 <- DimPlot(object = Human_UH, reduction = "pca",label=T)
plot_grid(p1,p2,p3,nrow=1)
saveRDS(Human_UH, file = "./Human_UH_tutorial.rds")

Human_UH <- FindClusters(object = Human_UH, resolution = c(0.1,0.2,0.4,0.6,0.8,1,1.2))
p1 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.0.1")
p2 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.0.2")
p3 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.0.4")
p4 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.0.6")
p5 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.0.8")
p6 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.1")
p7 <- DimPlot(object = Human_UH, reduction = "umap",label=T,group.by="RNA_snn_res.1.2")
plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow=3)


Human_UH.markers <- FindAllMarkers(object = Human_UH, only.pos = TRUE, min.pct = 0.25,thresh.use = 0.25)
Human_UH.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_logFC)
Human_UH.markers <- Human_UH.markers[Human_UH.markers$p_val_adj < 0.05,]
write.csv(Human_UH.markers,file="state_1_Human_UH.markers.csv")
top10 <- Human_UH.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = Human_UH, features = top10$gene) + NoLegend()



library(dplyr)
setwd(".")
top200 <- read.csv(file="./state_1_Human_UH.markers.csv")
reference_files_selec <- read.csv(file="./pbmc68k.markers.csv")

reference_files_selec <- arrange(reference_files_selec, cluster, desc(p_val))
cell_types <- unique(reference_files_selec$cluster)
all_score <- c()
for (i in cell_types){
y <- reference_files_selec$cluster
specific_symbols <- reference_files_selec[with(reference_files_selec,y==i),]
specific_symbols$order <- c(1:nrow(specific_symbols))
specific_score_all_cluster <- c()
number.group <- table(top200$cluster)
for (number_group in c(0:(nrow(number.group)-1))){
XY <- top200$cluster
cluster <- top200[with(top200,XY==number_group),]
cluster_specific_symbols <- merge(cluster,specific_symbols, by = "gene")
omit_cluster_specific <- na.omit(cluster_specific_symbols)
cluster_specific_score  <- (200*sum(omit_cluster_specific$order))/(nrow(specific_symbols)*(nrow(specific_symbols)+1))
cluster_specific_score <- as.data.frame(cluster_specific_score)
rownames(cluster_specific_score) <- paste("cluster",number_group,sep="_")
XUEYU <- sprintf("specific_score_cluster_%s <- cluster_specific_score",number_group)
eval(parse(text=XUEYU))
XUEYU <- sprintf("specific_score_all_cluster <- rbind(specific_score_all_cluster,specific_score_cluster_%s)",number_group)
eval(parse(text=XUEYU))}
colnames(specific_score_all_cluster) <- paste(i,"score",sep="_")
specific_score_all_cluster <- t(specific_score_all_cluster)
all_score <- rbind(all_score,specific_score_all_cluster)
}
all_score <- round(all_score,3)
write.csv(all_score,file="state_1_human_state_68k_intersect_score.csv")



Human_UH_T1.markers <- read.csv(file="./state_1_Human_UH.markers.csv")
top_marker <- c()
number.group <- length(unique(Human_UH_T1.markers$cluster))-1
for (i in c(0:number.group)){
  y <- Human_UH_T1.markers$cluster
  marker <- Human_UH_T1.markers[with(Human_UH_T1.markers,y==i),]
  j <- i+1
  top_marker[[j]] <- marker
  names(top_marker)[j] <- paste("clu",i,sep="_")
}
gcSampl <- c()
for (i in c(1:length(top_marker))){
t <- top_marker[[i]]
symbol <- as.character(t$gene)
DD <- symbol
t$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = DD,
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
t <- na.omit(t)
y <- t$p_val_adj
names <- t[with(t,y<0.05),]
entrez <- as.character(names$entrez)
gcSampl[[i]] <- entrez
names(gcSampl)[i] <- paste(names(top_marker)[i],"entrez",sep="_")
print(paste("cluster",i,"is done",sep = " "))
}

ck_OFF_HIGH_GO_BP <- compareCluster(geneCluster = gcSampl, fun = "enrichGO",OrgDb = org.Hs.eg.db,readable=T,ont = "BP")
ck_OFF_HIGH_KEGG <- compareCluster(geneCluster = gcSampl, fun = "enrichKEGG",organism="human")
ck_OFF_HIGH_DO <- compareCluster(geneCluster = gcSampl, fun = "enrichDO")

dotplot(ck_OFF_HIGH_GO_BP,showCategory=10,includeAll=FALSE) + theme(axis.text.x  = element_text(angle=45, vjust=1,size=8,hjust = 1)) + labs(title = "GO BP")
dotplot(ck_OFF_HIGH_KEGG,showCategory=10,includeAll=FALSE) + theme(axis.text.x  = element_text(angle=45, vjust=1,size=8,hjust = 1)) + labs(title = "KEGG")
dotplot(ck_OFF_HIGH_DO,showCategory=10,includeAll=FALSE) + theme(axis.text.x  = element_text(angle=45, vjust=1,size=8,hjust = 1)) + labs(title = "DO")
write.csv(ck_OFF_HIGH_GO_BP,file="./state_1_GO_BP.csv")
write.csv(ck_OFF_HIGH_KEGG,file="./state_1_KEGG.csv")
write.csv(ck_OFF_HIGH_DO,file="./state_1_DO.csv")


























