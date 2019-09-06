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

XY_FeaturePlot <- function (object, features, dims = c(1, 2), cells = NULL, cols = c("lightgrey",
      "blue"), pt.size = NULL, min.cutoff = NA, max.cutoff = NA,
      reduction = NULL, split.by = NULL, shape.by = NULL, blend = FALSE,
      blend.threshold = 0.5, order = NULL, label = FALSE, label.size = 4,
      ncol = NULL, combine = TRUE, coord.fixed = FALSE, ...)
 {
      no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
              size = 14, margin = margin(r = 7)))
      if (is.null(reduction)) {
          default_order <- c("umap", "tsne", "pca")
          reducs <- which(default_order %in% names(object@reductions))
          reduction <- default_order[reducs[1]]
      }
      if (length(x = dims) != 2 || !is.numeric(x = dims)) {
          stop("'dims' must be a two-length integer vector")
      }
      if (blend && length(x = features) != 2) {
          stop("Blending feature plots only works with two features")
      }
      dims <- paste0(Key(object = object[[reduction]]), dims)
      cells <- cells %||% colnames(x = object)
      data <- FetchData(object = object, vars = c(dims, features),
          cells = cells)
      features <- colnames(x = data)[3:ncol(x = data)]
      min.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = min(data[,
              feature]), no = cutoff))
      }, cutoff = min.cutoff, feature = features)
      max.cutoff <- mapply(FUN = function(cutoff, feature) {
          return(ifelse(test = is.na(x = cutoff), yes = max(data[,
              feature]), no = cutoff))
      }, cutoff = max.cutoff, feature = features)
      check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
          max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
      if (length(x = check.lengths) != 1) {
          stop("There must be the same number of minimum and maximum cuttoffs as there are features")
      }
      brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
          ]$maxcolors, no = length(x = cols))
      data[, 3:ncol(x = data)] <- sapply(X = 3:ncol(x = data),
          FUN = function(index) {
              data.feature <- as.vector(x = data[, index])
              min.use <- SetQuantile(cutoff = min.cutoff[index -
                  2], data.feature)
              max.use <- SetQuantile(cutoff = max.cutoff[index -
                  2], data.feature)
              data.feature[data.feature < min.use] <- min.use
              data.feature[data.feature > max.use] <- max.use
              if (brewer.gran == 2) {
                  return(data.feature)
              }
              data.cut <- if (all(data.feature == 0)) {
                  0
              }
              else {
                  as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                    breaks = brewer.gran)))
              }
              return(data.cut)
          })
      colnames(x = data)[3:ncol(x = data)] <- features
      rownames(x = data) <- cells
      data$split <- if (is.null(x = split.by)) {
          RandomName()
      }
      else {
          switch(EXPR = split.by, ident = Idents(object = object)[cells],
              object[[split.by, drop = TRUE]][cells])
      }
      if (!is.factor(x = data$split)) {
          data$split <- factor(x = data$split)
      }
      plots <- vector(mode = "list", length = ifelse(test = blend,
          yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
      xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[,
          dims[1]])))
      ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[,
          dims[2]])))
      if (blend) {
          ncol <- 4
          color.matrix <- BlendMatrix(col.threshold = blend.threshold)
          colors <- list(color.matrix[, 1], color.matrix[1, ],
              as.vector(x = color.matrix))
      }
      for (i in 1:length(x = levels(x = data$split))) {
          ident <- levels(x = data$split)[i]
          data.plot <- data[as.character(x = data$split) == ident,
              , drop = FALSE]
          if (blend) {
              data.plot <- cbind(data.plot[, dims], BlendExpression(data = data.plot[,
                  features[1:2]]))
              features <- colnames(x = data.plot)[3:ncol(x = data.plot)]
          }
          for (j in 1:length(x = features)) {
              feature <- features[j]
              if (blend) {
                  cols.use <- as.numeric(x = as.character(x = data.plot[,
                    feature])) + 1
                  cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
              }
              else {
                  cols.use <- NULL
              }
              data.plot <- data.plot[order(data.plot[,j+2],decreasing=F),]
              plot <- SingleDimPlot(data = data.plot[, c(dims,
                  feature)], dims = dims, col.by = feature, pt.size = pt.size,
                  cols = cols.use, label.size = label.size) +
                  scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
                  theme_cowplot()
              if (length(x = levels(x = data$split)) > 1) {
                  plot <- plot + theme(panel.border = element_rect(fill = NA,
                    colour = "black"))
                  plot <- plot + if (i == 1) {
                    labs(title = feature)
                  }
                  else {
                    labs(title = NULL)
                  }
                  if (j == length(x = features) && !blend) {
                    suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident)) +
                      no.right)
                  }
                  if (j != 1) {
                    plot <- plot + theme(axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                      axis.title.y.left = element_blank())
                  }
                  if (i != length(x = levels(x = data$split))) {
                    plot <- plot + theme(axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                      axis.title.x = element_blank())
                  }
              }
              else {
                  plot <- plot + labs(title = feature)
              }
              if (!blend) {
                  plot <- plot + guides(color = NULL)
                  if (length(x = cols) == 1) {
                    plot <- plot + scale_color_brewer(palette = cols)
                  }
                  else if (length(x = cols) > 1) {
                    plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols,
                      guide = "colorbar"))
                  }
              }
              if (coord.fixed) {
                  plot <- plot + coord_fixed()
              }
              plot <- plot
              plots[[(length(x = features) * (i - 1)) + j]] <- plot
          }
      }
      if (blend) {
          blend.legend <- BlendMap(color.matrix = color.matrix)
          for (i in 1:length(x = levels(x = data$split))) {
              suppressMessages(expr = plots <- append(x = plots,
                  values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                    1, yes = levels(x = data$split)[i], no = "")),
                    expand = c(0, 0)) + labs(x = features[1], y = features[2],
                    title = if (i == 1) {
                      paste("Color threshold:", blend.threshold)
                    } else {
                      NULL
                    }) + no.right), after = 4 * i - 1))
          }
      }
      plots <- Filter(f = Negate(f = is.null), x = plots)
      if (combine) {
          if (is.null(x = ncol)) {
              ncol <- 2
              if (length(x = features) == 1) {
                  ncol <- 1
              }
              if (length(x = features) > 6) {
                  ncol <- 3
              }
              if (length(x = features) > 9) {
                  ncol <- 4
              }
          }
          ncol <- ifelse(test = is.null(x = split.by) || blend,
              yes = ncol, no = length(x = features))
          legend <- if (blend) {
              "none"
          }
          else {
              split.by %iff% "none"
          }
          plots <- CombinePlots(plots = plots, ncol = ncol, legend = legend,
              nrow = split.by %iff% length(x = levels(x = data$split)))
      }
      return(plots)
}
environment(XY_FeaturePlot) <- asNamespace('Seurat')

HLH_T1 <- readRDS("./anno_HLH_T1.rds")
HLH_T2 <- readRDS("./anno_HLH_T2.rds")
HLH_T3 <- readRDS("./anno_HLH_T3.rds")
HLH_T4 <- readRDS("./anno_HLH_T4.rds")

Idents(HLH_T3,cells=WhichCells(HLH_T3,idents="CD 3+ CD8- CD4- T")) <- "CD4 T like"
Idents(HLH_T4,cells=WhichCells(HLH_T4,idents="CD 3+ CD8- CD4- T")) <- "Naive CD4 T"
Idents(HLH_T4,cells=WhichCells(HLH_T4,idents="FCGR3A+ Monocyte")) <- "CD16+ Monocyte"
Idents(HLH_T1,cells=WhichCells(HLH_T1,idents="CD14+ Monocyte")) <- "CD14+ Monocyte act"
Idents(HLH_T1,cells=WhichCells(HLH_T1,idents="CD4 T helper")) <- "CD4 T"
Idents(HLH_T1,cells=WhichCells(HLH_T1,idents="CD8B+ T")) <- "Naive T"
Idents(HLH_T1,cells=WhichCells(HLH_T1,idents="CD8AB+ exhausted T")) <- "CD8 T act/ex"
Idents(HLH_T2,cells=WhichCells(HLH_T2,idents="CD8B+ T")) <- "Naive T"
Idents(HLH_T2,cells=WhichCells(HLH_T2,idents="MKI67 high CD8AB+ T")) <- "MKI67+ CD8 T act/ex"
Idents(HLH_T2,cells=WhichCells(HLH_T2,idents="CD8AB+ T")) <- "CD8 T act/ex"
Idents(HLH_T3,cells=WhichCells(HLH_T3,idents="CD8AB+ T")) <- "CD8 T"
Idents(HLH_T3,cells=WhichCells(HLH_T3,idents="MKI67 high CD8AB+ T")) <- "MKI67+ CD8 T"
Idents(HLH_T3,cells=WhichCells(HLH_T3,idents="CD8A+ T")) <- "CD8 T"
Idents(HLH_T3,cells=WhichCells(HLH_T3,idents="CD4 T like")) <- "CD4 T"
Idents(HLH_T4,cells=WhichCells(HLH_T4,idents="CD8AB+ T")) <- "CD8 T"

HLH_T4$new_anno <- Idents(HLH_T4)
HLH_T3$new_anno <- Idents(HLH_T3)
HLH_T1$new_anno <- Idents(HLH_T1)
HLH_T2$new_anno <- Idents(HLH_T2)

p1 <- DimPlot(object = HLH_T1, reduction = "tsne", label = FALSE,repel=T,cols=c(colo_schem[14],colo_schem[1],colo_schem[4],
  colo_schem[12],colo_schem[11],colo_schem[5],colo_schem[8])) +labs(title="State 1")
p2 <- DimPlot(object = HLH_T2, reduction = "tsne", label = FALSE,repel=T,cols=c(colo_schem[14],colo_schem[3],colo_schem[1],
  colo_schem[7],colo_schem[6],colo_schem[5],colo_schem[8],colo_schem[10],colo_schem[4],colo_schem[9])) +labs(title="State 2")
p3 <- DimPlot(object = HLH_T3, reduction = "tsne", label = FALSE,repel=T,cols=c(colo_schem[4],colo_schem[2],colo_schem[15],colo_schem[6],
  colo_schem[5],colo_schem[9])) +labs(title="State 3")
p4 <- DimPlot(object = HLH_T4, reduction = "tsne", label = FALSE,repel=T,cols=c(colo_schem[2],colo_schem[4],colo_schem[16],
  colo_schem[11],colo_schem[9])) +labs(title="State 4")
aa <- plot_grid(p1, p2,p3,p4,nrow=1,ncol=4)
ggsave("./v3_HLH_ALL_TSNE.svg", plot=aa,width = 28, height = 5,dpi=1080)


#merge关注的细胞群
HLH_T1_sel <- subset(HLH_T1,cells=WhichCells(HLH_T1,idents="CD8 T act/ex"))
HLH_T2_sel <- subset(HLH_T2,cells=WhichCells(HLH_T2,idents=c("CD8 T act/ex","MKI67+ CD8 T act/ex")))
HLH_T3_sel <- subset(HLH_T3,cells=WhichCells(HLH_T3,idents=c("CD8 T","MKI67+ CD8 T")))
HLH_T4_sel <- subset(HLH_T4,cells=WhichCells(HLH_T4,idents="CD8 T"))

GO_LEUKOCYTE_DEGRANULATION <- c("ANXA3","CCL3","CHGA","CORO1A","CPLX2","GATA2","HCK","KIT","KLRF2","LAT","LAT2","LYN","MILR1","MRGPRX2","NR4A3","PI4K2A","PIK3CD","PIK3CG","PLA2G3","PPBP","RAB27A","RASGRP1","S100A13","SNAP23","STXBP2","STXBP3","UNC13D","VAMP2","VAMP7","VAMP8")
EXOCYTOSIS <- c("AKAP3","ARFGEF1","ARFGEF2","CADPS","CCL3","CCL5","CCL8","CPLX1","CPLX2","LIN7A","NLGN1","PKDREJ","RAB26","RABEPK","RIMS1","RPH3AL","SCIN","SCRN1","SEPT5","SNAP29","SYT1","SYTL4","VAMP3","VTI1B","YKT6")
costimulatory_molecule <- c("TNFRSF4", "ICOS","CD27", "ITGB2","CD2","TNFRSF14", "TNFRSF18", "TNFRSF8", "CD28")
CO_STIMULATORY <- c("CD40L", "TNFRSF9", "HAVCR1","CD28","ICOS","CD27","TNFRSF14","CD40L","CD40LG","TNFRSF9","TNFRSF4","TNFRSF25","TNFRSF18","TNFRSF8","SLAMF1","CD2","CD226","ITGB2")
CO_INHIBITORY <- c("LAG3","CTLA4","CD274","PDCD1","CD80","CD160","BTLA","VSIR","LAIR1","HAVCR2","TIMD4","CD244","TIGIT")
HLH_related_gene <- c("PRF1","UNC13D","STX11","STXBP2","XIAP","SH2D1A","RAB27A","MYO5A","AP3B1","LYST","ITK","SLC7A7","BLOC1S6","CD27","MAGT1","ADA","BTK","IL2RA","IL2RG","MVK","PNP","WAS","FAS","RECQL4","RAG1","RAG2")

HLH_related_gene1 <- c("PRF1")
HLH_related_gene2 <- c("UNC13D","STX11","STXBP2","XIAP","RAB27A","LYST","AP3B1")
HLH_related_gene3 <- c("NLRC4")
HLH_related_gene4 <- c("CD27","ITK","MAGT1","SH2D1A")

HLH_related_gene <- c("CD27","ITK","MAGT1","SH2D1A","PRF1","UNC13D","STX11","STXBP2","XIAP","RAB27A","LYST","AP3B1","NLRC4")

cc <- intersect(GO_LEUKOCYTE_DEGRANULATION,HLH_related_gene)
cc <- intersect(cc,EXOCYTOSIS)
cc <- intersect(cc,costimulatory_molecule)

HLH_T1_speci_raw <- FetchData(object = HLH_T1_sel, vars = HLH_related_gene,slot="data")
HLH_T2_speci_raw <- FetchData(object = HLH_T2_sel, vars = HLH_related_gene,slot="data")
HLH_T3_speci_raw <- FetchData(object = HLH_T3_sel, vars = HLH_related_gene,slot="data")
HLH_T4_speci_raw <- FetchData(object = HLH_T4_sel, vars = HLH_related_gene,slot="data")

cc <- intersect(colnames(HLH_T1_speci_raw),colnames(HLH_T2_speci_raw))
cc <- intersect(cc,colnames(HLH_T3_speci_raw))
cc <- intersect(cc,colnames(HLH_T4_speci_raw))

HLH_T1_speci_raw <- HLH_T1_speci_raw[,cc]
HLH_T2_speci_raw <- HLH_T2_speci_raw[,cc]
HLH_T3_speci_raw <- HLH_T3_speci_raw[,cc]
HLH_T4_speci_raw <- HLH_T4_speci_raw[,cc]

HLH_T1_speci_t <- as.data.frame(t(HLH_T1_speci_raw))
HLH_T2_speci_t <- as.data.frame(t(HLH_T2_speci_raw))
HLH_T3_speci_t <- as.data.frame(t(HLH_T3_speci_raw))
HLH_T4_speci_t <- as.data.frame(t(HLH_T4_speci_raw))

HLH_T1_t <- data.frame(T1=rowSums(HLH_T1_speci_t)/ncol(HLH_T1_speci_t))
HLH_T2_t <- data.frame(T2=rowSums(HLH_T2_speci_t)/ncol(HLH_T2_speci_t))
HLH_T3_t <- data.frame(T3=rowSums(HLH_T3_speci_t)/ncol(HLH_T3_speci_t))
HLH_T4_t <- data.frame(T4=rowSums(HLH_T4_speci_t)/ncol(HLH_T4_speci_t))

all_cbind <- cbind(HLH_T1_t,HLH_T2_t)
all_cbind <- cbind(all_cbind,HLH_T3_t)
all_cbind <- cbind(all_cbind,HLH_T4_t)

require(RColorBrewer)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

test <- all_cbind
test <- as.data.frame(t(apply(test, 1, function(x) (x-mean(x))/sd(x))))
test <- test[order(test$T1,decreasing=F),]
test[test > 1.2] <- 1.2
test[test < -1.2] <- -1.2

library(ggplot2)
library(pheatmap)
ph_HLH_related_gene <- pheatmap(test,cluster_rows = F,col=rev(cols),clustering_method="ward.D2",
         cluster_cols = F,show_rownames=T,show_colnames=T,main="HLH_related_gene in ALL CD8 T")

ggsave("./v3_CD8_T_HLH_related_HLH_related_gene_marker.svg", plot=ph_HLH_related_gene,width = 22, height = 5,dpi=1080)







