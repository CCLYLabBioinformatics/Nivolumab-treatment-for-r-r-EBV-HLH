XY_DimPlot <- function (object, dims = c(1, 2), cells = NULL, cols = NULL,
    pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
    shape.by = NULL, order = NULL, label = FALSE, label.size = 4,
    repel = FALSE, cells.highlight = NULL, cols.highlight = "#DE2D26",
    sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE)
{
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    reduction <- reduction %||% DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- Embeddings(object = object[[reduction]])[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    object[["ident"]] <- Idents(object = object)
    orig.groups <- group.by
    group.by <- group.by %||% "ident"
    data[, group.by] <- object[[group.by]][cells, , drop = TRUE]
    data <- data[order(data[,group.by]),]
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- SingleDimPlot(data = data[, c(dims, x, split.by,
            shape.by)], dims = dims, col.by = x, cols = cols,
            pt.size = pt.size, shape.by = shape.by, order = order,
            label = FALSE, cells.highlight = cells.highlight,
            cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
            na.value = na.value)
        if (label) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel,
                size = label.size, split.by = split.by)
        }
        if (!is.null(x = split.by)) {
            plot <- plot + FacetTheme() + facet_wrap(facets = vars(!!sym(x = split.by)),
                ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                  length(x = unique(x = data[, split.by]))
                }
                else {
                  ncol
                })
        }
        return(plot)
    })
    if (!is.null(x = split.by)) {
        ncol <- 1
    }
    if (combine) {
        plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
    }
    return(plots)
}
environment(XY_DimPlot) <- asNamespace('Seurat')

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
          default_order <- c("umap", "tsne", "pca","dm")
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
      xlims <- c(min(data[, dims[1]]), max(data[,dims[1]]))
      ylims <- c(min(data[, dims[2]]), max(data[,dims[2]]))
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
