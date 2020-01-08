### quick rename of any cluster in a list of clusters
renameClusters = function(clusters, From, To) {
  newClusters = plyr::mapvalues(x = clusters,
                                from = From,
                                to = To)
  return(newClusters)
}



######
geneSetSum = function(seurat, geneList, listName){
  geneList = intersect(rownames(seurat@scale.data), geneList)
  seurat@meta.data[, listName] = colSums(seurat@scale.data[geneList, ])
  return(seurat)
}


VlnPlotMarkers = function(seurat, markers, outDir, clusterOrder = NULL, plotColours = NULL ){
  markers = markers[markers$Gene %in% rownames(seurat@scale.data), , drop = F]
  ### Looking at marker genes
  for (cellType in unique(markers$CellType)){
    genes =  as.character(markers$Gene[markers$CellType == cellType])
    seurat@meta.data[, cellType] = colSums(seurat@scale.data[genes, ])
    png(paste0(outDir, "/", cellType, ".png"), height = 800, width = 700)
    print(VlnPlot(object = seurat,
                  features.plot = cellType,
                  point.size.use = 0,
                  remove.legend = T,
                  x.lab.rot = T,
                  size.x.use = 12,
                  size.y.use = 12,
                  size.title.use = 14) + ylab(label ="Total Expression") +
            scale_fill_manual(values = plotColours))
    dev.off()
  }
  return(seurat)
}

#### Wilcoxon test between clusters
WILCOX_DE = function(x, cellType, clusters = combined@ident){
  cells.1 = clusters == cellType
  cells.2 = !cells.1
  v1 = as.numeric(x[cells.1])
  v2 =  as.numeric(x[cells.2])
  pVal = wilcox.test(v1,v2)$p.value
  medCluster = median(v1)
  medRest = median(v2)
  LFC = medCluster  - medRest
  return(c(pVal, medCluster, medRest, LFC))
}

calcPVals = function(seurat,  colNames){
  markerSums = as.data.frame(t(seurat@meta.data[, colNames]))
  ###log expression values
  markerSums = apply(markerSums, 1, function(x){
    minVal = min(x)
    x = x + abs(minVal) + 0.01
    x = log(x)
    return(x)
  })




  res = lapply(levels(seurat@ident), function(x)
  {tmp = apply(markerSums, 2, WILCOX_DE, cellType = x, clusters = seurat@ident)}
  )
  names(res) = levels(seurat@ident)
  res = lapply(res, function(DF) {rownames(DF) <-  c("pVal", "log_meanCluster", "log_meanRest", "LFC");
  DF})
  resDF = data.frame(t(data.frame(res)))
  resDF$p.adjust= p.adjust(resDF$pVal)
  resDF$cellType = rep(levels(seurat@ident), each = 2)

  rNames  = paste(rep(levels(seurat@ident), each = 2),
                  rep(c("Early", "Late"), length(levels(seurat@ident))),
                  sep = " ")

  rownames(resDF) = rNames

  return(resDF)
}



##### convert seurat to sce for sc3 and scmap
convertToSCE = function(seurat){
  sce = as.SingleCellExperiment(seurat)
  colData(sce)$cell_type1 = seurat@ident
  colData(sce)$feature_symbol = rownames(sce)
  return(sce)}



######## SEURAT FUNCTIONS #########
LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0,
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t,
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 +
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05,
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t,
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05,
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t,
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}


###### CONVERT CDS to SINGLE CELL EXPERIMENT ######
exportCDS = function(cds){
  sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(cds@assayData$exprs)),
    colData = pData(cds),
    rowData = fData(cds)
  )
  return(sce)
}


