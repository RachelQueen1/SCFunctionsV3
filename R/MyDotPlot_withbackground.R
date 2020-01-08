#' My Dot Plot function for plotting a panel f marker genes calculated by Seurat. 
#' More options for plotting different cluster names
#' 
#' @param sObj seurat object
#' @param marker_panel Dataframe with a column called "Cell_Type" and a column  called "Gene"
#' @param assay Which assay to use. Default is "RNA"
#' @param expression_data Which expression data to use. Default is "scale.data". Could use counts
#' @param dot.scale maximum dot size
#' @param scale.by 'size' or 'radius'. Default is "radius"
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' clusters = seuratObj@active.ident
#' MyDotPlot(sObj = seurat_filtered, marker_panel = retina_sanes



MyDotPlot <- function(sObj, 
                       marker_panel,
                       assay = "RNA",
                       expression_data = "scale.data",
                       dot.scale = 6,
                       scale.by = "radius"
                       ){
    ### get expression data
    exprs = slot(sObj, "assays")[[assay]] %>% slot(expression_data) %>% data.frame()
    
    ###binary matrix gene expression
    is.exp = slot(sObj, "assays")[["RNA"]] %>% slot("counts") %>% .[] > 0   
    is.exp =  is.exp %>% data.frame() %>%  .[] * 1
    
    ### check all genes are present in seurat object
    marker_panel = marker_panel %>% filter(marker_panel$Gene %in% rownames(exprs))
    gene.list = marker_panel$Gene
    
    
    ### extract just the marker genes
    exprs = exprs %>% 
        filter(rownames(exprs) %in% gene.list) %>% 
        t() %>% 
        data.frame() %>%
        magrittr::set_colnames(rownames(exprs)[rownames(exprs) %in% gene.list])
    
    is.exp = is.exp %>%
        filter(rownames(is.exp) %in% gene.list) %>%
        t() %>%
        data.frame() %>%
        magrittr::set_colnames(rownames(is.exp)[rownames(is.exp) %in% gene.list])
    
    ### add clusters
    exprs$Clusters <- is.exp$Clusters <- sObj@active.ident
    
    
    df_gene= exprs %>% group_by(Clusters) %>% summarise_all(mean) %>% data.frame()
    df_gene[,-1] = df_gene[,-1] %>% scale() %>% MinMax(min = -2.5, max = 2.5)
    
    
    
    df_percent= is.exp %>% group_by(Clusters) %>% summarise_all(mean) %>% data.frame()
    df_percent[, 2:ncol(df_percent)] = df_percent[, 2:ncol(df_percent)] * 100
    # percentages <- df_percent[,-1]
    # percentages[percentages == 0] <- NA
    # df_percent[,-1] <-percentages
    
    ## update gene list
    gene.list <- gene.list %>% gsub("-", ".", .)
    ## reorder data frames
    df_gene <- df_gene[, c("Clusters", gene.list)]
    df_percent <- df_percent[, c("Clusters", gene.list)]
    
    ### create ggplot df
    ggdf = reshape2::melt(df_gene)
    colnames(ggdf) = c("cluster", "gene", "expression")
    ggdf$percent = reshape2::melt(df_percent)$value  
    
    
    #ggdf<-ggdf %>% purrr::map_df(rev)
    
    ### number of clusters
    noClusters = length(unique(ggdf$cluster)) + 0.5
    
    ### add in option to supply custom colours
    
    #### Plot Colours #####
    plotCols = c(RColorBrewer::brewer.pal(8, "Pastel1"),
                 RColorBrewer::brewer.pal(7, "Pastel2"))
    
    # plotCols = c(RColorBrewer::brewer.pal(5, "Set1"),
    #              RColorBrewer::brewer.pal(5, "Dark2"),
    #              RColorBrewer::brewer.pal(5, "Dark2"))
    # 
    ### add in error message if catergories greater than 15 (reduce numer of groups or use custom colours)
    
    
    marker_panel$Cell_Type <- factor(marker_panel$Cell_Type, levels = unique(marker_panel$Cell_Type))
    pLegend = ggplot(marker_panel, aes(Cell_Type, 
                                       Gene, 
                                       color = Cell_Type)) + 
        geom_point(shape = 15) + 
        scale_color_manual(values = plotCols[1:length(unique(marker_panel$Cell_Type))])+
        theme(legend.position="bottom")+ guides(colour = guide_legend(nrow = 1))
    
    pLegend <- get_legend(pLegend)
    
    
    
    ### scale function
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    
    #### MAKE PLOT #####
    p1 = ggplot(ggdf, aes(gene, cluster)) +
        geom_point(size = 0) 
    
    
    ### add background colours
    for (ct in unique(marker_panel$Cell_Type)){
        i = which(unique(marker_panel$Cell_Type) %in% ct)
        rng = which(marker_panel$Cell_Type %in% ct)
        xStart = min(rng) -0.5
        xStop = max(rng) + 0.5
        col=plotCols[i]
        p1 = p1 + geom_rect(data=NULL,aes_string(xmin=xStart,
                                                 xmax=xStop,
                                                 ymin=-Inf,
                                                 ymax=Inf),
                            fill=col,
                            alpha =0.1) 
            
        }
    
   
    ### add points to plot
    p1 = p1 +  geom_point(aes(colour = expression, size = percent)) +
        #scale_colour_gradient(low="#f0f0f0", high="#000000") +
        scale_colour_gradient(low="#ffffff", high="#4100aa") +
        scale_alpha(range = c(0.1, 1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale.func(range = c(0, dot.scale))
    p2 = plot_grid(p1, pLegend, nrow = 2, rel_heights = c(1, 0.1))
    return(p2)
}





### to do
### make into function
### push to github

