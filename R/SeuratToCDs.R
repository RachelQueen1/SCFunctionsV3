#' Convert Seurat object to CDs for Monocle 3
#'
#' @param sObj 
#'
#' @return cds object
#' @export
#'
#' @examples cds <- SeuratToCDs(sobj)

SeuratToCDs <- function(sObj){
    genesUse <- rownames(sObj)
    expression_matrix = sObj@assays$RNA@counts[genesUse, ]
    cell_metadata = sObj@meta.data[colnames(expression_matrix), ]
    
    gene_annotation = data.frame(gene_short_name = rownames(expression_matrix))
    rownames(gene_annotation) = rownames(expression_matrix)
    
    cds <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)
    return(cds)
}