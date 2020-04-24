seurat_v3_hvg <- function(x, nfeatures = 5000) {
  hvg.df <- Seurat::FindVariableFeatures(x, selection.method = "vst",
                                         verbose = FALSE)
  return(order(hvg.df$vst.variance.standardized, decreasing = TRUE)[1:nfeatures])
}

