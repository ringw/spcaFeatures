RunUMAPCenteredCosine <- function(seurat, spca, dims=1:50, seed.use=1, umap_transform=diag(2)) {
  seurat[['spca']] = spca
  seurat[['spca.centered']] <- CreateDimReducObject(
    embeddings = spca@cell.embeddings %>% scale(scale = FALSE, center = TRUE),
    key = "SPCACEN_",
    assay = "RNA"
  )
  seurat <- seurat %>%
    RunUMAP(
      reduction = "spca.centered", reduction.name = "umap.spca", dims = dims,
      seed.use = seed.use
    )
  dimnames(umap_transform) <- rep(list(colnames(seurat[['umap.spca']])), 2)
  seurat[['umap.spca']]@cell.embeddings <- (
    seurat[['umap.spca']]@cell.embeddings %*% umap_transform
  )
  seurat
}