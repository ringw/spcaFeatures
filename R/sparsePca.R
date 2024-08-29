runif_uint64 <- function() {
  # Julia will accept an UInt64 random seed.
  julia_seed = do.call(
    paste0,
    append(
      list("0x"),
      as.character(
        as.hexmode(
          runif(2, max = bitwShiftL(1, 30) * 4)
          + bitwShiftL(-1, 30) * 2)
      )
    )
  )
}

#' @importFrom Seurat GetAssayData
#' @importFrom Seurat VariableFeatures
#'
seurat_spca <- function(seurat, assay='RNA', nfeatures=1000, do.correct.elbow = F, ...) {
  if (exists("LayerData", where="package:Seurat"))
    scale.data <- Seurat::LayerData(seurat, assay=assay, layer="scale.data")
  else
    scale.data <- GetAssayData(scale.data, assay=assay, slot="scale.data")
  # Suppress choosing head of VariableFeatures by passing NA, NULL, or 0.
  if (is.numeric(nfeatures) && as.logical(nfeatures))
    scale.data <- scale.data[head(VariableFeatures(seurat[[assay]]), nfeatures), ]
  covar <- tcrossprod(scale.data) / (nrow(scale.data) - 1)
  psd_spca_feature_loadings(covar, ...) %>%
    seurat_spca_from_feature_loadings(seurat, assay, do.correct.elbow)
}

psd_spca_feature_loadings <- function(
  covar, varnum = 5, npcs = 50, search_cap = 500000, eigen_gap = 0.01
) {
  julia_seed = runif_uint64()
  feature_loadings = run_optimal_spca(
    covar, K=varnum, D=npcs, search_cap=search_cap, eigen_gap=eigen_gap,
    uint64_seed=julia_seed, cgroup=cgroup
  )
  # Fix the signs of feature loadings using the median sign.
  feature_loadings_heatmap = feature_loadings %>% apply(1, \(v) v %>% subset(. != 0)) %>% t
  feature_loadings = feature_loadings * rowMedians(sign(feature_loadings_heatmap))
  t(feature_loadings)
}

seurat_spca_from_feature_loadings <- function(
  feature_loadings, seurat, assay, do.correct.elbow,
  do.rename.features = TRUE
) {
  if (do.rename.features)
    colnames(feature_loadings) <- paste0('SPARSE_', seq(ncol(feature_loadings)))
  dataMeans = rowMeans(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  dataSds = rowSds(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  # data should be nonnegative, so we will update sd to a pseudo-sd of "1" if
  # the mean is 0.
  dataSds = dataSds %>% replace(dataMeans == 0, 1)
  scale_data_command_name <- str_glue("ScaleData.{assay}")
  if (scale_data_command_name %in% names(seurat@commands))
    scale_data_command <- seurat@commands[[scale_data_command_name]]
  else
    scale_data_command <- list(do.center=TRUE, do.scale=TRUE)
  scale_data_from_zero = (
    seurat[[assay]]@scale.data[rownames(feature_loadings), ]
    + (
      (
        if (scale_data_command$do.center)
        dataMeans
        else 0
      ) / (
        if (scale_data_command$do.scale)
        dataSds
        else 1
      )
    )
  )
  cell_embeddings = t(scale_data_from_zero) %*% feature_loadings
  stdev = colSds(cell_embeddings)
  obj = CreateDimReducObject(
    cell_embeddings,
    feature_loadings,
    stdev=stdev,
    assay=assay,
    key='SPARSE_'
  )
  if (do.correct.elbow) {
    obj.perm = order(obj@stdev, decreasing=T)
    obj.names = colnames(obj@cell.embeddings)
    # un-permuted feature loadings and their permutation
    obj@misc$stoch.correction = obj.perm
    obj@misc$stoch.feature.loadings = obj@feature.loadings
    obj@cell.embeddings = obj@cell.embeddings[, obj.perm, drop=FALSE]
    colnames(obj@cell.embeddings) = obj.names
    obj@feature.loadings = obj@feature.loadings[, obj.perm, drop=FALSE]
    colnames(obj@feature.loadings) = obj.names
    obj@stdev = obj@stdev[obj.perm]
    names(obj@stdev) = obj.names
  }
  obj
}

seurat_spca_from_feature_loadings_nocenter <- function(
  feature_loadings, seurat, assay, do.correct.elbow
) {
  colnames(feature_loadings) = paste0('SPARSE_', seq(ncol(feature_loadings)))
  dataMeans = rowMeans(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  dataSds = rowSds(
    seurat[[assay]]@data[rownames(feature_loadings), ]
  )
  # data should be nonnegative, so we will update sd to a pseudo-sd of "1" if
  # the mean is 0.
  dataSds = dataSds %>% replace(dataMeans == 0, 1)
  cell_embeddings = t(seurat[[assay]]@scale.data[rownames(feature_loadings), ]) %*% feature_loadings
  stdev = colSds(cell_embeddings)
  obj = CreateDimReducObject(
    cell_embeddings,
    feature_loadings,
    stdev=stdev,
    assay=assay,
    key='SPARSE_'
  )
  if (do.correct.elbow) {
    obj.perm = order(obj@stdev, decreasing=T)
    obj.names = colnames(obj@cell.embeddings)
    # un-permuted feature loadings
    obj@misc$search.feature.loadings = obj@feature.loadings
    obj@cell.embeddings = obj@cell.embeddings[, obj.perm]
    colnames(obj@cell.embeddings) = obj.names
    obj@feature.loadings = obj@feature.loadings[, obj.perm]
    colnames(obj@feature.loadings) = obj.names
    obj@stdev = obj@stdev[obj.perm]
    names(obj@stdev) = obj.names
  }
  obj
}

# Function on Seurat object. To be released later in an spcaFeatures library.
RunSparsePCA <- function(seurat, ...) {
  more_args <- list(...)
  if (!("matrix_name" %in% names(more_args))) {
    # User did not create a "covar" matrix so we will create it.
    scale.data <- seurat[[
      if (is.character(more_args$assay)) more_args$assay else "RNA"
    ]]@scale.data
    if ("nfeatures" %in% names(more_args)) {
      scale.data <- scale.data[head(VariableFeatures(seurat), more_args$nfeatures), ]
      more_args <- more_args[names(more_args) != "nfeatures"]
    }
    seurat@misc$covar <- tcrossprod(scale.data) / (nrow(scale.data) - 1)
    more_args <- more_args %>% append(list(matrix_name="covar"))
  }
  seurat[['spca']] = do.call(seurat_spca, append(list(seurat), more_args))
  seurat
}
