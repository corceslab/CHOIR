# ---------------------------------------------------------------------------
# Helper functions for constructing the hierarchical clustering tree
# ---------------------------------------------------------------------------

# Run dimensionality reduction ---------------------------
#
# Generates a dimensionality reduction using the supplied object and method.
# Recursive for multi-modal data, using vectors for input parameter values.
#
# Returns a list containing the cell embeddings coordinates of the reduction
# ('reduction_coords'), a vector containing the variable feature names
# ('var_features'), and if 'return_full' = TRUE, the full dimensionality
# reduction with all associated metadata ('full_reduction').
#
# Parameters:
# object -- An object of class Seurat, SingleCellExperiment, or ArchRProject
# normalization_method -- String or vector indicating which normalization method to use
# reduction_method -- String or vector indicating which dimensionality reduction method to use
# reduction_params -- A list of additional parameters to be passed to the selected dimensionality reduction method
# n_var_features -- A numerical value indicating how many variable features to identify
# batch_correction_method -- String or vector indicating which batch correction method to use
# batch_correction_params -- A list of additional parameters to be passed to the selected batch correction method
# batch_labels -- A string indicating the column containing the batch labels
# use_assay -- For Seurat or SingleCellExperiment objects, a string or vector indicating the assay(s) to use in the provided object
# use_slot -- For Seurat objects, a string or vector indicating the slot/layer to use in the provided object
# ArchR_matrix -- For ArchR objects, a string or vector indicating which matrix to use in the provided object
# ArchR_depthcol -- For ArchR objects, a string or vector indicating which column to use for correlation with sequencing depth
# atac -- A boolean value indicating whether the provided data is ATAC-seq data
# use_cells -- A vector of cell names/IDs to subset the object by, prior to dimensionality reduction
# return_full -- A boolean value indicating whether to return the full dimensionality reduction with all associated metadata
# n_cores -- A numeric value indicating the number of cores to use for parallelization
# random_seed -- A numeric value indicating the random seed used
# verbose -- A boolean value indicating whether to use verbose output during the execution of this function

.runDimReduction <- function(object,
                             normalization_method,
                             reduction_method,
                             reduction_params,
                             n_var_features,
                             batch_correction_method,
                             batch_correction_params,
                             batch_labels,
                             use_assay,
                             use_slot,
                             ArchR_matrix,
                             ArchR_depthcol,
                             atac,
                             use_cells = NULL,
                             return_full = FALSE,
                             n_cores,
                             random_seed,
                             verbose) {

  # Set full reduction to NULL by default
  full_reduction <- NULL

  ## If there are multiple data modalities from the same cells (e.g., RNA + ATAC)
  if (methods::is(object, "ArchRProject")) {
    n_modalities <- length(ArchR_matrix)
  } else {
    n_modalities <- length(use_assay)
  }

  if (n_modalities > 1) {
    # Initiate lists
    var_features_list <- vector("list", n_modalities)
    reduction_coords_list <- vector("list", n_modalities)
    reduction_full_list <- vector("list", n_modalities)
    # Run dimensionality reduction recursively on each
    for (i in 1:n_modalities) {
      # Match input arguments
      normalization_method_i <- .matchArg(normalization_method, i)
      reduction_method_i <- .matchArg(reduction_method, i)
      n_var_features_i <- .matchArg(n_var_features, i)
      batch_correction_method_i <- .matchArg(batch_correction_method, i)
      use_assay_i <- .matchArg(use_assay, i)
      use_slot_i <- .matchArg(use_slot, i)
      ArchR_matrix_i <- .matchArg(ArchR_matrix, i)
      ArchR_depthcol_i <- .matchArg(ArchR_depthcol, i)
      atac_i <- .matchArg(atac, i)

      if (length(reduction_params) > 0) {
        if (methods::is(object, "ArchRProject") & names(reduction_params)[i] == ArchR_matrix_i) {
          reduction_params_i <- reduction_params[[i]]
        } else if (!methods::is(object, "ArchRProject") & names(reduction_params)[i] == use_assay_i) {
          reduction_params_i <- reduction_params[[i]]
        } else {
          reduction_params_i <- reduction_params
        }
      } else {
        reduction_params_i <- reduction_params
      }

      if (length(batch_correction_params) > 0) {
        if (methods::is(object, "ArchRProject") & names(batch_correction_params)[i] == ArchR_matrix_i) {
          batch_correction_params_i <- batch_correction_params[[i]]
        } else if (!methods::is(object, "ArchRProject") & names(batch_correction_params)[i] == use_assay_i) {
          batch_correction_params_i <- batch_correction_params[[i]]
        } else {
          batch_correction_params_i <- batch_correction_params
        }
      } else {
        batch_correction_params_i <- batch_correction_params
      }

      # Run dimensionality reduction
      reduction_output_i <- .runDimReduction(object = object,
                                             normalization_method = normalization_method_i,
                                             reduction_method = reduction_method_i,
                                             reduction_params = reduction_params_i,
                                             n_var_features = n_var_features_i,
                                             batch_correction_method = batch_correction_method_i,
                                             batch_correction_params = batch_correction_params_i,
                                             batch_labels = batch_labels,
                                             use_assay = use_assay_i,
                                             use_slot = use_slot_i,
                                             ArchR_matrix = ArchR_matrix_i,
                                             ArchR_depthcol = ArchR_depthcol_i,
                                             atac = atac_i,
                                             use_cells = use_cells,
                                             return_full = methods::is(object, "ArchRProject"),
                                             n_cores = n_cores,
                                             random_seed = random_seed,
                                             verbose = verbose)
      var_features_list[[i]] <- reduction_output_i[["var_features"]]
      reduction_coords_list[[i]] <- reduction_output_i[["reduction_coords"]]
      reduction_full_list[[i]] <- reduction_output_i[["full_reduction"]]
    }
    # For ArchR objects, combine reductions
    if (methods::is(object, "ArchRProject")) {
      for (i in 1:length(reduction_full_list)) {
        object@reducedDims[[paste0("DR_", i)]] <- reduction_full_list[[i]]
      }
      object <- ArchR::addCombinedDims(object, reducedDims = paste0("DR_", seq(1:length(reduction_full_list))), name =  "LSI_Combined")
      if (return_full == TRUE) {
        full_reduction <- object@reducedDims$LSI_Combined
      }
      reduction_output <- list("reduction_coords" = object@reducedDims$LSI_Combined$matRD,
                               "var_features" = var_features_list,
                               "full_reduction" = full_reduction)
    } else {
      # Create output list
      reduction_output <- list("reduction_coords" = reduction_coords_list,
                               "var_features" = var_features_list,
                               "full_reduction" = full_reduction)
    }

  } else {

    ## Run dimensionality reduction based on object type and selected method

    # Number of variable features to use
    if (is.null(n_var_features)) {
      if (atac == FALSE) {
        n_var_features <- 2000
      } else if (atac == TRUE) {
        n_var_features <- 25000
      }
    }
    # For Seurat and SingleCellExperiment objects
    if (methods::is(object, "Seurat") | methods::is(object, "SingleCellExperiment")) {
      # Extract normalized feature x cell matrix
      feature_matrix <- .getMatrix(object,
                                   use_assay = use_assay,
                                   use_slot = use_slot,
                                   use_cells = use_cells,
                                   verbose = verbose)
      n_cells <- ncol(feature_matrix)
      # If batch correction, extract metadata
      if (batch_correction_method != "none") {
        if (methods::is(object, "Seurat")) {
          metadata <- object@meta.data
        } else if (methods::is(object, "SingleCellExperiment")) {
          metadata <- object@colData
        }
        if (!is.null(use_cells)) {
          metadata <- metadata[use_cells,]
        }
        metadata$CellID <- rownames(metadata)
        metadata <- metadata[match(colnames(feature_matrix), metadata$CellID), ]
        if (!identical(rownames(metadata), colnames(feature_matrix))) {
          stop("Cannot match cell IDs to batch correction group.")
        }
      }
      # Non-ATAC vs ATAC data
      if (atac == FALSE) {
        # Identify variable features or run SCTransform
        if (normalization_method == "SCTransform" & (is.null(use_assay) || use_assay != "SCT")) {
          feature_matrix <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(feature_matrix),
                                                                 variable.features.n = n_var_features,
                                                                 seed.use = random_seed,
                                                                 verbose = FALSE)@assays$SCT@scale.data)
          var_features <- rownames(feature_matrix)
        } else if (normalization_method == "none" & (is.null(use_assay) || use_assay != "SCT")) {
          var_features <- Seurat::FindVariableFeatures(feature_matrix, verbose = FALSE)
          if ("vst.variance.standardized" %in% colnames(var_features)) {
            var_features <- var_features %>%
              dplyr::arrange(-vst.variance.standardized) %>%
              utils::head(n_var_features) %>%
              rownames()
          } else {
            var_features <- var_features %>%
              dplyr::arrange(-variance.standardized) %>%
              utils::head(n_var_features) %>%
              rownames()
          }
        } else if (!is.null(use_assay) && use_assay == "SCT") {
          var_features <- rownames(feature_matrix)
          if (length(var_features) != n_var_features) {
            warning("Number of variable features pulled from provided 'SCT' assay does not match provided input to parameter 'n_var_features'.")
          }
        }

        # If no method is provided, run with default method
        if (is.null(reduction_method)) {
          reduction_method <- "PCA"
        } else {
          # Check if provided method is among allowable methods
          if (!(reduction_method %in% c("PCA", "LSI"))) {
            stop("Allowable dimensionality reduction methods for Seurat and SingleCellExperiment objects include 'PCA' and 'LSI'. Please supply valid input!")
          }
        }
        # By method
        if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Running ", reduction_method, " with ", n_var_features, " variable features..")
        if (reduction_method == "PCA") {

          if (normalization_method == "SCTransform" | (!is.null(use_assay) && use_assay == "SCT")) {
            scaled_features <- feature_matrix[var_features,]
          } else {
            scaled_features <- suppressWarnings(Seurat::ScaleData(feature_matrix[var_features,], verbose = FALSE))
          }

          # If parameters are provided, check 'em
          if (any(names(reduction_params) %in% c("object", "assay", "features", "seed.use"))) {
            stop("Parameter inputs for 'reduction_params' conflict with parameters that are necessarily set by CHOIR. Please supply valid input!")
          }
          # If provided parameters need to be supplemented, do so
          if (!any(names(reduction_params) == "verbose")) {
            reduction_params$verbose <- FALSE
          }
          if (!any(names(reduction_params) == "npcs")) {
            reduction_params$npcs <- min(50, n_cells-1)
          } else {
            reduction_params$npcs <- min(reduction_params$npcs, n_cells-1)
          }
          # Run PCA
          reduction_coords <- suppressWarnings(do.call(Seurat::RunPCA, c(list("object" = scaled_features,
                                                                              "assay" = "CHOIR",
                                                                              "features" = var_features,
                                                                              "seed.use" = random_seed),
                                                                         reduction_params))@cell.embeddings)
        } else if (reduction_method == "LSI") {

          if (normalization_method == "SCTransform" | (!is.null(use_assay) && use_assay == "SCT")) {
            scaled_features <- feature_matrix[var_features,]
          } else {
            scaled_features <- suppressWarnings(Seurat::ScaleData(feature_matrix[var_features,], verbose = FALSE))
          }

          # If parameters are provided, check 'em
          if (any(names(reduction_params) %in% c("object", "assay", "features"))) {
            stop("Parameter inputs for 'reduction_params' conflict with parameters that are necessarily set by CHOIR. Please supply valid input!")
          }
          # If provided parameters need to be supplemented, do so
          if (!any(names(reduction_params) == "verbose")) {
            reduction_params$verbose <- FALSE
          }
          .requirePackage("Signac", source = "cran")
          reduction_coords <- do.call(Signac::RunSVD, c(list("object" = scaled_features[var_features,]),
                                                        reduction_params))@cell.embeddings
        }
      } else if (atac == TRUE) {
        # Find "top" features instead of variable features
        .requirePackage("Signac", source = "cran")
        var_features <- Signac::FindTopFeatures(feature_matrix, verbose = FALSE) %>%
          dplyr::arrange(-percentile) %>%
          utils::head(n_var_features) %>%
          rownames()
        # If no method is provided, run with default method
        if (is.null(reduction_method)) {
          reduction_method <- "LSI"
        } else {
          # Check if provided method is among allowable methods
          if (!(reduction_method %in% c("LSI"))) {
            stop("Allowable dimensionality reduction methods for Seurat and SingleCellExperiment objects containing ATAC-seq data include 'LSI' only. Please supply valid input!")
          }
        }
        # By method
        if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Running ", reduction_method, " for ATAC-seq data with ", n_var_features, " variable features..")
        if (reduction_method == "LSI") {
          # If parameters are provided, check 'em
          if (any(names(reduction_params) %in% c("object", "assay", "features"))) {
            stop("Parameter inputs for 'reduction_params' conflict with parameters that are necessarily set by CHOIR. Please supply valid input!")
          }
          # If provided parameters need to be supplemented, do so
          if (!any(names(reduction_params) == "verbose")) {
            reduction_params$verbose <- FALSE
          }
          reduction_coords <- do.call(Signac::RunSVD, c(list("object" = feature_matrix[var_features,]),
                                                        reduction_params))@cell.embeddings
        }
      }

      # Harmony batch correction
      if (batch_correction_method == "Harmony") {
        # Check number of batches
        n_batches <- dplyr::n_distinct(metadata[,batch_labels])
        if (n_batches > 1) {
          # Check provided parameters
          if (any(names(batch_correction_params) %in% c("data_mat", "meta_data", "vars_use", "do_pca"))) {
            stop("Parameter inputs for 'batch_correction_params' conflict with parameters that are necessarily set by CHOIR. Please supply valid input!")
          }
          # If provided parameters need to be supplemented, do so
          if (!any(names(batch_correction_params) == "verbose")) {
            batch_correction_params$verbose <- FALSE
          }
          # Run Harmony
          if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Running Harmony batch correction using column '", batch_labels, "'..")
          if (!("character" %in% methods::is(metadata[,batch_labels]))) {
            metadata[,batch_labels] <- as.character(metadata[,batch_labels])
          }
          reduction_coords <- suppressWarnings(do.call(harmony::HarmonyMatrix, c(list("data_mat" = reduction_coords,
                                                                                      "meta_data" = metadata,
                                                                                      "vars_use" = batch_labels),
                                                                                 batch_correction_params)))
        } else {
          message(format(Sys.time(), "%Y-%m-%d %X"), " : Only one batch present. Skipped Harmony batch correction.")
        }
      }
    } else if (methods::is(object, "ArchRProject")) {
      # Subset object if cell names are provided
      if (!is.null(use_cells)) {
        object <- object[use_cells,]
      }
      n_cells <- nrow(object@cellColData)
      # If no matrix name and/or depthcol is provided
      if (atac == FALSE) {
        if (is.null(ArchR_matrix)) {
          ArchR_matrix <- "GeneExpressionMatrix"
        }
        if (is.null(ArchR_depthcol)) {
          ArchR_depthcol <- "Gex_nUMI"
        }
      } else if (atac == TRUE) {
        if (is.null(ArchR_matrix)) {
          ArchR_matrix <- "TileMatrix"
        }
        if (is.null(ArchR_depthcol)) {
          ArchR_depthcol <- "nFrags"
        }
      }

      # If no method is provided, run with default method
      if (is.null(reduction_method)) {
        reduction_method <- "IterativeLSI"
      } else {
        # Check if provided method is among allowable methods
        if (reduction_method != "IterativeLSI") {
          stop("Allowable dimensionality reduction methods for ArchRProject objects include 'IterativeLSI'. Please supply valid input!")
        }
      }
      # By method
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Running ", reduction_method, " with ", n_var_features, " variable features..")
      if (reduction_method == "IterativeLSI") {
        # Check provided parameters
        if (any(names(reduction_params) %in% c("ArchRProj", "name", "varFeatures", "saveIterations", "useMatrix", "depthCol", "force", "threads", "seed"))) {
          stop("Parameter inputs for 'reduction_params' conflict with parameters that are necessarily set by CHOIR. Please supply valid input!")
        }
        # If provided parameters need to be supplemented, do so
        if (!any(names(reduction_params) == "verbose")) {
          reduction_params$verbose <- FALSE
        }
        if(!any(names(reduction_params) == "dimsToUse")) {
          reduction_params$dimsToUse <- 1:min(30, n_cells - 1)
        }
        # Run iterative LSI
        if (ArchR_matrix == "GeneScoreMatrix") {
          iterativeLSI_matrix <- "TileMatrix"
        } else {
          iterativeLSI_matrix <- ArchR_matrix
        }

        if (ArchR_matrix == "GeneExpressionMatrix") {
          object <- do.call(ArchR::addIterativeLSI, c(list("ArchRProj" = object,
                                                           "name" = "CHOIR_IterativeLSI",
                                                           "varFeatures" = 2000,
                                                           "saveIterations" = FALSE,
                                                           "useMatrix" = iterativeLSI_matrix,
                                                           "depthCol" = ArchR_depthcol,
                                                           "force" = TRUE,
                                                           "seed" = random_seed,
                                                           "threads" = n_cores),
                                                      reduction_params))
        } else {
          object <- do.call(ArchR::addIterativeLSI, c(list("ArchRProj" = object,
                                                           "name" = "CHOIR_IterativeLSI",
                                                           "varFeatures" = 25000,
                                                           "saveIterations" = FALSE,
                                                           "useMatrix" = iterativeLSI_matrix,
                                                           "depthCol" = ArchR_depthcol,
                                                           "force" = TRUE,
                                                           "seed" = random_seed,
                                                           "threads" = n_cores),
                                                      reduction_params))
        }

        # Extract variable features (dataframe)
        if (ArchR_matrix != "GeneScoreMatrix") {
          var_features <- object@reducedDims$CHOIR_IterativeLSI$LSIFeatures
        } else {
          # Extract GeneScoreMatrix ### FIX LATER ###
          feature_matrix <- ArchR::getMatrixFromProject(object, useMatrix = "GeneScoreMatrix")
          print("check1")
          feature_names <- feature_matrix@elementMetadata$name
          print("check2")
          feature_matrix <- feature_matrix@assays@data$GeneScoreMatrix
          print("check3")
          rownames(feature_matrix) <- feature_names
          print("check4")
          # Subset to current cells
          if (!is.null(use_cells)) {
            print(use_cells[1:5])
            print(colnames(feature_matrix)[1:5])
            feature_matrix <- feature_matrix[,use_cells]
          }
          print("check5")
          feature_matrix <- as.matrix(feature_matrix)
          print(feature_matrix[1:5,1:5])
          print("check6")
          # Find variable features
          var_features <-  Seurat:::FindVariableFeatures.V3Matrix(feature_matrix, verbose = TRUE)
          print("check7")
          print(head(var_features))
          var_features <- data.frame(var_features)
          print("check8")
          if ("vst.variance.standardized" %in% colnames(var_features)) {
            var_features <- var_features %>%
              dplyr::arrange(-vst.variance.standardized) %>%
              utils::head(n_var_features) %>%
              rownames()
          } else {
            var_features <- var_features %>%
              dplyr::arrange(-variance.standardized) %>%
              utils::head(n_var_features) %>%
              rownames()
          }
          print("check9")
        }

        # Harmony batch correction
        # Check number of batches
        if (batch_correction_method == "Harmony") {
          n_batches <- dplyr::n_distinct(object@cellColData[, batch_labels])
        } else {
          n_batches <- 1
        }
        if (batch_correction_method == "Harmony" & n_batches > 1) {
          # Check provided parameters
          if (any(names(batch_correction_params) %in% c("ArchRProj", "reducedDims", "name", "groupBy", "force"))) {
            stop("Parameter inputs for 'batch_correction_params' conflict with parameters that are necessarily set by CHOIR. Please supply valid input!")
          }
          # If provided parameters need to be supplemented, do so
          if (!any(names(batch_correction_params) == "verbose")) {
            batch_correction_params$verbose <- FALSE
          }
          if (!("character" %in% methods::is(object@cellColData[, batch_labels]))) {
            object@cellColData[, batch_labels] <- as.character(object@cellColData[, batch_labels])
          }
          if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Running Harmony batch correction using column '", batch_labels, "'..")
          object <- do.call(ArchR::addHarmony, c(list("ArchRProj" = object,
                                                      "reducedDims" = "CHOIR_IterativeLSI",
                                                      "name" = "CHOIR_Harmony",
                                                      "groupBy" = batch_labels,
                                                      "force" = TRUE),
                                                 batch_correction_params))

          # Extract dimensionality reduction coordinates
          reduction_coords <- object@reducedDims$CHOIR_Harmony$matDR
          # If part of multiple modalities, also extract full ArchR reduction
          if (return_full == TRUE) {
            full_reduction <- object@reducedDims$CHOIR_Harmony
          }
        } else {
          if (batch_correction_method == "Harmony") {
            message(format(Sys.time(), "%Y-%m-%d %X"), " : Only one batch present. Skipped Harmony batch correction.")
          }
          # Extract dimensionality reduction coordinates
          reduction_coords <- object@reducedDims$CHOIR_IterativeLSI$matSVD
          # If part of multiple modalities, also extract full ArchR reduction
          if (return_full == TRUE) {
            full_reduction <- object@reducedDims$CHOIR_IterativeLSI
          }
        }
      }
    }
    # Check whether reduction_coords have row and column names
    if (is.null(rownames(reduction_coords))) {
      if (!is.null(use_cells)) {
        rownames(reduction_coords) <- use_cells
      } else {
        rownames(reduction_coords) <- .getCellIDs(object = object, use_assay = use_assay)
      }
    }
    if (is.null(colnames(reduction_coords))) {
      colnames(reduction_coords) <- paste0(reduction_method, "_", seq(1, ncol(reduction_coords)))
    }
    # Output
    reduction_output <- list("reduction_coords" = reduction_coords,
                             "var_features" = var_features,
                             "full_reduction" = full_reduction)
  }
  # Returns a list containing:
  # - Dimensionality reduction coordinates
  # - Variable features
  # - Full reduction if applicable
  return(reduction_output)
}

# Get multi-modal distance ---------------------------
#
# Generate a distance matrix from multi-modal data
#
# object -- Either a Seurat or SingleCellExperiment object with multiple modalities
# reduction_list -- List of dimensionality reductions for each modality
# dim_list -- List indicating which dimensions to use for each reduction
.getMultiModalDistance <- function(object,
                                   reduction_list,
                                   dim_list) {
  # Number of cells
  n_cells <- ncol(object)
  # Modality weights
  modality_weights <- Seurat:::FindModalityWeights(
    object = object,
    reduction.list = reduction_list,
    dims.list = dim_list)
  # NN object
  weighted_nn <- Seurat:::MultiModalNN(
    object = object,
    modality.weight = modality_weights,
    k.nn = n_cells - 1,
    knn.range = n_cells - 1)
  # Distance matrix
  dist_mat <- weighted_nn@nn.dist
  dist_mat <- cbind(dist_mat, matrix(data = NA, nrow = n_cells, ncol = 1))
  # Index matrix
  idx_mat <- weighted_nn@nn.idx
  idx_mat <- cbind(idx_mat, matrix(data = seq(1, n_cells, 1), nrow = n_cells, ncol = 1))
  # Reorder
  dist_list <- lapply(seq(1:n_cells),
                      FUN = function(x) {
                        dist_mat[x,] <- dist_mat[x, order(idx_mat[x,])]
                      })
  distance_matrix <- do.call(rbind, dist_list)
  rownames(distance_matrix) <- weighted_nn@cell.names
  colnames(distance_matrix) <- weighted_nn@cell.names
  # Return distance matrix
  return(distance_matrix)
}

# Get starting resolution ---------------------------
#
# Find the minimum resolution with more clusters than resolution = 0
#
# snn_matrix -- A shared nearest neighbor adjacency matrix
# cluster_params -- A list of additional parameters to be passed to Seurat::FindClusters()
# random_seed -- A numeric value indicating the random seed used
# verbose -- A boolean value indicating whether to use verbose output during the execution of this function
.getStartingResolution <- function(snn_matrix,
                                   cluster_params = cluster_params,
                                   random_seed = 1,
                                   verbose = TRUE) {

  # Find number of cells
  n_cells <- ncol(snn_matrix)

  # If using Leiden & data set is large, use "igraph"
  if (cluster_params$algorithm == 4 & !any(names(cluster_params) == "method") & n_cells > 30000) {
    cluster_params$method <- "igraph"
  }

  # Find initial cluster results
  res0_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix,
                                                                         "resolution" = 0,
                                                                         "random.seed" = random_seed),
                                                                    cluster_params)))
  # If singletons are grouped
  if (cluster_params$group.singletons == TRUE) {
    res0_clusters[,1] <- as.numeric(as.factor(res0_clusters[,1]))
  } else {
    res0_clusters <- .relabelSingletons(res0_clusters)
  }
  n_clust_res0 <- dplyr::n_distinct(res0_clusters[,1])

  try(res1_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix,
                                                                             "resolution" = 1,
                                                                             "random.seed" = random_seed),
                                                                        cluster_params))), silent = TRUE)
  if (!exists("res1_clusters")) {
    starting_res <- 0.1
    decimal_places <- 1
  } else {
    # If singletons are grouped
    if (cluster_params$group.singletons == TRUE) {
      res0_clusters[,1] <- as.numeric(as.factor(res0_clusters[,1]))
      res1_clusters[,1] <- as.numeric(as.factor(res1_clusters[,1]))
    } else {
      res0_clusters <- .relabelSingletons(res0_clusters)
      res1_clusters <- .relabelSingletons(res1_clusters)
    }
    n_clust_res1 <- dplyr::n_distinct(res1_clusters[,1])

    # Number of decimal places to use for rounding
    decimal_places <- 1

    # Progress bar
    pb <- progress::progress_bar$new(format = "                      [[ Current tree: :current iterations in :elapsed ]] ",
                                     total = NA)
    pb$tick(0)
    if (n_clust_res0 == n_clust_res1) {
      # From 1-4, +0.5 until there are more clusters than at resolution = 1
      # From 4-10, +1 until there are more clusters than at resolution = 1
      # From 10 & up, +5 until there are more clusters than at resolution = 1
      res <- 1.5
      new_n_clust <- n_clust_res1
      stop <- FALSE

      while (new_n_clust <= n_clust_res1 & stop == FALSE) {
        if (verbose) pb$tick()
        if (exists("new_clusters")) rm(new_clusters)
        try(new_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix,
                                                                                  "resolution" = res,
                                                                                  "random.seed" = random_seed),
                                                                             cluster_params))), silent = TRUE)
        # Stop if resolution is too high
        if (!exists("new_clusters")) {
          stop <- TRUE
          if (res < 4) {
            res <- res - 0.5
          } else if (res < 10) {
            res <- res - 1
          } else {
            res <- res - 5
          }
        } else {
          # If singletons are grouped
          if (cluster_params$group.singletons == TRUE) {
            new_clusters[,1] <- as.numeric(as.factor(new_clusters[,1]))
          } else {
            new_clusters <- .relabelSingletons(new_clusters)
          }
          new_n_clust <- dplyr::n_distinct(new_clusters[,1])
          if (new_n_clust <= n_clust_res1) {
            if (res < 4) {
              res <- res + 0.5
            } else if (res < 10) {
              res <- res + 1
            } else {
              res <- res + 5
            }
          }
        }
      }
      starting_res <- res
      if (verbose) message("                      Starting resolution: ", starting_res)
    } else if (n_clust_res1 > n_clust_res0) {
      # From 1 to 0.1, -0.3 until there is the same # of clusters as res=0
      # From 0.1 on, go down by a power of 10 until there is the same # of clusters as res=0
      res <- 0.7
      new_n_clust <- n_clust_res1
      while (new_n_clust > n_clust_res0) {
        if (verbose) pb$tick()
        new_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix,
                                                                              "resolution" = res,
                                                                              "random.seed" = random_seed),
                                                                         cluster_params)))
        # If singletons are grouped
        if (cluster_params$group.singletons == TRUE) {
          new_clusters[,1] <- as.numeric(as.factor(new_clusters[,1]))
        } else {
          new_clusters <- .relabelSingletons(new_clusters)
        }
        new_n_clust <- dplyr::n_distinct(new_clusters[,1])

        if (new_n_clust > n_clust_res0) {
          if (res > 0.1) {
            res <- res - 0.3
          } else {
            res <- res/10
            decimal_places <- decimal_places + 1
          }
        }
        res <- round(res, decimal_places)
      }
      # Our starting resolution is the resolution prior to the first one that
      # has the same # of clusters as res = 0
      if (res >= 0.1) {
        starting_res <- res + 0.3
      } else {
        starting_res <- res*10
        decimal_places <- decimal_places - 1
      }
      if (verbose) message("                      Starting resolution: ", starting_res)
    }
  }
  decimal_places <- decimal_places + 1
  # Return clusters at res = 0, starting resolution, & parameters
  return(list("starting_resolution" = starting_res,
              "res0_clusters" = res0_clusters,
              "decimal_places" = decimal_places))
}


# Relabel singleton clusters ---------------------------
#
# Will assign a unique cluster ID to each singleton cluster
#
# cluster_df -- A dataframe with 1 column containing the cluster IDs, where all singletons are labeled "singleton"
.relabelSingletons <- function(cluster_df) {
  # Relabel singletons as separate clusters
  colnames(cluster_df) <- "cluster"
  cluster_df$cluster <- as.character(cluster_df$cluster)
  cluster_df$index <- seq(1:nrow(cluster_df))
  cluster_df <- cluster_df %>%
    dplyr::mutate(cluster = ifelse(cluster== "singleton", paste0(cluster, index), cluster)) %>%
    dplyr::select(-index)

  cluster_df$cluster <- as.numeric(as.factor(cluster_df$cluster))

  # Return cluster dataframe
  return(cluster_df)
}

# Optimize tree ---------------------------
#
# Optimize the hierarchical structure of a clustering tree
#
# cluster_tree -- Dataframe containing the cluster IDs of each cell across the levels of a hierarchical clustering tree
# dist_matrix -- A cell to cell distance matrix
# reduction -- A dimensionality reduction matrix
# distance_approx -- Whether to use distance approximations
.optimizeTree <- function(cluster_tree,
                          dist_matrix = NULL,
                          reduction = NULL,
                          distance_approx = TRUE) {
  if (distance_approx == FALSE) {
    .requirePackage("clv", source = "cran")
  }
  n_levels <- ncol(cluster_tree)
  new_tree <- data.frame("CellID" = rownames(cluster_tree),
                         "new" = cluster_tree[,1])
  colnames(new_tree)[2] <- colnames(cluster_tree)[1]
  # For each cluster split, if there are > 3 resulting child clusters,
  # Create a hierarchy based on the cluster distances
  for (lvl in 2:n_levels) {
    # For each parent cluster
    unique_parent_IDs <- unique(cluster_tree[, lvl - 1])
    # Initiate list of subtrees
    subtree_list <- vector("list", length(unique_parent_IDs))
    max_columns <- 1
    for (parent in 1:length(unique_parent_IDs)) {
      parent_inds <- which(cluster_tree[, lvl - 1] == unique_parent_IDs[parent])
      subtree <- data.frame(CellID = rownames(cluster_tree[parent_inds,]),
                            original = cluster_tree[parent_inds, lvl])
      n_child_clusters <- dplyr::n_distinct(cluster_tree[parent_inds, lvl])
      # If > 3 child clusters
      if (n_child_clusters > 3) {
        # Create key for child clusters
        child_clusters <- data.frame(original = unique(cluster_tree[parent_inds, lvl]),
                                     key = paste0("c", unique(as.integer(as.numeric(as.factor(cluster_tree[parent_inds, lvl]))))))
        if (distance_approx == TRUE) {
          # Subset dimensionality reduction
          reduction_parent <- reduction[parent_inds, ]
          # Calculate centroid distances
          intercluster_distances <- stats::as.dist(.getCentroidDistance(reduction = reduction,
                                                                        clusters = paste0("c", as.integer(as.numeric(as.factor(cluster_tree[parent_inds, lvl]))))))
        } else {
          # Subset distance matrix
          dist_matrix_parent <- as.matrix(dist_matrix)[parent_inds, parent_inds]
          # Calculate distances between child clusters
          distances <- clv::cls.scatt.diss.mx(diss.mx = dist_matrix_parent,
                                              clust = as.integer(as.numeric(as.factor(cluster_tree[parent_inds, lvl]))))
          intercluster_distances <- stats::as.dist(distances$intercls.average)
        }
        # Create hierarchy
        hierarchy = stats::hclust(intercluster_distances)

        # For each level of the hierarchy, create a new column for our clustering tree
        for (k in 2:n_child_clusters) {
          k_clusters <- stats::cutree(hierarchy, k = k)
          convert <- merge(child_clusters,
                           data.frame(key = names(k_clusters), new = k_clusters),
                           by = "key")
          convert$new <- paste0("P", parent, "_", convert$new)
          colnames(convert) <- c("key", "original", paste0(colnames(cluster_tree)[lvl-1], ".", k))
          subtree <- merge(subtree, dplyr::select(convert, -key), by = "original", all.x = TRUE)
        }
        subtree <- subtree %>%
          dplyr::select(-original)
        max_columns <- max(max_columns, n_child_clusters)
      } else {
        subtree <- subtree %>%
          dplyr::mutate(new = paste0("P", parent, "_", 1)) %>%
          dplyr::select(-original)
        colnames(subtree) <- c("CellID", paste0(colnames(cluster_tree)[lvl-1], ".", 2))
      }
      subtree_list[[parent]] <- subtree
    }

    if (max_columns > 1) {
      # Create shell dataframe for subtree clusters
      subtrees <- data.frame(matrix(rep(NA, max_columns), nrow = 1, ncol = max_columns))
      colnames(subtrees) <- c("CellID", paste0(colnames(cluster_tree)[lvl-1], ".", seq(2, max_columns)))

      # Populate cluster IDs up to max_columns
      for (i in 1:length(unique_parent_IDs)) {
        subtree_i <- data.frame(subtree_list[[i]])
        if (ncol(subtree_i) < max_columns) {
          for (col in seq(ncol(subtree_i) + 1, max_columns)) {
            subtree_i[, paste0(colnames(cluster_tree)[lvl-1], ".", col)] <- subtree_i[, ncol(subtree_i)]
          }
        }
        subtrees <- rbind(subtrees, subtree_i)
      }
      # Cluster IDs to integers
      for (col in 2:ncol(subtrees)) {
        subtrees[, col] <- as.numeric(as.factor(subtrees[, col]))
      }

      # Insert into new tree
      subtrees <- subtrees[-1,]
      subtrees <- subtrees[match(new_tree$CellID, subtrees$CellID), ]
      new_tree <- cbind(new_tree, dplyr::select(subtrees, -CellID))
    }
    # Add current level to new tree
    new_tree$new <- cluster_tree[,lvl]
    colnames(new_tree)[ncol(new_tree)] <- colnames(cluster_tree)[lvl]
  }
  rownames(new_tree) <- new_tree$CellID
  new_tree <- new_tree %>% dplyr::select(-CellID)

  n_levels <- ncol(new_tree)
  # For each level of the cluster tree (except the first)
  # Check whether level is identical to preceding level
  remove <- c()
  for (lvl in 2:n_levels) {
    if (dplyr::n_distinct(new_tree[, lvl - 1]) ==
        dplyr::n_distinct(paste0(new_tree[, lvl], new_tree[, lvl - 1], sep = "..."))) {
      remove <- c(remove, lvl)
    }
  }
  # Remove duplicate levels
  if (length(remove) > 0) {
    if (length(remove) == (n_levels - 1)) {
      old_tree <- new_tree
      new_tree <- data.frame(lvl1 = old_tree[,1])
      colnames(new_tree) <- colnames(old_tree)[1]
      rownames(new_tree) <- rownames(old_tree)
    } else {
      new_tree <- new_tree[, -remove]
    }
  }
  # Return
  return(new_tree)
}

# Check hierarchy of cluster tree ---------------------------
#
# Check whether provided clustering tree is strictly hierarchical
#
# cluster_tree - Dataframe containing the cluster IDs of each cell across the levels of a hierarchical clustering tree
.checkHierarchy <- function(cluster_tree) {
  # Number of clusters at each level of clustering tree
  n_clusters = apply(cluster_tree, 2, function(x) length(unique(x[!is.na(x)])))
  # Check if levels are in order of the # of clusters
  if (is.unsorted(n_clusters)) {
    ord = order(n_clusters, decreasing = F)
    cluster_tree = cluster_tree[, ord]
    n_clusters = n_clusters[ord]
  }

  ## Construct tree & retrieve non-hierarchical cell movements
  ## Code adapted from R package mrtree (Minshi Peng et al. 2021)
  # Construct tree
  n_levels <- ncol(cluster_tree)
  tree <- NULL
  for (l in 1:(n_levels - 1)) {
    tab <- table(cluster_tree[, l], cluster_tree[, l + 1])
    rowsum <- rowSums(tab)
    colsum <- colSums(tab)
    norm <- matrix(rep(rowsum, ncol(tab)), ncol = ncol(tab)) +
      matrix(rep(colsum,nrow(tab)), nrow = nrow(tab), byrow = T)
    norm[norm == 0] <- Inf
    tab <- tab/norm
    nonzero_ind <- which(tab > 0)
    edgelist <- data.frame(start = rownames(tab)[row(tab)[nonzero_ind]],
                           end = colnames(tab)[col(tab)[nonzero_ind]],
                           count = tab[nonzero_ind], cost = Inf)
    edgelist$start <- paste(l, edgelist$start, sep = ";")
    edgelist$end <- paste(l + 1, edgelist$end, sep = ";")
    tree <- rbind(tree, edgelist)
  }
  # Retrieve non-hierarchical cell movements ("bad nodes")
  nb_in_edge_out <- stats::aggregate(tree$start, by = list(tree$end), FUN = function(x) length(unique(x)))
  bad_nodes <- nb_in_edge_out$Group.1[nb_in_edge_out$x > 1]

  # Report error if tree is not strictly hierarchical
  if (length(bad_nodes) > 0) {
    stop("The provided clustering tree is not strictly hierarchical at ",
         length(bad_nodes),
         " nodes. Please run buildTree() to generate a hierarchical clustering tree,",
         " or see package 'mrtree' for methods of reconciling non-hierarchical cluster trees.")
  }

  # Return cluster_tree (with sorted columns)
  return(cluster_tree)
}

# Check cluster tree labels ---------------------------
#
# Check and, if necessary, correct cluster labels of provided cluster tree
#
# cluster_tree -- Dataframe containing the cluster IDs of each cell across the levels of a hierarchical clustering tree
.checkClusterLabels <- function(cluster_tree) {
  # Loop across columns
  for (i in 1:ncol(cluster_tree)) {
    if (colnames(cluster_tree)[i] != paste0("L", i)) {
      # Change column name to signify level
      colnames(cluster_tree)[i] <- paste0("L", i)
    }
    if (!all(grepl(paste0("^P\\d*_L", i, "_\\d*$"), cluster_tree[, i]))) {
      if (all(grepl("^P\\d*_L", cluster_tree[, i]))) {
        warning("Detected problem in labels of clustering tree. Relabeling cluster names. If you used subtrees, this information will be lost.")
      }
      # Change cluster labels, first to just numerical
      cluster_tree[, i] <- as.numeric(as.factor(cluster_tree[, i]))
      # Now add labels for root tree ("P0") and level
      cluster_tree[, i] <- paste0("P0_L", i, "_", cluster_tree[, i])
    }
  }
  return(cluster_tree)
}


# Generate clustering tree ---------------------------
#
# Construct a clustering tree
#
# snn_matrix -- A shared nearest neighbor adjacency matrix
# nn_matrix -- A nearest neighbor adjacency matrix
# dist_matrix -- A cell to cell distance matrix
# reduction -- A dimensionality reduction matrix
# input_matrix -- Matrix containing the feature x cell data on which to train the random forest classifiers
# distance_approx -- Whether to use distance approximations
# tree_type -- A string indicating which type of tree to construct: 'silhouette', 'subtree', or 'full'
# max_clusters -- The extent to which the hierarchical clustering tree will be expanded: numerical or 'auto'
# cluster_params -- A list of additional parameters to be passed to Seurat::FindClusters()
# starting_resolution -- A numeric value indicating the resolution at which to begin clustering from getStartingResolution()
# res0_clusters -- A vector of cluster IDs at resolution = 0 from getStartingResolution()
# decimal_places -- A numeric value indicating the number of decimal places for rounding the resolution
# min_cluster_depth -- A numeric value indicating maximum cluster size at the bottom of the tree
# alpha -- A numerical value indicating the significance level used for the permutation test comparisons
# exclude_features -- A character vector indicating features that should be excluded from the input matrix
# n_iterations -- A numeric value indicating the number of iterations run for each random forest classifier comparison
# n_trees -- A numeric value indicating the number of trees in each random forest
# use_variance -- A boolean value indicating whether to use variance in permutation test
# min_accuracy -- A numeric value indicating the minimum accuracy below which clusters will be automatically merged
# min_connections -- A numeric value indicating the minimum number of nearest neighbors between two clusters for them to be considered "adjacent"
# max_repeat_errors -- A numeric value indicating the maximum number of cells that will be considered as repeated errors
# sample_max -- A numeric value indicating maximum number of cells used per cluster to train/test random forest classifier
# downsampling_rate -- A numeric value indicating the proportion of cells used per cluster to train/test random forest classifier
# batch_correction_method -- Character string or vector indicating which batch correction method to use
# batches -- Character vector of batch labels for each cell
# min_reads -- A numeric used to filter out features that do not have more than 1 read for this many cells in at least one of the clusters
# tree_records -- A dataframe comprising records from tree generation
# tree_id -- Name of tree
# n_cores -- A numeric value indicating the number of cores to use for parallelization
# random_seed -- A numeric value indicating the random seed used
.getTree <- function(snn_matrix,
                     nn_matrix = NULL,
                     dist_matrix = NULL,
                     reduction = NULL,
                     input_matrix = NULL,
                     distance_approx = TRUE,
                     tree_type = "silhouette",
                     max_clusters = NULL,
                     cluster_params,
                     starting_resolution = NULL,
                     res0_clusters = NULL,
                     decimal_places = NULL,
                     min_cluster_depth = 2000,
                     alpha = 0.05,
                     exclude_features = NULL,
                     n_iterations = 100,
                     n_trees = 50,
                     use_variance = TRUE,
                     min_accuracy = 0.5,
                     min_connections = 1,
                     max_repeat_errors = 20,
                     sample_max = Inf,
                     downsampling_rate = NULL,
                     min_reads = NULL,
                     batch_correction_method = NULL,
                     batches = NULL,
                     batch_LOO = NULL,
                     tree_records = NULL,
                     tree_id = "P0",
                     n_cores,
                     random_seed) {
  if (tree_type == "silhouette") {
    tree_output <- .getTree.silhouette(snn_matrix = snn_matrix,
                                 dist_matrix = dist_matrix,
                                 reduction = reduction,
                                 distance_approx = distance_approx,
                                 cluster_params = cluster_params,
                                 starting_resolution = starting_resolution,
                                 res0_clusters = res0_clusters,
                                 decimal_places = decimal_places,
                                 tree_records = tree_records,
                                 tree_id = tree_id,
                                 n_cores = n_cores,
                                 random_seed = random_seed)
  } else if (tree_type == "full") {
    tree_output <- .getTree.full(snn_matrix = snn_matrix,
                                 max_clusters = max_clusters,
                                 cluster_params = cluster_params,
                                 starting_resolution = starting_resolution,
                                 res0_clusters = res0_clusters,
                                 decimal_places = decimal_places,
                                 tree_records = tree_records,
                                 n_cores = n_cores,
                                 random_seed = random_seed)
  } else if (tree_type == "subtree") {
    tree_output <- .getTree.subtree(snn_matrix = snn_matrix,
                                    nn_matrix = nn_matrix,
                                    dist_matrix = dist_matrix,
                                    reduction = reduction,
                                    input_matrix = input_matrix,
                                    distance_approx = distance_approx,
                                    cluster_params = cluster_params,
                                    min_cluster_depth = min_cluster_depth,
                                    alpha = alpha,
                                    exclude_features = exclude_features,
                                    n_iterations = n_iterations,
                                    n_trees = n_trees,
                                    use_variance = use_variance,
                                    min_accuracy = min_accuracy,
                                    min_connections = min_connections,
                                    max_repeat_errors = max_repeat_errors,
                                    sample_max = sample_max,
                                    downsampling_rate = downsampling_rate,
                                    min_reads = min_reads,
                                    batch_correction_method = batch_correction_method,
                                    batches = batches,
                                    batch_LOO = batch_LOO,
                                    tree_records = tree_records,
                                    tree_id = tree_id,
                                    n_cores = n_cores,
                                    random_seed = random_seed)
  }

  return(list("cluster_tree" = tree_output[["cluster_tree"]],
              "tree_records" = tree_output[["tree_records"]]))
}

.getTree.silhouette <- function(snn_matrix,
                          dist_matrix,
                          reduction,
                          distance_approx,
                          cluster_params,
                          starting_resolution,
                          res0_clusters,
                          decimal_places,
                          tree_records,
                          tree_id = "P0",
                          n_cores,
                          random_seed) {
  # Find number of cells
  n_cells <- ncol(snn_matrix)

  # If using Leiden & data set is large, use "igraph"
  if (cluster_params$algorithm == 4 & !any(names(cluster_params) == "method") & n_cells > 30000) {
    cluster_params$method <- "igraph"
  }

  # Check whether package cluster is imported if distance_approx == FALSE
  if (distance_approx == FALSE) {
    .requirePackage("cluster", source = "cran")
  }

  # Initialize dataframe for multi-level clustering results across range of resolutions
  multi_level_clusters <- res0_clusters
  colnames(multi_level_clusters) <- "L0"
  # Set up values
  old_res <- 0
  res <- starting_resolution
  gap <- starting_resolution/2
  n_clust <- dplyr::n_distinct(res0_clusters[,1])

  # Starting silhouette score
  if (n_clust > 1) {
    if (distance_approx == TRUE) {
      new_silhouette <- mean(bluster::approxSilhouette(as.matrix(reduction), res0_clusters[,1])$width)
    } else {
      new_silhouette <- mean(unlist(cluster::silhouette(as.numeric(as.factor(res0_clusters[,1]))-1, dist_matrix)[,3]))
    }
    max_silhouette <- new_silhouette
  } else {
    new_silhouette <- -1
    max_silhouette <- -2
  }

  # Add to records
  tree_records <- rbind(tree_records, data.frame(tree_type = "silhouette",
                                                 tree_name = tree_id,
                                                 num_cells = n_cells,
                                                 resolution = 0,
                                                 num_clusters = n_clust,
                                                 silhouette = ifelse(n_clust > 1, round(new_silhouette, 5), NA),
                                                 neighbors_distance = NA,
                                                 neighbors_mean_accuracy = NA,
                                                 neighbors_var_accuracy = NA,
                                                 neighbors_percentile_accuracy = NA,
                                                 neighbors_percentile_variance = NA,
                                                 neighbors_decision = NA,
                                                 stop_branching_reason = NA))

  stop <- FALSE
  stop_reason <- "Further subclustering impossible."
  stop_counter <- 0
  level <- 0
  max_sil_lvl <- paste0("L", level)
  max_sil_res <- old_res
  gap_counter <- 0
  reset_counter <- 0

  # Progress bar
  pb <- progress::progress_bar$new(format = "                      [[ Current tree: :current iterations in :elapsed ]] ",
                                   total = NA)
  pb$tick(0)

  while (stop == FALSE) {
    pb$tick()
    if (exists("new_clusters")) rm(new_clusters)
    try(new_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix,
                                                                              "resolution" = res,
                                                                              "random.seed" = random_seed),
                                                                         cluster_params))), silent = TRUE)
    if (!exists("new_clusters")) {
      stop <- TRUE
      # Add to records
      tree_records <- rbind(tree_records, data.frame(tree_type = "silhouette",
                                                     tree_name = tree_id,
                                                     num_cells = n_cells,
                                                     resolution = res,
                                                     num_clusters = n_clust,
                                                     silhouette = NA,
                                                     neighbors_distance = NA,
                                                     neighbors_mean_accuracy = NA,
                                                     neighbors_var_accuracy = NA,
                                                     neighbors_percentile_accuracy = NA,
                                                     neighbors_percentile_variance = NA,
                                                     neighbors_decision = NA,
                                                     stop_branching_reason = stop_reason))
    } else {
      # If singletons are grouped
      if (cluster_params$group.singletons == TRUE) {
        new_clusters[,1] <- as.numeric(as.factor(new_clusters[,1]))
      } else {
        new_clusters <- .relabelSingletons(new_clusters)
      }
      # Number of new clusters
      new_n_clust <- dplyr::n_distinct(new_clusters[,1])

      # Stop prior to reaching all singletons / if max cluster size is 2
      if (max(table(new_clusters)) <= 2) {
        stop <- TRUE
        stop_reason <- paste0("Cannot subdivide further.")
        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "silhouette",
                                                       tree_name = tree_id,
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = new_n_clust,
                                                       silhouette = NA,
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = stop_reason))
      } else if (new_n_clust < n_clust) {
        # Stop if number of clusters starts decreasing
        stop <- TRUE
        stop_reason <- paste0("Number of clusters is decreasing.")
        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "silhouette",
                                                       tree_name = tree_id,
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = new_n_clust,
                                                       silhouette = NA,
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = stop_reason))
      } else if ((new_n_clust > (n_clust + 0.25*n_clust)) &
                 (new_n_clust - n_clust > 3) &
                 reset_counter < 3) {
        # If number of clusters has increased too dramatically (>25% increase)
        # Halve the gap & reset to previous resolution
        res <- res - gap
        gap <- gap/2
        gap_counter <- 0
        # Don't do this more than three times in a row
        reset_counter <- reset_counter + 1
      } else if (new_n_clust > n_clust) {
        # If number of clusters has increased sufficiently
        # Add new level to clustering tree
        level = level + 1
        multi_level_clusters[, paste0("L", level)] <- new_clusters[,1]
        # Update values for next level
        gap <- res - old_res
        gap_counter <- 0
        reset_counter <- 0
        old_res <- res
        n_clust <- new_n_clust

        # Calculate overall silhouette
        if (distance_approx == TRUE) {
          new_silhouette <- mean(bluster::approxSilhouette(as.matrix(reduction), new_clusters[,1])$width)
        } else {
          new_silhouette <- mean(unlist(cluster::silhouette(as.numeric(as.factor(new_clusters[,1]))-1, dist_matrix)[,3]))
        }

        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "silhouette",
                                                       tree_name = tree_id,
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = new_n_clust,
                                                       silhouette = round(new_silhouette, 5),
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = NA))

        # Compare silhouette to previous max
        if (new_silhouette >= max_silhouette) {
          max_silhouette <- max(max_silhouette, new_silhouette)
          max_sil_lvl <- paste0("L", level)
          max_sil_res <- res
          # Reset stop countdown
          stop_counter <- 0
        } else if (new_silhouette < max_silhouette) {
          # When silhouette has peaked, impose threshold for silhouette tree
          stop_counter <- stop_counter + 1
          # Stop 5 rounds after max silhouette
          if (stop_counter >= 5) {
            stop <- TRUE
            stop_reason <- paste0("Identified resolution with maximum silhouette.")
            # Add stop reason to records
            tree_records[nrow(tree_records), "stop_branching_reason"] <- stop_reason
          }
        }
      }
      res <- res + gap
      gap_counter <- gap_counter + 1
      if (gap_counter >= 5) {
        gap <- gap*2
        gap_counter <- 0
      }
      # Round, otherwise R has addition issues
      res <- round(res, decimal_places)
    }
  }

  # Clean up
  rm(snn_matrix)
  rm(dist_matrix)
  rm(reduction)

  # After max silhouette has been established, we'll remove any subsequent resolutions from the tree
  if (which(colnames(multi_level_clusters) == max_sil_lvl) == 1) {
    multi_level_clusters <- data.frame(col = multi_level_clusters[, 1])
    colnames(multi_level_clusters) <- max_sil_lvl
  } else {
    multi_level_clusters <- multi_level_clusters[,1:which(colnames(multi_level_clusters) == max_sil_lvl)]
  }
  ## We now have our multi-level clustering results across a range of resolutions
  # If necessary, reconcile into strictly hierarchical clustering tree using MRtree
  if ((ncol(multi_level_clusters) == 2 & dplyr::n_distinct(multi_level_clusters[,1]) == 1) | ncol(multi_level_clusters) == 1) {
    cluster_tree <- multi_level_clusters
  } else {
    multi_level_clusters <- as.matrix(multi_level_clusters)
    # Convert to numerical for MRtree
    multi_level_clusters <- apply(multi_level_clusters, 2, function(x) as.numeric(as.factor(x)))
    # Run MRtree
    suppressMessages(invisible(utils::capture.output(mrtree_output <- mrtree::mrtree(multi_level_clusters, n.cores = n_cores, consensus = FALSE,
                                                                                     augment.path = FALSE, verbose = FALSE))))
    cluster_tree <- data.frame(mrtree_output$labelmat.mrtree)
    # Add to records
    tree_records <- rbind(tree_records, data.frame(tree_type = "silhouette",
                                                   tree_name = tree_id,
                                                   num_cells = n_cells,
                                                   resolution = NA,
                                                   num_clusters = dplyr::n_distinct(cluster_tree[,ncol(cluster_tree)]),
                                                   silhouette = NA,
                                                   neighbors_distance = NA,
                                                   neighbors_mean_accuracy = NA,
                                                   neighbors_var_accuracy = NA,
                                                   neighbors_percentile_accuracy = NA,
                                                   neighbors_percentile_variance = NA,
                                                   neighbors_decision = NA,
                                                   stop_branching_reason = "Tree at maximum silhouette score after MRtree"))
  }
  # Add a level where all cells belong to the same cluster (if necessary)
  if (dplyr::n_distinct(cluster_tree[,1]) > 1) {
    cluster_tree <- cbind(data.frame(K1 = rep(0, nrow(cluster_tree))), cluster_tree)
  }
  # Add cell names as row names
  rownames(cluster_tree) <- rownames(res0_clusters)

  return(list("cluster_tree" = cluster_tree,
              "tree_records" = tree_records))
}

.getTree.full <- function(snn_matrix,
                          max_clusters,
                          cluster_params,
                          starting_resolution,
                          res0_clusters,
                          decimal_places,
                          tree_records,
                          n_cores,
                          random_seed) {
  # Find number of cells
  n_cells <- ncol(snn_matrix)

  # If using Leiden & data set is large, use "igraph"
  if (cluster_params$algorithm == 4 & !any(names(cluster_params) == "method") & n_cells > 30000) {
    cluster_params$method <- "igraph"
  }

  # Initialize dataframe for multi-level clustering results across range of resolutions
  multi_level_clusters <- res0_clusters
  colnames(multi_level_clusters) <- "L0"
  # Set up values
  old_res <- 0
  res <- starting_resolution
  gap <- starting_resolution/2
  n_clust <- dplyr::n_distinct(res0_clusters[,1])
  stop <- FALSE
  stop_reason <- "Further subclustering impossible."
  stop_counter <- 0
  level <- 0
  gap_counter <- 0
  reset_counter <- 0

  # Add to records
  tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                 tree_name = "P0",
                                                 num_cells = n_cells,
                                                 resolution = 0,
                                                 num_clusters = n_clust,
                                                 silhouette = NA,
                                                 neighbors_distance = NA,
                                                 neighbors_mean_accuracy = NA,
                                                 neighbors_var_accuracy = NA,
                                                 neighbors_percentile_accuracy = NA,
                                                 neighbors_percentile_variance = NA,
                                                 neighbors_decision = NA,
                                                 stop_branching_reason = NA))

  # Progress bar
  pb <- progress::progress_bar$new(format = "[[ Current tree: :current in :elapsed ]] ",
                                   total = NA)
  pb$tick(0)

  while (stop == FALSE) {
    pb$tick()
    if (exists("new_clusters")) rm(new_clusters)
    try(new_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix,
                                                                              "resolution" = res,
                                                                              "random.seed" = random_seed),
                                                                         cluster_params))), silent = TRUE)
    if (!exists("new_clusters")) {
      stop <- TRUE
      # Add to records
      tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                     tree_name = "P0",
                                                     num_cells = n_cells,
                                                     resolution = res,
                                                     num_clusters = n_clust,
                                                     silhouette = NA,
                                                     neighbors_distance = NA,
                                                     neighbors_mean_accuracy = NA,
                                                     neighbors_var_accuracy = NA,
                                                     neighbors_percentile_accuracy = NA,
                                                     neighbors_percentile_variance = NA,
                                                     neighbors_decision = NA,
                                                     stop_branching_reason = stop_reason))
    } else {
      # If singletons are grouped
      if (cluster_params$group.singletons == TRUE) {
        new_clusters[,1] <- as.numeric(as.factor(new_clusters[,1]))
      } else {
        new_clusters <- .relabelSingletons(new_clusters)
      }
      # Number of new clusters
      new_n_clust <- dplyr::n_distinct(new_clusters[,1])

      # Stop prior to reaching all singletons / if max cluster size is 2
      if (max(table(new_clusters)) <= 2) {
        stop <- TRUE
        stop_reason <- paste0("Cannot subdivide further.")
        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                       tree_name = "P0",
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = new_n_clust,
                                                       silhouette = NA,
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = stop_reason))
      } else if (new_n_clust < n_clust) {
        # Stop if number of clusters starts decreasing
        stop <- TRUE
        stop_reason <- paste0("Number of clusters is decreasing.")
        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                       tree_name = "P0",
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = new_n_clust,
                                                       silhouette = NA,
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = stop_reason))
      } else if ((new_n_clust > (n_clust + 0.25*n_clust)) &
                 (new_n_clust - n_clust > 3) &
                 reset_counter < 3) {
        # If number of clusters has increased too dramatically (>25% increase)
        # Halve the gap & reset to previous resolution
        res <- res - gap
        gap <- gap/2
        gap_counter <- 0
        # Don't do this more than three times in a row
        reset_counter <- reset_counter + 1
      } else if (new_n_clust > n_clust) {
        # If number of clusters has increased sufficiently
        # Add new level to clustering tree
        level = level + 1
        multi_level_clusters[, paste0("L", level)] <- new_clusters[,1]
        # Update values for next level
        gap <- res - old_res
        gap_counter <- 0
        reset_counter <- 0
        old_res <- res
        n_clust <- new_n_clust
        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                       tree_name = "P0",
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = new_n_clust,
                                                       silhouette = NA,
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = NA))
        if (new_n_clust > max_clusters) {
          stop <- TRUE
          stop_reason <- paste0("Reached maximum number of clusters.")
          # Add stop reason to records
          tree_records[nrow(tree_records), "stop_branching_reason"] <- stop_reason
        }
        res <- res + gap
        gap_counter <- gap_counter + 1
        if (gap_counter >= 5) {
          gap <- gap*2
          gap_counter <- 0
        }
        # Round, otherwise R has addition issues
        res <- round(res, decimal_places)
      }
    }
  }

  # Clean up
  rm(snn_matrix)

  ## We now have our multi-level clustering results across a range of resolutions
  # If necessary, reconcile into strictly hierarchical clustering tree using MRtree
  if ((ncol(multi_level_clusters) == 2 & dplyr::n_distinct(multi_level_clusters[,1]) == 1) | ncol(multi_level_clusters) == 1) {
    cluster_tree <- multi_level_clusters
  } else {
    multi_level_clusters <- as.matrix(multi_level_clusters)
    # Convert to numerical for MRtree
    multi_level_clusters <- apply(multi_level_clusters, 2, function(x) as.numeric(as.factor(x)))
    # Run MRtree
    suppressMessages(invisible(utils::capture.output(mrtree_output <- mrtree::mrtree(multi_level_clusters, n.cores = n_cores, consensus = FALSE,
                                                                                     augment.path = FALSE, verbose = FALSE))))
    cluster_tree <- data.frame(mrtree_output$labelmat.mrtree)
    # Add to records
    tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                   tree_name = "P0",
                                                   num_cells = n_cells,
                                                   resolution = NA,
                                                   num_clusters = dplyr::n_distinct(cluster_tree[,ncol(cluster_tree)]),
                                                   silhouette = NA,
                                                   neighbors_distance = NA,
                                                   neighbors_mean_accuracy = NA,
                                                   neighbors_var_accuracy = NA,
                                                   neighbors_percentile_accuracy = NA,
                                                   neighbors_percentile_variance = NA,
                                                   neighbors_decision = NA,
                                                   stop_branching_reason = "Full tree after MRtree."))
  }
  # Add a level where all cells belong to the same cluster (if necessary)
  if (dplyr::n_distinct(cluster_tree[,1]) > 1) {
    cluster_tree <- cbind(data.frame(K1 = rep(0, nrow(cluster_tree))), cluster_tree)
  }
  # Add cell names as row names
  rownames(cluster_tree) <- rownames(res0_clusters)

  return(list("cluster_tree" = cluster_tree,
              "tree_records" = tree_records))
}


.getTree.subtree <- function(snn_matrix,
                             nn_matrix,
                             dist_matrix,
                             reduction,
                             input_matrix,
                             distance_approx,
                             cluster_params,
                             min_cluster_depth,
                             alpha,
                             exclude_features,
                             n_iterations,
                             n_trees,
                             use_variance,
                             min_accuracy,
                             min_connections,
                             max_repeat_errors,
                             sample_max,
                             downsampling_rate,
                             min_reads,
                             batch_correction_method,
                             batches,
                             batch_LOO,
                             tree_records,
                             tree_id,
                             n_cores,
                             random_seed) {
  # Find number of cells
  n_cells <- ncol(snn_matrix)
  cell_IDs <- colnames(snn_matrix)

  # If using Leiden & data set is large, use "igraph"
  if (cluster_params$algorithm == 4 & !any(names(cluster_params) == "method") & n_cells > 30000) {
    cluster_params$method <- "igraph"
  }

  if (distance_approx == FALSE) {
    .requirePackage("clv", source = "cran")
  }

  # Data frame to gather comparison records
  all_metrics <- c('comparison', 'cluster1_size', 'cluster2_size', 'sample_size',
                   'mean_accuracy', 'var_accuracy', 'mean_errors',
                   'mean_permuted_accuracy', 'var_permuted_accuracy',
                   'percentile_accuracy', 'percentile_variance',
                   'n_repeat_errors1', 'n_repeat_errors2',
                   'mean_repeat_errors1', 'mean_repeat_errors2',
                   'mean_modified_accuracy', 'var_modified_accuracy',
                   'percentile_modified_accuracy', 'percentile_modified_variance',
                   'batches_used', 'batch_mean_accuracies', 'batch_mean_variances',
                   'batch_LOO_mean_accuracies', 'batch_LOO_var_accuracies', 'batch_LOO_mean_errors',
                   'batch_LOO_mean_permuted_accuracies', 'batch_LOO_var_permuted_accuracies',
                   'batch_LOO_percentile_accuracies', 'batch_LOO_percentile_variances',
                   'connectivity', 'root_distance', 'subtree_distance', 'time',
                   'decision')
  selected_metrics <- all_metrics[c(1:11,
                                    `if`(max_repeat_errors > 0, 12:15, NULL),
                                    `if`(max_repeat_errors > 0, 16:19, NULL),
                                    `if`(batch_correction_method == "Harmony", 20:22, NULL),
                                    `if`(batch_LOO == TRUE, 23:29, NULL),
                                    `if`(min_connections > 0, 30, NULL),
                                    31:32)]

  comparison_records <- data.frame(matrix(ncol = length(selected_metrics), nrow = 0))
  colnames(comparison_records) <- selected_metrics

  # Initialize dataframe for multi-level clustering results across range of resolutions
  multi_level_clusters <- data.frame(L0 = rep(0, n_cells))
  rownames(multi_level_clusters) <- colnames(snn_matrix)

  if (nrow(tree_records) > 0) {
    current_tree_name <- paste0(tree_id, "_", as.numeric(gsub(".*?(\\d+)$", "\\1", tree_records[nrow(tree_records), "tree_name"])) + 1)
  } else {
    current_tree_name <- tree_id
  }

  unresolved_cells <- colnames(snn_matrix)

  # Progress bar
  pb <- progress::progress_bar$new(format = "[[ Current: :percent in :elapsed ]] ",
                                   total = n_cells)
  pb$tick(0)

  # While there are cells belonging to clusters that have not been overclustered, continue to subcluster
  while (length(unresolved_cells) > 0) {
    pb$tick(0)
    # Clusters left to subcluster (clusters containing unresolved cells)
    to_subcluster <- unique(multi_level_clusters[unresolved_cells, ncol(multi_level_clusters)])

    # Subcluster first eligible cluster
    current_cells <- rownames(multi_level_clusters)[which(multi_level_clusters[, ncol(multi_level_clusters)] == to_subcluster[1])]
    # Get starting resolution
    starting_resolution_list <- .getStartingResolution(snn_matrix[current_cells, current_cells],
                                                       cluster_params = cluster_params,
                                                       random_seed = random_seed,
                                                       verbose = FALSE)
    res0_clusters <- starting_resolution_list[["res0_clusters"]]
    decimal_places <- starting_resolution_list[["decimal_places"]]

    # Set up values
    old_res <- 0
    res <- starting_resolution_list[["starting_resolution"]]
    gap <- starting_resolution_list[["starting_resolution"]]/2
    n_clust <- dplyr::n_distinct(res0_clusters[,1])
    stop <- FALSE
    stop_reason <- "Further subclustering impossible."
    level <- as.numeric(sub("L", "", colnames(multi_level_clusters)[ncol(multi_level_clusters)]))
    gap_counter <- 0
    reset_counter <- 0
    level_counter <- 0

    # If >1 cluster at res = 0, add to levels
    if (n_clust > 1) {
      level = level + 1
      multi_level_clusters[, paste0("L", level)] <- multi_level_clusters[, paste0("L", level - 1)]
      multi_level_clusters[current_cells, paste0("L", level)] <- res0_clusters[,1] + max(multi_level_clusters[, paste0("L", level - 1)]) + 1
      # Stop if any cluster is below or equal to max
      if (min(table(res0_clusters[,1])) <= min_cluster_depth) {
        stop <- TRUE
      }
    }

    # Subcluster until overclustered (or until alternative stop criteria are met)
    while (stop == FALSE) {
      pb$tick(0)
      if (exists("new_clusters")) rm(new_clusters)
      try(new_clusters <- suppressWarnings(do.call(Seurat::FindClusters, c(list("object" = snn_matrix[current_cells, current_cells],
                                                                                "resolution" = res,
                                                                                "random.seed" = random_seed),
                                                                           cluster_params))), silent = TRUE)
      if (!exists("new_clusters")) { # If clustering fails, stop
        stop <- TRUE
        # Add to records
        tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                       tree_name = current_tree_name,
                                                       num_cells = n_cells,
                                                       resolution = res,
                                                       num_clusters = n_clust,
                                                       silhouette = NA,
                                                       neighbors_distance = NA,
                                                       neighbors_mean_accuracy = NA,
                                                       neighbors_var_accuracy = NA,
                                                       neighbors_percentile_accuracy = NA,
                                                       neighbors_percentile_variance = NA,
                                                       neighbors_decision = NA,
                                                       stop_branching_reason = stop_reason))
        # Resolve cells
        unresolved_cells <- unresolved_cells[!(unresolved_cells %in% current_cells)]
        # Progress
        pb$tick(length(current_cells))
      } else {
        # If singletons are grouped
        if (cluster_params$group.singletons == TRUE) {
          new_clusters[,1] <- as.numeric(as.factor(new_clusters[,1]))
        } else {
          new_clusters <- .relabelSingletons(new_clusters)
        }
        # Number of new clusters
        new_n_clust <- dplyr::n_distinct(new_clusters[,1])
        # Stop prior to reaching all singletons / if max cluster size is 2
        if (max(table(new_clusters)) <= 2) {
          stop <- TRUE
          stop_reason <- "Cannot subdivide further."
          # Add to records
          tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                         tree_name = current_tree_name,
                                                         num_cells = n_cells,
                                                         resolution = res,
                                                         num_clusters = new_n_clust,
                                                         silhouette = NA,
                                                         neighbors_distance = NA,
                                                         neighbors_mean_accuracy = NA,
                                                         neighbors_var_accuracy = NA,
                                                         neighbors_percentile_accuracy = NA,
                                                         neighbors_percentile_variance = NA,
                                                         neighbors_decision = NA,
                                                         stop_branching_reason = stop_reason))
          # Resolve cells
          unresolved_cells <- unresolved_cells[!(unresolved_cells %in% current_cells)]
          # Progress
          pb$tick(length(current_cells))
        } else if (new_n_clust < n_clust) {
          # Stop if number of clusters starts decreasing
          stop <- TRUE
          stop_reason <- "Number of clusters is decreasing, will subcluster unresolved cells from previous level."
          # Add to records
          tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                         tree_name = current_tree_name,
                                                         num_cells = n_cells,
                                                         resolution = res,
                                                         num_clusters = n_clust,
                                                         silhouette = NA,
                                                         neighbors_distance = NA,
                                                         neighbors_mean_accuracy = NA,
                                                         neighbors_var_accuracy = NA,
                                                         neighbors_percentile_accuracy = NA,
                                                         neighbors_percentile_variance = NA,
                                                         neighbors_decision = NA,
                                                         stop_branching_reason = stop_reason))
          # Check which cells are resolved (from previous level)
          # Identify non-singleton clusters
          non_singleton_clusters <- names(table(multi_level_clusters[current_cells, ncol(multi_level_clusters)]))[table(multi_level_clusters[current_cells, ncol(multi_level_clusters)]) > 1]
          if (length(non_singleton_clusters) > 1) {
            # For each cluster, find the nearest cluster
            # Then find the maximum (i.e., the cluster with the most distant next door neighbors)
            if (distance_approx == TRUE) {
              # Use centroid linkage distance
              intercluster_distances <- .getCentroidDistance(reduction = reduction,
                                                             clusters = as.integer(as.numeric(as.factor(multi_level_clusters[current_cells, ncol(multi_level_clusters)]))))
            } else {
              # Use average linkage distance
              distances <- clv::cls.scatt.diss.mx(diss.mx = as.matrix(dist_matrix[current_cells, current_cells]),
                                                  clust = as.integer(as.numeric(as.factor(multi_level_clusters[current_cells, ncol(multi_level_clusters)]))))
              intercluster_distances <- distances$intercls.average[paste0("c", non_singleton_clusters), paste0("c", non_singleton_clusters)]
            }
            # Resolve clusters with < 4 cells
            resolve_clusters <- names(table(multi_level_clusters[current_cells, ncol(multi_level_clusters)]))[table(multi_level_clusters[current_cells, ncol(multi_level_clusters)]) < 4]
            if (length(resolve_clusters) > 0) {
              remove_cells <- rownames(multi_level_clusters[unresolved_cells, ])[multi_level_clusters[unresolved_cells,
                                                                                                      ncol(multi_level_clusters)] %in% resolve_clusters]
              unresolved_cells <- unresolved_cells[!(unresolved_cells %in% remove_cells)]
              # Progress
              pb$tick(length(remove_cells))
            }
            # Test neighbor pairs of each remaining cluster
            neighbor_pairs <- do.call(rbind, lapply(
              seq(1:nrow(intercluster_distances)),
              FUN = function(x)
                data.frame(cluster1 = as.numeric(sub("c", "", rownames(intercluster_distances)[x])),
                           cluster2 = as.numeric(sub("c", "", colnames(intercluster_distances)[intercluster_distances[x, ] ==
                                                                                                 min(intercluster_distances[x, ][intercluster_distances[x, ] > 0])])),
                           min_dist = min(intercluster_distances[x, ][intercluster_distances[x, ] > 0])))) %>%
              dplyr::mutate(comp = ifelse(cluster1 < cluster2, paste0(cluster1, "_", cluster2), paste0(cluster2, "_", cluster1))) %>%
              dplyr::group_by(comp) %>% dplyr::slice(1)
            while (nrow(neighbor_pairs) > 0) {
              # Test first neighbor pair
              cluster1_name <- neighbor_pairs$cluster1[1]
              cluster2_name <- neighbor_pairs$cluster2[1]
              cluster1_cells <- rownames(multi_level_clusters[current_cells,])[multi_level_clusters[current_cells, ncol(multi_level_clusters)] == cluster1_name]
              cluster2_cells <- rownames(multi_level_clusters[current_cells,])[multi_level_clusters[current_cells, ncol(multi_level_clusters)] == cluster2_name]
              # Check if this pair of clusters should remain split
              comparison <- .runPermutationTest(cluster1_name = "Cluster1",
                                                cluster1_cells = cluster1_cells,
                                                cluster1_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                             batches[cluster1_cells],
                                                                             NULL),
                                                cluster2_name  = "Cluster2",
                                                cluster2_cells = cluster2_cells,
                                                cluster2_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                             batches[cluster2_cells],
                                                                             NULL),
                                                alpha = alpha,
                                                n_iterations = n_iterations,
                                                n_trees = n_trees,
                                                use_variance = use_variance,
                                                min_accuracy = min_accuracy,
                                                min_connections = min_connections,
                                                max_repeat_errors = max_repeat_errors,
                                                collect_all_metrics = FALSE,
                                                sample_max = sample_max,
                                                downsampling_rate = downsampling_rate,
                                                min_reads = min_reads,
                                                input_matrix = input_matrix[current_cells, ],
                                                nn_matrix = nn_matrix[current_cells, current_cells],
                                                comparison_records = comparison_records,
                                                feature_importance_records = NULL,
                                                P0_distance = NA,
                                                P_i_distance = NA,
                                                n_cores = n_cores,
                                                random_seed = random_seed)
              # Add to records
              tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                             tree_name = current_tree_name,
                                                             num_cells = n_cells,
                                                             resolution = res,
                                                             num_clusters = n_clust,
                                                             silhouette = NA,
                                                             neighbors_distance = neighbor_pairs$min_dist[1],
                                                             neighbors_mean_accuracy = comparison[["comparison_records"]][1,"mean_accuracy"],
                                                             neighbors_var_accuracy = comparison[["comparison_records"]][1,"var_accuracy"],
                                                             neighbors_percentile_accuracy = comparison[["comparison_records"]][1,"percentile_accuracy"],
                                                             neighbors_percentile_variance = comparison[["comparison_records"]][1,"percentile_variance"],
                                                             neighbors_decision = comparison[["result"]],
                                                             stop_branching_reason = NA))
              # Remove row
              neighbor_pairs <- neighbor_pairs[-1,]
              # If merge
              if (comparison[["result"]] == "merge") {
                tree_records[nrow(tree_records), "stop_branching_reason"] <- "Current neighbor pair is sufficiently overclustered"
                # Resolve cells
                remove_cells <- c(cluster1_cells, cluster2_cells)[c(cluster1_cells, cluster2_cells) %in% unresolved_cells]

                unresolved_cells <- unresolved_cells[!(unresolved_cells %in% remove_cells)]
                # Remove neighbor pair of cluster2 if present
                neighbor_pairs <- neighbor_pairs %>% dplyr::filter(cluster1 != cluster2_name)
                # Progress
                pb$tick(length(remove_cells))
              }
            }
          } else {
            # If all but 1 cluster is a singleton, stop
            stop <- TRUE
            stop_reason <- "All but 1 cluster are singletons"
            # Add to records
            tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                           tree_name = current_tree_name,
                                                           num_cells = n_cells,
                                                           resolution = res,
                                                           num_clusters = n_clust,
                                                           silhouette = NA,
                                                           neighbors_distance = NA,
                                                           neighbors_mean_accuracy = NA,
                                                           neighbors_var_accuracy = NA,
                                                           neighbors_percentile_accuracy = NA,
                                                           neighbors_percentile_variance = NA,
                                                           neighbors_decision = NA,
                                                           stop_branching_reason = stop_reason))
            # Resolve clusters with < 4 cells
            resolve_clusters <- names(table(multi_level_clusters[current_cells, ncol(multi_level_clusters)]))[table(multi_level_clusters[current_cells, ncol(multi_level_clusters)]) < 4]
            remove_cells <- rownames(multi_level_clusters[unresolved_cells, ])[multi_level_clusters[unresolved_cells,
                                                                                                    ncol(multi_level_clusters)] %in% resolve_clusters]
            unresolved_cells <- unresolved_cells[!(unresolved_cells %in% remove_cells)]
            # Progress
            pb$tick(length(remove_cells))
          }
        } else if ((new_n_clust > (n_clust + 0.25*n_clust)) &
                   (new_n_clust - n_clust > 3) &
                   reset_counter < 3) {
          # If number of clusters has increased too dramatically (>25% increase)
          # Halve the gap & reset to previous resolution
          res <- res - gap
          gap <- gap/2
          gap_counter <- 0
          # Don't do this more than three times in a row
          reset_counter <- reset_counter + 1
        } else if (new_n_clust > n_clust) {
          # If number of clusters has increased sufficiently
          # Add new level to clustering tree
          level = level + 1
          multi_level_clusters[, paste0("L", level)] <- multi_level_clusters[, paste0("L", level - 1)]
          multi_level_clusters[current_cells, paste0("L", level)] <- new_clusters[,1] + max(multi_level_clusters[, paste0("L", level - 1)]) + 1
          # Update values for next level
          gap <- res - old_res
          gap_counter <- 0
          reset_counter <- 0
          old_res <- res
          n_clust <- new_n_clust
          # Check thresholds for stopping point
          # Identify non-singleton clusters
          non_singleton_clusters <- names(table(new_clusters[,1]))[table(new_clusters[,1]) > 1]
          if (length(non_singleton_clusters) > 1) {
            # For each cluster, find the nearest cluster
            # Then find the maximum (i.e., the cluster with the most distant next door neighbors)
            if (distance_approx == TRUE) {
              # Use centroid linkage distance
              intercluster_distances <- .getCentroidDistance(reduction = reduction,
                                                             clusters = as.integer(as.numeric(as.factor(new_clusters[,1]))))
            } else {
              # Use average linkage distance
              distances <- clv::cls.scatt.diss.mx(diss.mx = as.matrix(dist_matrix[current_cells, current_cells]),
                                                  clust = as.integer(as.numeric(as.factor(new_clusters[,1]))))
              intercluster_distances <- distances$intercls.average[paste0("c", non_singleton_clusters), paste0("c", non_singleton_clusters)]
            }
            furthest_neighbors <-
              dplyr::slice_max(do.call(rbind, lapply(
                seq(1:nrow(intercluster_distances)),
                FUN = function(x)
                  data.frame(cluster1 = as.numeric(sub("c", "", rownames(intercluster_distances)[x])),
                             cluster2 = as.numeric(sub("c", "", colnames(intercluster_distances)[intercluster_distances[x, ] ==
                                                                                                   min(intercluster_distances[x, ][intercluster_distances[x, ] > 0])])),
                             min_dist = min(intercluster_distances[x, ][intercluster_distances[x, ] > 0])))), min_dist)
            cluster1_cells <- rownames(new_clusters)[new_clusters[,1] == furthest_neighbors$cluster1[1]]
            cluster2_cells <- rownames(new_clusters)[new_clusters[,1] == furthest_neighbors$cluster2[1]]
            # Check if this pair of clusters should remain split
            comparison <- .runPermutationTest(cluster1_name = "Cluster1",
                                              cluster1_cells = cluster1_cells,
                                              cluster1_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                           batches[cluster1_cells],
                                                                           NULL),
                                              cluster2_name  = "Cluster2",
                                              cluster2_cells = cluster2_cells,
                                              cluster2_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                           batches[cluster2_cells],
                                                                           NULL),
                                              alpha = alpha,
                                              n_iterations = n_iterations,
                                              n_trees = n_trees,
                                              use_variance = use_variance,
                                              min_accuracy = min_accuracy,
                                              min_connections = min_connections,
                                              max_repeat_errors = max_repeat_errors,
                                              collect_all_metrics = FALSE,
                                              sample_max = sample_max,
                                              downsampling_rate = downsampling_rate,
                                              min_reads = min_reads,
                                              input_matrix = input_matrix[current_cells,],
                                              nn_matrix = nn_matrix[current_cells, current_cells],
                                              comparison_records = comparison_records,
                                              feature_importance_records = NULL,
                                              P0_distance = NA,
                                              P_i_distance = NA,
                                              n_cores = n_cores,
                                              random_seed = random_seed)
            # Add to records
            tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                           tree_name = current_tree_name,
                                                           num_cells = n_cells,
                                                           resolution = res,
                                                           num_clusters = new_n_clust,
                                                           silhouette = NA,
                                                           neighbors_distance = furthest_neighbors$min_dist[1],
                                                           neighbors_mean_accuracy = comparison[["comparison_records"]][1,"mean_accuracy"],
                                                           neighbors_var_accuracy = comparison[["comparison_records"]][1,"var_accuracy"],
                                                           neighbors_percentile_accuracy = comparison[["comparison_records"]][1,"percentile_accuracy"],
                                                           neighbors_percentile_variance = comparison[["comparison_records"]][1,"percentile_variance"],
                                                           neighbors_decision = comparison[["result"]],
                                                           stop_branching_reason = NA))
            if (comparison[["result"]] == "merge") {
              stop <- TRUE
              stop_reason <- "Sufficiently overclustered based on furthest neighbors"
              tree_records[nrow(tree_records), "stop_branching_reason"] <- stop_reason
              # Resolve cells
              unresolved_cells <- unresolved_cells[!(unresolved_cells %in% current_cells)]
              # Progress
              pb$tick(length(current_cells))
            } else {
              level_counter <- level_counter + 1

              # After 5 new levels without a stop, test all neighbors
              if (level_counter >=5 & length(current_cells) > min_cluster_depth) {
                stop <- TRUE
                stop_reason <- "Will subcluster unresolved cells."
                tree_records[nrow(tree_records), "stop_branching_reason"] <- stop_reason
                # Resolve clusters with < 4 cells
                resolve_clusters <- names(table(new_clusters[,1]))[table(new_clusters[,1]) < 4]
                if (length(resolve_clusters) > 0) {
                  remove_cells <- rownames(multi_level_clusters[unresolved_cells, ])[multi_level_clusters[unresolved_cells,
                                                                                                          ncol(multi_level_clusters)] %in% resolve_clusters]
                  unresolved_cells <- unresolved_cells[!(unresolved_cells %in% remove_cells)]
                  # Progress
                  pb$tick(length(remove_cells))
                }
                # Test neighbor pairs of each remaining cluster
                neighbor_pairs <- do.call(rbind, lapply(
                  seq(1:nrow(intercluster_distances)),
                  FUN = function(x)
                    data.frame(cluster1 = as.numeric(sub("c", "", rownames(intercluster_distances)[x])),
                               cluster2 = as.numeric(sub("c", "", colnames(intercluster_distances)[intercluster_distances[x, ] ==
                                                                                                     min(intercluster_distances[x, ][intercluster_distances[x, ] > 0])])),
                               min_dist = min(intercluster_distances[x, ][intercluster_distances[x, ] > 0])))) %>%
                  dplyr::filter(!(cluster1 == furthest_neighbors$cluster1[1] & cluster2 == furthest_neighbors$cluster2[1]),
                                !(cluster1 == furthest_neighbors$cluster2[1] & cluster2 == furthest_neighbors$cluster1[1])) %>%
                  dplyr::mutate(comp = ifelse(cluster1 < cluster2, paste0(cluster1, "_", cluster2), paste0(cluster2, "_", cluster1))) %>%
                  dplyr::group_by(comp) %>% dplyr::slice(1)
                while (nrow(neighbor_pairs) > 0) {
                  # Test first neighbor pair
                  cluster1_name <- neighbor_pairs$cluster1[1]
                  cluster2_name <- neighbor_pairs$cluster2[1]
                  cluster1_cells <- rownames(new_clusters)[new_clusters[,1] == cluster1_name]
                  cluster2_cells <- rownames(new_clusters)[new_clusters[,1] == cluster2_name]
                  # Check if this pair of clusters should remain split
                  comparison <- .runPermutationTest(cluster1_name = "Cluster1",
                                                    cluster1_cells = cluster1_cells,
                                                    cluster1_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                 batches[cluster1_cells],
                                                                                 NULL),
                                                    cluster2_name  = "Cluster2",
                                                    cluster2_cells = cluster2_cells,
                                                    cluster2_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                 batches[cluster2_cells],
                                                                                 NULL),
                                                    alpha = alpha,
                                                    n_iterations = n_iterations,
                                                    n_trees = n_trees,
                                                    use_variance = use_variance,
                                                    min_accuracy = min_accuracy,
                                                    min_connections = min_connections,
                                                    max_repeat_errors = max_repeat_errors,
                                                    collect_all_metrics = FALSE,
                                                    sample_max = sample_max,
                                                    downsampling_rate = downsampling_rate,
                                                    min_reads = min_reads,
                                                    input_matrix = input_matrix[current_cells, ],
                                                    nn_matrix = nn_matrix[current_cells, current_cells],
                                                    comparison_records = comparison_records,
                                                    feature_importance_records = NULL,
                                                    P0_distance = NA,
                                                    P_i_distance = NA,
                                                    n_cores = n_cores,
                                                    random_seed = random_seed)
                  # Add to records
                  tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                                 tree_name = current_tree_name,
                                                                 num_cells = n_cells,
                                                                 resolution = res,
                                                                 num_clusters = new_n_clust,
                                                                 silhouette = NA,
                                                                 neighbors_distance = neighbor_pairs$min_dist[1],
                                                                 neighbors_mean_accuracy = comparison[["comparison_records"]][1,"mean_accuracy"],
                                                                 neighbors_var_accuracy = comparison[["comparison_records"]][1,"var_accuracy"],
                                                                 neighbors_percentile_accuracy = comparison[["comparison_records"]][1,"percentile_accuracy"],
                                                                 neighbors_percentile_variance = comparison[["comparison_records"]][1,"percentile_variance"],
                                                                 neighbors_decision = comparison[["result"]],
                                                                 stop_branching_reason = NA))
                  # Remove row
                  neighbor_pairs <- neighbor_pairs[-1,]
                  # If merge
                  if (comparison[["result"]] == "merge") {
                    tree_records[nrow(tree_records), "stop_branching_reason"] <- "Current neighbor pair is sufficiently overclustered"
                    # Resolve cells
                    remove_cells <- c(cluster1_cells, cluster2_cells)[c(cluster1_cells, cluster2_cells) %in% unresolved_cells]

                    unresolved_cells <- unresolved_cells[!(unresolved_cells %in% remove_cells)]
                    # Remove neighbor pair of cluster2 if present
                    neighbor_pairs <- neighbor_pairs %>% dplyr::filter(cluster1 != cluster2_name)
                    # Progress
                    pb$tick(length(remove_cells))
                  }
                }
              }
            }
          } else {
            # If all but 1 cluster is a singleton, stop
            stop <- TRUE
            stop_reason <- "All but 1 cluster are singletons"
            # Add to records
            tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                           tree_name = current_tree_name,
                                                           num_cells = n_cells,
                                                           resolution = res,
                                                           num_clusters = new_n_clust,
                                                           silhouette = NA,
                                                           neighbors_distance = NA,
                                                           neighbors_mean_accuracy = NA,
                                                           neighbors_var_accuracy = NA,
                                                           neighbors_percentile_accuracy = NA,
                                                           neighbors_percentile_variance = NA,
                                                           neighbors_decision = NA,
                                                           stop_branching_reason = stop_reason))
            # Resolve clusters with < 4 cells
            resolve_clusters <- names(table(new_clusters[,1]))[table(new_clusters[,1]) < 4]
            remove_cells <- rownames(multi_level_clusters[unresolved_cells, ])[multi_level_clusters[unresolved_cells,
                                                                                                    ncol(multi_level_clusters)] %in% resolve_clusters]
            unresolved_cells <- unresolved_cells[!(unresolved_cells %in% remove_cells)]
            # Progress
            pb$tick(length(remove_cells))
          }
        }
        res <- res + gap
        gap_counter <- gap_counter + 1
        if (gap_counter >= 5) {
          gap <- gap*2
          gap_counter <- 0
        }
        # Round, otherwise R has addition issues
        res <- round(res, decimal_places)
      }
    }
  }
  # Clean up
  rm(snn_matrix)
  rm(nn_matrix)
  rm(dist_matrix)
  rm(reduction)
  rm(input_matrix)

  ## We now have our multi-level clustering results across a range of resolutions
  # If necessary, reconcile into strictly hierarchical clustering tree using MRtree
  if ((ncol(multi_level_clusters) == 2 & dplyr::n_distinct(multi_level_clusters[,1]) == 1) | ncol(multi_level_clusters) == 1) {
    cluster_tree <- multi_level_clusters
  } else {
    multi_level_clusters <- as.matrix(multi_level_clusters)
    # Convert to numerical for MRtree
    multi_level_clusters <- apply(multi_level_clusters, 2, function(x) as.numeric(as.factor(x)))
    # Run MRtree
    suppressMessages(invisible(utils::capture.output(mrtree_output <- mrtree::mrtree(multi_level_clusters, n.cores = n_cores, consensus = FALSE,
                                                                                     augment.path = FALSE, verbose = FALSE))))
    cluster_tree <- data.frame(mrtree_output$labelmat.mrtree)
    # Add to records
    tree_records <- rbind(tree_records, data.frame(tree_type = "subtree",
                                                   tree_name = current_tree_name,
                                                   num_cells = n_cells,
                                                   resolution = NA,
                                                   num_clusters = dplyr::n_distinct(cluster_tree[, ncol(cluster_tree)]),
                                                   silhouette = NA,
                                                   neighbors_distance = NA,
                                                   neighbors_mean_accuracy = NA,
                                                   neighbors_var_accuracy = NA,
                                                   neighbors_percentile_accuracy = NA,
                                                   neighbors_percentile_variance = NA,
                                                   neighbors_decision = NA,
                                                   stop_branching_reason = "Subtree after MRtree."))
  }
  # Remove level where all cells belong to the same cluster (if present)
  # Unless it's the only column
  if (dplyr::n_distinct(cluster_tree[,1]) == 1 & ncol(cluster_tree) > 1) {
    cluster_tree <- cluster_tree[, -1, drop = FALSE]
  }
  # Add cell names as row names
  rownames(cluster_tree) <- cell_IDs

  return(list("cluster_tree" = cluster_tree,
              "tree_records" = tree_records))
}

#' Infer clustering tree ---------------------------
#'
#' Generate clustering tree from provided, pre-generated clusters. Provide a set
#' of cluster labels and either a dimensionality reduction or distance matrix.
#' If a dimensionality reduction is provided, centroid distances will be
#' calculated and used.
#'
#' @param cluster_labels  A named vector of cluster IDs. Names must correspond
#' to cell IDs.
#' @param dist_matrix An optional distance matrix of cell to cell distances
#' (based on dimensionality reduction cell embeddings).
#' @param reduction An optional matrix of dimensionality reduction cell
#' embeddings to be used for distance calculations.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Can be set to \code{FALSE} for a
#' cleaner output.
#'
#' @return A clustering tree as a dataframe.
#' @export
#'
inferTree <- function(cluster_labels,
                      dist_matrix = NULL,
                      reduction = NULL,
                      verbose = TRUE) {
  .validInput(cluster_labels, name = "cluster_labels")
  .validInput(dist_matrix, name = "dist_matrix", names(cluster_labels))
  .validInput(reduction, name = "reduction", list("inferTree", names(cluster_labels)))
  .validInput(verbose, name = "verbose")

  if (is.null(reduction) & is.null(dist_matrix)) {
    stop("Please provide input to either 'reduction' or 'dist_matrix'.")
  } else if (!is.null(reduction) & !is.null(dist_matrix)) {
    warning("Input provided to both 'reduction' and 'dist_matrix'. Only input to 'dist_matrix' was used.")
    distance_description <- "from provided 'dist_matrix'"
    distance_approx <- FALSE
  } else if (is.null(reduction) & !is.null(dist_matrix)) {
    distance_description <- "from provided 'dist_matrix'"
    distance_approx <- FALSE
  } else if (!is.null(reduction) & is.null(dist_matrix)) {
    distance_description <- "calculated from provided 'reduction'"
    distance_approx <- TRUE
  }
  n_clusters <- dplyr::n_distinct(cluster_labels)

  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Inferring clustering tree from ", n_clusters, " provided clusters, using distances ", distance_description, "..")
  # Create initial clustering tree dataframe
  initial_tree <- data.frame(L1 = 1,
                             L2 = as.numeric(as.factor(cluster_labels)))
  # Create hierarchy
  cluster_tree <- .optimizeTree(cluster_tree = initial_tree,
                                reduction = reduction,
                                dist_matrix = dist_matrix,
                                distance_approx = distance_approx)
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Inferred clustering tree has ", ncol(cluster_tree), " levels.")
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Labeling clusters according to CHOIR conventions..")
  cluster_tree <- .checkClusterLabels(cluster_tree)

  # Add cell IDs as rownames to cluster tree
  rownames(cluster_tree) <- names(cluster_labels)

  return(cluster_tree)
}
