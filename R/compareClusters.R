#' Compare any two clusters using CHOIR's random forest classifier permutation
#' testing approach
#'
#' This function will take two provided clusters and assess whether they are
#' distinguishable by a permutation test using random forest classifier
#' prediction accuracies.
#'
#' @param object An object of class 'Seurat', 'SingleCellExperiment', or
#' 'ArchRProject'. Not used if values are provided for parameters
#' 'input_matrix' and 'nn_matrix'.
#' @param key The name under which CHOIR-related data is retrieved from the
#' object. Defaults to 'CHOIR'. Not used if values are provided for parameters
#' 'input_matrix' and 'nn_matrix'.
#' @param cluster1_cells A character vector of cell names belonging to cluster
#' 1.
#' @param cluster2_cells A character vector of cell names belonging to cluster
#' 2.
#' @param ident1 A string indicating the label for cluster 1.
#' @param ident2 A string indicating the label for cluster 2.
#' @param group_by A string indicating the column of cluster labels that
#' 'ident1' and 'ident2' belong to.
#' @param alpha A numeric value indicating the significance level used for
#' permutation test comparisons of cluster prediction accuracies. Defaults to
#' 0.05.
#' @param feature_set A string indicating whether to train random forest
#' classifiers on 'all' features or only variable ('var') features. Defaults to
#' 'var'.
#' @param exclude_features A character vector indicating features that should be
#' excluded from input to the random forest classifier. Default = \code{NULL}
#' will not exclude any features.
#' @param n_iterations A numeric value indicating the number of iterations run
#' for each permutation test comparison. Defaults to 100.
#' @param n_trees A numeric value indicating the number of trees in each random
#' forest. Defaults to 50.
#' @param use_variance A boolean value indicating whether to use the variance of
#' the random forest accuracy scores as part of the permutation test threshold.
#' Defaults to \code{TRUE}.
#' @param min_accuracy A numeric value indicating the minimum accuracy required
#' of the random forest classifier, below which clusters will be automatically
#' merged. Defaults to 0.5 (chance).
#' @param min_connections A numeric value indicating the minimum number of
#' nearest neighbors between two clusters for them to be considered 'adjacent'.
#' Non-adjacent clusters will not be merged. Defaults to 1.
#' @param max_repeat_errors Used to account for situations in which random
#' forest classifier errors are concentrated among a few cells that are
#' repeatedly misassigned. A numeric value indicating the maximum number of such
#' 'repeat errors' that will be taken into account. If set to 0, 'repeat errors'
#' will not be evaluated. Defaults to 20.
#' @param collect_all_metrics A boolean value indicating whether to collect and
#' save additional metrics from the random forest classifier comparisons,
#' including feature importances and tree depth. Defaults to \code{FALSE}.
#' @param sample_max A numeric value indicating the maximum number of cells used
#' per cluster to train/test each random forest classifier. Default = \code{Inf}
#' does not cap the number of cells used.
#' @param downsampling_rate A numeric value indicating the proportion of cells used
#' per cluster to train/test each random forest classifier. Default = "auto" sets
#' the downsampling rate according to the dataset size, for efficiency.
#' @param min_reads A numeric value used to filter out features prior to input
#' to the random forest classifier. Default = \code{NULL} will filter out
#' features with 0 counts for the current clusters being compared. Numeric input
#' values should be used only with count input matrices, e.g., ATAC tile
#' matrices, whereby at least 1 read will be required for the provided number
#' of cells.
#' @param normalization_method A character string or vector indicating which
#' normalization method to use. In general, input data should be supplied to
#' CHOIR after normalization, except in cases when the user wishes to use
#' \code{Seurat::SCTransform()} normalization. Permitted values are 'none' or
#' 'SCTransform'. Defaults to 'none'.
#' @param batch_labels If applying batch correction, a character string or
#' vector indicating the name of the column containing the batch labels.
#' Defaults to \code{NULL}.
#' @param use_assay For Seurat or SingleCellExperiment objects, a character
#' string or vector indicating the assay(s) to use in the provided object.
#' Default = \code{NULL} will choose the current active assay for Seurat objects
#' and the \code{logcounts} assay for SingleCellExperiment objects.
#' @param use_slot For Seurat objects, a character string or vector indicating
#' the layers(s) — previously known as slot(s) — to use in the provided object.
#' Default = \code{NULL} will choose a layer/slot based on the selected assay.
#' If a non-standard assay is provided, do not leave \code{use_slot} as
#' \code{NULL}.
#' @param ArchR_matrix For ArchR objects, a character string or vector
#' indicating which matri(ces) to use in the provided object. Default =
#' \code{NULL} will use the 'TileMatrix' for ATAC-seq data or the
#' 'GeneExpressionMatrix' for RNA-seq data.
#' @param atac A boolean value or vector indicating whether the provided data is
#' ATAC-seq data. Defaults to \code{FALSE}. For multi-omic datasets containing
#' ATAC-seq data, it is important to supply this parameter as a vector
#' corresponding to each modality in order.
#' @param input_matrix An optional matrix containing the feature x cell data on
#' which to train the random forest classifiers. Default = \code{NULL} will use
#' the feature x cell matri(ces) indicated by function \code{buildTree()}.
#' @param nn_matrix An optional matrix containing the nearest neighbor adjacency
#' of the cells. Default = \code{NULL} will look for the adjacency matri(ces)
#' generated by function \code{buildTree()}.
#' @param var_features An optional character vector of variable features to be
#' used for subsequent clustering steps. Default = \code{NULL} will use
#' the variable features identified by function \code{buildTree()}.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param random_seed A numeric value indicating the random seed to be used.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Can be set to \code{FALSE} for a
#' cleaner output.
#'
#' @return Returns a list containing the following elements: \describe{
#'   \item{comparison_result}{A string, either "merge" or "split", indicating
#'   the result of the comparison.}
#'   \item{comparison_records}{A dataframe including the metrics recorded for
#'   the comparison}
#'   \item{feature_importances}{If 'collect_all_metrics' is true, a dataframe
#'   containing the feature importance scores for each gene in the comparison}
#'   }
#'
#' @export
#'
compareClusters <- function(object = NULL,
                            key = "CHOIR",
                            cluster1_cells = NULL,
                            cluster2_cells = NULL,
                            ident1 = NULL,
                            ident2 = NULL,
                            group_by = NULL,
                            alpha = 0.05,
                            feature_set = "var",
                            exclude_features = NULL,
                            n_iterations = 100,
                            n_trees = 50,
                            use_variance = TRUE,
                            min_accuracy = 0.5,
                            min_connections = 1,
                            max_repeat_errors = 20,
                            collect_all_metrics = FALSE,
                            sample_max = Inf,
                            downsampling_rate = "auto",
                            min_reads = NULL,
                            normalization_method = "none",
                            batch_labels = NULL,
                            use_assay = NULL,
                            use_slot = NULL,
                            ArchR_matrix = NULL,
                            atac = FALSE,
                            input_matrix = NULL,
                            nn_matrix = NULL,
                            var_features = NULL,
                            n_cores = NULL,
                            random_seed = 1,
                            verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "compareClusters")
  .validInput(key, "key", list("compareClusters", object))
  .validInput(cluster1_cells, "cluster1_cells", cluster2_cells)
  .validInput(cluster2_cells, "cluster2_cells", cluster1_cells)
  .validInput(group_by, "group_by", object)
  .validInput(ident1, "ident1", list(object, group_by, ident2))
  .validInput(ident2, "ident2", list(object, group_by, ident1))
  .validInput(alpha, "alpha")
  .validInput(feature_set, "feature_set", list("compareClusters"))
  .validInput(var_features, "var_features", input_matrix)
  .validInput(exclude_features, "exclude_features")
  .validInput(n_iterations, "n_iterations")
  .validInput(n_trees, "n_trees")
  .validInput(use_variance, "use_variance")
  .validInput(min_accuracy, "min_accuracy")
  .validInput(min_connections, "min_connections")
  .validInput(max_repeat_errors, "max_repeat_errors")
  .validInput(sample_max, "sample_max")
  .validInput(downsampling_rate, "downsampling_rate")
  .validInput(min_reads, "min_reads")
  .validInput(batch_labels, "batch_labels", object)
  .validInput(collect_all_metrics, "collect_all_metrics")
  .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
  .validInput(use_slot, "use_slot", list(object, use_assay, FALSE, NULL))
  .validInput(ArchR_matrix, "ArchR_matrix", list(object, FALSE, NULL))
  # Number of modalities & object type
  if (methods::is(object, "ArchRProject")) {
    n_modalities <- max(length(ArchR_matrix), 1)
    object_type <- "ArchRProject"
    .requirePackage("ArchR", installInfo = "Instructions at archrproject.com")
  } else {
    n_modalities <- max(length(use_assay), 1)
    if (methods::is(object, "Seurat")) {
      object_type <- "Seurat"

      if (is.null(use_assay)) {
        if ("Assay5" %in% methods::is(object[[Seurat::DefaultAssay(object)]])) {
          seurat_version <- "v5"
        } else if ("Assay" %in% methods::is(object[[Seurat::DefaultAssay(object)]])) {
          seurat_version <- "v4"
        } else {
          stop("Assay '", Seurat::DefaultAssay(object), "' provided for parameter 'use_assay' is not of class 'Assay' or 'Assay5', please supply valid input!")
        }
      } else {
        if ("Assay5" %in% methods::is(object[[use_assay[1]]])) {
          seurat_version <- "v5"
        } else if ("Assay" %in% methods::is(object[[use_assay[1]]])) {
          seurat_version <- "v4"
        } else {
          stop("Assay '", use_assay[1], "' provided for parameter 'use_assay' is not of class 'Assay' or 'Assay5', please supply valid input!")
        }
      }
    } else if (methods::is(object, "SingleCellExperiment")) {
      object_type <- "SingleCellExperiment"
      .requirePackage("SingleCellExperiment", source = "bioc")
    }
  }
  .validInput(normalization_method, "normalization_method", list(object, n_modalities, use_assay))
  .validInput(atac, "atac", n_modalities)
  .validInput(input_matrix, "input_matrix", c(cluster1_cells, cluster2_cells))
  .validInput(nn_matrix, "nn_matrix", c(cluster1_cells, cluster2_cells))
  .validInput(n_cores, "n_cores")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # Get cell IDs
  cell_IDs <- .getCellIDs(object, use_assay = use_assay)

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }
  # Set downsampling rate
  if (downsampling_rate == "auto") {
    downsampling_rate <- min(1, (1/2)^(log10(length(cell_IDs)/5000)))
  }
  # Random seed reproducibility
  if (n_cores > 1) {
    RNGkind("L'Ecuyer-CMRG")
  }

  # ---------------------------------------------------------------------------
  # Prepare object for random forest comparison
  # ---------------------------------------------------------------------------
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 1/2) Preparing object..")

  # Extract batch labels
  if (!is.null(batch_labels)) {
    batches <- .retrieveData(object = object, key = key, type = "cell_metadata", name = batch_labels)
    names(batches) <- cell_IDs
  }

  # Cluster names / cell IDs
  if (!is.null(cluster1_cells)) {
    if (!is.null(ident1)) {
      warning("Input values for 'ident1', 'ident2', and 'group_by' are not used when values for 'cluster1_cells' and 'cluster2_cells' are provided.")
    }
    ident1 <- "Cluster1"
    ident2 <- "Cluster2"
  } else if (!is.null(ident1)) {
    if (methods::is(object, "Seurat")) {
      cluster1_cells <- rownames(object@meta.data)[object@meta.data[, group_by] == ident1]
      cluster2_cells <- rownames(object@meta.data)[object@meta.data[, group_by] == ident2]
    } else if (methods::is(object, "SingleCellExperiment")) {
      cluster1_cells <- rownames(object@colData)[object@colData[, group_by] == ident1]
      cluster2_cells <- rownames(object@colData)[object@colData[, group_by] == ident2]
    } else if (methods::is(object, "ArchRProject")) {
      cluster1_cells <- rownames(object@cellColData)[object@cellColData[, group_by] == ident1]
      cluster2_cells <- rownames(object@cellColData)[object@cellColData[, group_by] == ident2]
    }
  } else {
    stop("Must provide input values for either 'cluster1_cells' & 'cluster2_cells' OR 'ident1', 'ident2', & 'group_by'.")
  }

  # Extract input matrix/matrices
  if (!is.null(input_matrix)) {
    input_matrix_provided <- TRUE
    input_matrix <- .getMatrix(use_matrix = input_matrix,
                               use_features = var_features,
                               exclude_features = exclude_features,
                               use_cells = c(cluster1_cells, cluster2_cells),
                               verbose = FALSE)
    input_matrix <- BiocGenerics::t(input_matrix)
    # Feature names
    features <- colnames(input_matrix)
  } else {
    input_matrix_provided <- FALSE
    # Determine which input & nn matrices to use
    roots <- unlist(stringr::str_extract_all(c(ident1, ident2), "P\\d*"))
    if (dplyr::n_distinct(roots) == 1) {
      if (n_input_matrices > 1) {
        use_input_matrix <- roots[1]
      } else {
        use_input_matrix <- "P0"
      }
      use_nn_matrix <- roots[1]
    } else {
      use_input_matrix <- "P0"
      use_nn_matrix <- "P0"
    }
    # Use all / var features
    if (feature_set == "all") {
      use_features <- NULL
    } else {
      use_features <- .retrieveData(object, key, "var_features", paste0(use_input_matrix, "_var_features"))
      if (is.null(use_features)) {
        stop("Could not find variable features under key '", key,
             "'. Please provide as input to parameter 'var_features'.")
      }
    }
    # Extract input matrix for random forest comparisons
    if (n_modalities > 1) {
      input_matrix_list <- vector("list", n_modalities)
      for (m in 1:n_modalities) {
        # Match input arguments
        use_assay_m <- .matchArg(use_assay, m)
        use_slot_m <- .matchArg(use_slot, m)
        ArchR_matrix_m <- .matchArg(ArchR_matrix, m)
        # Features
        if (length(use_features_subtree) > 1) {
          use_features_m <- use_features[[m]]
        } else {
          use_features_m <- use_features
        }
        input_matrix_list[[m]] <- .getMatrix(object = object,
                                             use_assay = use_assay_m,
                                             use_slot = use_slot_m,
                                             ArchR_matrix = ArchR_matrix_m,
                                             use_features = use_features_m,
                                             exclude_features = exclude_features,
                                             use_cells = c(cluster1_cells, cluster2_cells),
                                             verbose = FALSE)
        if (is.null(input_matrix_list[[m]])) {
          stop("Could not find input matrix under key '", key,
               "'. Please provide as input to parameter 'input_matrix'.")
        }
        if (normalization_method_m == "SCTransform") {
          input_matrix_list[[m]] <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(input_matrix_list[[m]]),
                                                                         return.only.var.genes = FALSE,
                                                                         seed.use = random_seed,
                                                                         verbose = FALSE)@assays$SCT@scale.data)
        }
      }
      input_matrix <- do.call(rbind, input_matrix_list)
      input_matrix <- BiocGenerics::t(input_matrix)
    } else {
      input_matrix <- .getMatrix(object = object,
                                 use_assay = use_assay,
                                 use_slot = use_slot,
                                 ArchR_matrix = ArchR_matrix,
                                 use_features = use_features,
                                 exclude_features = exclude_features,
                                 use_cells = c(cluster1_cells, cluster2_cells),
                                 verbose = FALSE)
      if (is.null(input_matrix)) {
        stop("Could not find input matrix under key '", key,
             "'. Please provide as input to parameter 'input_matrix'.")
      }
      if (normalization_method == "SCTransform") {
        input_matrix <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(input_matrix),
                                                             return.only.var.genes = FALSE,
                                                             seed.use = random_seed,
                                                             verbose = FALSE)@assays$SCT@scale.data)
      }
      input_matrix <- BiocGenerics::t(input_matrix)
    }
    # Add features to vector
    features <- colnames(input_matrix)
  }

  # Extract nearest neighbor matrix/matrices
  if (!is.null(nn_matrix)) {
    nn_matrix_provided <- TRUE
  } else if (min_connections == 0 & collect_all_metrics == FALSE) {
    nn_matrix <- NULL
  } else {
    nn_matrix_provided <- FALSE
    nn_matrix <- .retrieveData(object, key, "graph", paste0(use_input_matrix, "_graph_nn"))
    if (is.null(nn_matrix)) {
      stop("Could not find nearest neighbor adjacency matrix under key '", key,
           "'. Please provide as input to parameter 'nn_matrix'.")
    }
  }

  # Create data frame to record comparison
  # Data frame to track comparisons to avoid redundancy
  all_metrics <- c('comparison', 'cluster1_size', 'cluster2_size', 'sample_size',
                   'mean_accuracy', 'var_accuracy', 'mean_errors',
                   'mean_permuted_accuracy', 'var_permuted_accuracy',
                   'percentile_accuracy', 'percentile_variance',
                   'n_repeat_errors1', 'n_repeat_errors2',
                   'mean_repeat_errors1', 'mean_repeat_errors2',
                   'mean_modified_accuracy', 'var_modified_accuracy',
                   'percentile_modified_accuracy', 'percentile_modified_variance',
                   'batches_used', 'batch_mean_accuracies',
                   'connectivity', 'time',
                   'decision')
  selected_metrics <- all_metrics[c(1:11,
                                    `if`(collect_all_metrics == TRUE | max_repeat_errors > 0, 12:15, NULL),
                                    `if`(max_repeat_errors > 0, 16:19, NULL),
                                    `if`(!is.null(batch_labels), 20:21, NULL),
                                    `if`(collect_all_metrics == TRUE | min_connections > 0, 22, NULL),
                                    23:24)]
  comparison_records <- data.frame(matrix(ncol = length(selected_metrics), nrow = 0))
  colnames(comparison_records) <- selected_metrics

  # Initialize feature importance records
  if (collect_all_metrics == TRUE) {
    feature_importance_records <- data.frame(matrix(ncol = (ncol(input_matrix)+2), nrow = 0))
    colnames(feature_importance_records) <- c('cluster1', 'cluster2', colnames(input_matrix))
  } else {
    feature_importance_records <- NULL
  }

  # ---------------------------------------------------------------------------
  # Run comparison
  # ---------------------------------------------------------------------------
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/2) Comparing clusters..")

  # Run comparison
  comparison_output <- .runPermutationTest(cluster1_name = ident1,
                                           cluster1_cells = cluster1_cells,
                                           cluster1_cell_batches = `if`(!is.null(batch_labels),
                                                                        batches[cluster1_cells],
                                                                        NULL),
                                           cluster2_name = ident2,
                                           cluster2_cells = cluster2_cells,
                                           cluster2_cell_batches = `if`(!is.null(batch_labels),
                                                                        batches[cluster2_cells],
                                                                        NULL),
                                           alpha = alpha,
                                           n_iterations = n_iterations,
                                           n_trees = n_trees,
                                           use_variance = use_variance,
                                           min_accuracy = min_accuracy,
                                           min_connections = min_connections,
                                           max_repeat_errors = max_repeat_errors,
                                           collect_all_metrics = collect_all_metrics,
                                           sample_max = sample_max,
                                           downsampling_rate = downsampling_rate,
                                           min_reads = min_reads,
                                           input_matrix = input_matrix,
                                           nn_matrix = nn_matrix,
                                           comparison_records = comparison_records,
                                           feature_importance_records = feature_importance_records,
                                           P0_distance = NA,
                                           P_i_distance = NA,
                                           comparison_start_time = Sys.time(),
                                           n_cores = n_cores,
                                           random_seed = random_seed)

  if (verbose) {
    if (comparison_output[["result"]] == "merge") {
      message(format(Sys.time(), "%Y-%m-%d %X"), " : Using current thresholds, clusters should be merged.")
    } else {
      message(format(Sys.time(), "%Y-%m-%d %X"), " : Using current thresholds, clusters should remain split.")
    }
  }

  if (max_repeat_errors > 0 & comparison_output[["comparison_records"]]$decision[1] == "split: repeat error") {
    message(" - Repeatedly misassigned cells affected the result of this comparison. Set 'max_repeat_errors' parameter to 0 if you would like to disable this setting.")
  }
  if (min_connections > 0 & comparison_output[["comparison_records"]]$decision[1] == "split: min connections") {
    message(" - Clusters were split due to the minimum number of nearest neighbor connections.")
  }
  if (comparison_output[["comparison_records"]]$decision[1] == "merge: min accuracy") {
    message(" - Clusters were merged due to the minimum accuracy threshold.")
  }
  if (comparison_output[["comparison_records"]]$decision[1] == "merge: batch-dependent") {
    message(" - Clusters were merged due to batch-dependence.")
  }

  # Output as named list
  if (collect_all_metrics == TRUE) {
    output_list <- list("comparison_result" = comparison_output[["result"]],
                        "comparison_records" = comparison_output[["comparison_records"]],
                        "feature_importances" = comparison_output[["feature_importances"]])
  } else {
    output_list <- list("comparison_result" = comparison_output[["result"]],
                        "comparison_records" = comparison_output[["comparison_records"]])
  }
  return(output_list)
}
