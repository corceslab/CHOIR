#' Compare any two clusters using CHOIR's random forest classifier permutation
#' testing approach
#'
#' This function will take two provided clusters and assess whether they are
#' distinguishable by a permutation test using random forest classifier
#' prediction accuracies.
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{ArchRProject}. For multi-omic data, we recommend using
#' \code{ArchRProject} objects. Not used if values are provided for parameters
#' \code{input_matrix} and \code{nn_matrix}.
#' @param key The name under which CHOIR-related data for this run is stored in
#' the object. Defaults to “CHOIR”. Not used if input is provided for
#' parameters \code{input_matrix} and \code{nn_matrix}.
#' @param cluster1_cells A character vector of cell names belonging to cluster
#' 1.
#' @param cluster2_cells A character vector of cell names belonging to cluster
#' 2.
#' @param ident1 A string indicating the label for cluster 1.
#' @param ident2 A string indicating the label for cluster 2.
#' @param group_by A string indicating the column of cluster labels that
#' 'ident1' and 'ident2' belong to.
#' @param alpha A numerical value indicating the significance level used for
#' permutation test comparisons of cluster distinguishability. Defaults to 0.05.
#' @param feature_set A string indicating whether to train random forest
#' classifiers on “all” features or only variable (“var”) features. Defaults to
#' “var”. Computational time and memory required may increase if more features
#' are used. Using all features instead of variable features may result in more
#' conservative cluster calls.
#' @param exclude_features A character vector indicating features that should be
#' excluded from input to the random forest classifier. Defaults to \code{NULL},
#' which means that no features will be excluded. This parameter can be used,
#' for example, to exclude features correlated with cell quality, such as
#' mitochondrial genes. Failure to exclude problematic features could result in
#' clusters driven by cell quality, while over-exclusion of features could
#' reduce the ability of CHOIR to distinguish cell populations that differ by
#' those features.
#' @param n_iterations A numerical value indicating the number of iterations run
#' for each permutation test comparison. Increasing the number of iterations
#' will approximately linearly increase the computational time required but
#' provide a more accurate estimation of the significance of the permutation
#' test. Decreasing the number of iterations runs the risk of leading to
#' underclustering due to lack of statistical power. The default value, 100
#' iterations, was selected because it avoids underclustering, while minimizing
#' computational time and the diminishing returns from running CHOIR with
#' additional iterations.
#' @param n_trees A numerical value indicating the number of trees in each
#' random forest. Defaults to 50. Increasing the number of trees is likely to
#' increase the computational time required. Though not entirely predictable,
#' increasing the number of trees up to a point may enable more nuanced
#' distinctions, but is likely to provide diminishing returns.
#' @param use_variance A Boolean value indicating whether to use the variance of
#' the random forest accuracy scores as part of the permutation test threshold.
#' Defaults to \code{TRUE}. Setting this parameter to \code{FALSE} will make
#' CHOIR considerably less conservative, identifying more clusters, particularly
#' on large datasets.
#' @param min_accuracy A numerical value indicating the minimum accuracy
#' required of the random forest classifier, below which clusters will be
#' automatically merged. Defaults to 0.5, representing the random chance
#' probability of assigning correct cluster labels; therefore, decreasing the
#' minimum accuracy is not recommended. Increasing the minimum accuracy will
#' lead to more conservative cluster assignments and will often decrease the
#' computational time required, because fewer cluster comparisons may be needed.
#' @param min_connections A numerical value indicating the minimum number of
#' nearest neighbors between two clusters for those clusters to be considered
#' adjacent. Non-adjacent clusters will not be merged. Defaults to 1. This
#' threshold allows CHOIR to avoid running the full permutation test comparison
#' for clusters that are highly likely to be distinct, saving computational
#' time. The intent of this parameter is only to avoid running permutation test
#' comparisons between clusters that are so different that they should not be
#' merged. Therefore, we do not recommend increasing this parameter value
#' beyond 10, as higher values may result in instances of overclustering.
#' @param max_repeat_errors A numerical value indicating the maximum number of
#' repeatedly mislabeled cells that will be taken into account during the
#' permutation tests. This parameter is used to account for situations in which
#' random forest classifier errors are concentrated among a few cells that are
#' repeatedly misassigned. If set to 0, such repeat errors will not be
#' evaluated. Defaults to 20. These situations are relatively infrequent, but
#' setting this parameter to lower values (especially 0) may result in
#' underclustering due to a small number of intermediate cells. Setting this
#' parameter to higher values may lead to instances of overclustering and is not
#' recommended.
#' @param collect_all_metrics A Boolean value indicating whether to collect and
#' save additional metrics from the random forest classifiers, including feature
#' importances for every comparison. Defaults to \code{FALSE}. Setting this
#' parameter to \code{TRUE} will slightly increase the computational time
#' required. This parameter has no effect on the final cluster calls.
#' @param sample_max A numerical value indicating the maximum number of cells to
#' be sampled per cluster to train/test each random forest classifier. Defaults
#' to \code{Inf} (infinity), which does not cap the number of cells used, so all
#' cells will be used in all comparisons. Decreasing this parameter may decrease
#' the computational time required, but may result in instances of
#' underclustering. If input is provided to both the \code{downsampling_rate}
#' and \code{sample_max} parameters, the minimum resulting cell number is
#' calculated and used for each comparison.
#' @param downsampling_rate A numerical value indicating the proportion of cells
#' to be sampled per cluster to train/test each random forest classifier. For
#' efficiency, the default value, "auto", sets the downsampling rate according
#' to the dataset size. Decreasing this parameter may decrease the computational
#' time required, but may also make the final cluster calls more conservative.
#' If input is provided to both \code{downsampling_rate} and
#' \code{sample_max parameters}, the minimum resulting cell number is calculated
#' and used for each comparison.
#' @param min_reads A numeric value used to filter out features prior to input
#' to the random forest classifier. The default value, \code{NULL}, will filter
#' out features with 0 counts for the current clusters being compared. Higher
#' values should be used with caution, but may increase the signal-to-noise
#' ratio encountered by the random forest classifiers.
#' @param normalization_method A character string or vector indicating which
#' normalization method to use. In general, input data should be supplied to
#' CHOIR after normalization, except when the user wishes to use
#' \code{Seurat SCTransform} normalization. Permitted values are “none” or
#' “SCTransform”. Defaults to “none”. Because CHOIR has not been tested
#' thoroughly with \code{SCTransform} normalization, we do not recommend this
#' approach at this time. For multi-omic datasets, provide a vector with a value
#' corresponding to each provided value of \code{use_assay} or
#' \code{ArchR_matrix} in the same order.
#' @param batch_labels A character string that, if applying batch correction,
#' specifies the name of the column in the input object metadata containing the
#' batch labels. Defaults to \code{NULL}.
#' @param use_assay For \code{Seurat} or \code{SingleCellExperiment} objects, a
#' character string or vector indicating the assay(s) to use in the provided
#' object. The default value, \code{NULL}, will choose the current active assay
#' for \code{Seurat} objects and the \code{logcounts} assay for
#' \code{SingleCellExperiment} objects.
#' @param use_slot For \code{Seurat} objects, a character string or vector
#' indicating the layers(s)—previously known as slot(s)—to use in the provided
#' object. The default value, \code{NULL}, will choose a layer/slot based on the
#' selected assay. If an assay other than "RNA", "sketch”, "SCT”, or
#' "integrated" is provided, you must specify a value for \code{use_slot}. For
#' multi-omic datasets, provide a vector with a value corresponding to each
#' provided value of \code{use_assay} in the same order.
#' @param ArchR_matrix For \code{ArchR} objects, a character string or vector
#' indicating which matrix or matrices to use in the provided object. The
#' default value, \code{NULL}, will use the “GeneScoreMatrix” for ATAC-seq data
#' or the “GeneExpressionMatrix” for RNA-seq data. For multi-omic datasets,
#' provide a vector with a value corresponding to each modality.
#' @param atac A Boolean value or vector indicating whether the provided data is
#' ATAC-seq data. For multi-omic datasets, provide a vector with a value
#' corresponding to each provided value of \code{use_assay} or
#' \code{ArchR_matrix} in the same order. Defaults to \code{FALSE}.
#' @param input_matrix An optional matrix containing the feature x cell data
#' provided by the user, on which to train the random forest classifiers. By
#' default, this parameter is set to \code{NULL}, and CHOIR will look for the
#' feature x cell matri(ces) indicated by function \code{buildTree}.
#' @param nn_matrix An optional matrix containing the nearest neighbor adjacency
#' of the cells, provided by the user. By default, this parameter is set to
#' \code{NULL}, and CHOIR will look for the adjacency matri(ces) generated by
#' function \code{buildTree}.
#' @param var_features An optional character vector of names of variable
#' features to be used for subsequent clustering steps. By default, this
#' parameter is set to \code{NULL}, and variable features previously identified
#' by function \code{buildTree} will be used. Input to this parameter is
#' required when a dimensionality reduction is supplied to parameter
#' \code{reduction}. For multi-omic datasets, concatenate feature names for all
#' modalities.
#' @param n_cores A numerical value indicating the number of cores to use for
#' parallelization. By default, CHOIR will use the number of available cores
#' minus 2. CHOIR is parallelized at the computation of permutation test
#' iterations. Therefore, any number of cores up to the number of iterations
#' will theoretically decrease the computational time required. In practice,
#' 8–16 cores are recommended for datasets up to 500,000 cells.
#' @param random_seed A numerical value indicating the random seed to be used.
#' Defaults to 1. CHOIR uses randomization throughout the generation and pruning
#' of the clustering tree. Therefore, changing the random seed may yield slight
#' differences in the final cluster assignments.
#' @param verbose A Boolean value indicating whether to use verbose output
#' during the execution of CHOIR. Defaults to \code{TRUE}, but can be set to
#' \code{FALSE} for a cleaner output.
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
    batches <- as.character(batches)
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
                   'batches_used', 'batch_mean_accuracies', 'batch_mean_variances',
                   'connectivity', 'time',
                   'decision')
  selected_metrics <- all_metrics[c(1:11,
                                    `if`(collect_all_metrics == TRUE | max_repeat_errors > 0, 12:15, NULL),
                                    `if`(max_repeat_errors > 0, 16:19, NULL),
                                    `if`(!is.null(batch_labels), 20:22, NULL),
                                    `if`(collect_all_metrics == TRUE | min_connections > 0, 23, NULL),
                                    24:25)]
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
                        "feature_importances" = comparison_output[["feature_importance_records"]])
  } else {
    output_list <- list("comparison_result" = comparison_output[["result"]],
                        "comparison_records" = comparison_output[["comparison_records"]])
  }
  return(output_list)
}
