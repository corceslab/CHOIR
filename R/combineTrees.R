#' Compile separately pruned clustering trees and standardize thresholds
#'
#' This function is designed for the use of CHOIR on atlas-scale data.
#' It will take the output from buildParentTree, and the output from subsequent
#' applications of CHOIR to each resulting subtree, and combine these results
#' into a single compiled tree, which is then further pruned to standardize
#' thresholds across each subtree.
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{ArchRProject} that was output from function \code{buildParentTree}.
#' For multi-omic data, we recommend using \code{ArchRProject} objects.
#' @param subtree_list A list containing the CHOIR records from each subtree.
#' @param key The name under which CHOIR-related data for this run is stored in
#' the object. Defaults to “CHOIR”.
#' @param alpha A numerical value indicating the significance level used for
#' permutation test comparisons of cluster distinguishability. Defaults to 0.05.
#' Decreasing the alpha value will yield more conservative clusters (fewer
#' clusters) and will often decrease the computational time required, because
#' fewer cluster comparisons may be needed.
#' @param p_adjust A string indicating which multiple comparison adjustment
#' method to use. Permitted values are “bonferroni”, “fdr”, and “none”. Defaults
#' to “bonferroni”. Other correction methods may be less conservative,
#' identifying more clusters, as CHOIR applies filters that reduce the total
#' number of tests performed.
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
#' time. Therefore, setting this parameter to 0 will increase the number of
#' permutation test comparisons run and, thus, the computational time. The
#' intent of this parameter is only to avoid running permutation test
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
#' @param distance_awareness A numerical value representing the distance
#' threshold above which a cluster will not merge with another cluster and
#' significance testing will not be used. Specifically, this value is a
#' multiplier applied to the distance between a cluster and its closest
#' distinguishable neighbor based on random forest comparison. Defaults to 2,
#' which sets this threshold at a two-fold increase in distance over the closest
#' distinguishable neighbor. This threshold allows CHOIR to avoid running the
#' full permutation test comparison for clusters that are highly likely to be
#' distinct, saving computational time. To omit all distance calculations and
#' perform permutation testing on all comparisons, set this parameter to
#' \code{Inf}. Setting this parameter to \code{Inf} or increasing the input
#' value will increase the number of permutation test comparisons run and, thus,
#' the computational time. In rare cases, very small distant clusters may be
#' erroneously merged when distance thresholds are not used. The intent of this
#' parameter is only to avoid running permutation test comparisons between
#' clusters that are so different that they should not be merged. We do not
#' recommend decreasing this parameter value below 1.5, as lower values may
#' result in instances of overclustering.
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
#' @param batch_correction_method A character string indicating which batch
#' correction method to use. Permitted values are “Harmony” and “none”. Defaults
#' to “none”. Batch correction should only be used when the different batches
#' are not expected to also have unique cell types or cell states. Using batch
#' correction would ensure that clusters do not originate from a single batch,
#' thereby making the final cluster calls more conservative.
#' @param batch_labels A character string that, if applying batch correction,
#' specifies the name of the column in the input object metadata containing the
#' batch labels. Defaults to \code{NULL}.
#' @param max_n_batch A numeric value used if applying batch correction,
#' indicating the maximum number batches to use in each permutation test.
#' Batches are selected in order from largest to smallest. Defaults to
#' \code{Inf}, which will use all batches that pass cell number thresholds. In
#' datasets that contain many batches (>10), where there is no logical way to
#' group batches prior to applying CHOIR, setting \code{max_n_batch} to a
#' smaller value (we suggest between 1 and 5) may help avoid excessive
#' downsampling when running the permutation tests, and thereby help avoid
#' underclustering.
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
#' @param countsplit A Boolean value indicating whether or not to use count
#' split input data (see \code{countsplit} package), such that one matrix of
#' counts is used for clustering tree generation and a separate matrix is used
#' for all random forest classifier permutation testing. Defaults to
#' \code{FALSE}. Enabling count splitting is likely to result in more
#' conservative final cluster calls and is likely to perform best in datasets
#' with high read depths.
#' @param countsplit_suffix A character vector indicating the suffixes that
#' distinguish the two count split matrices to be used. Suffixes are appended
#' onto the input string/vector for parameter \code{use_slot} for \code{Seurat}
#' objects, \code{use_assay} for \code{SingleCellExperiment} objects, or
#' \code{ArchR_matrix} for \code{ArchR} objects. When count splitting is
#' enabled, the default value \code{NULL} uses suffixes "_1" and "_2".
#' @param input_matrix An optional matrix containing the feature x cell data
#' provided by the user, on which to train the random forest classifiers. By
#' default, this parameter is set to \code{NULL}, and CHOIR will look for the
#' feature x cell matri(ces) indicated by function \code{buildParentTree}.
#' @param nn_matrix An optional matrix containing the nearest neighbor adjacency
#' of the cells, provided by the user. By default, this parameter is set to
#' \code{NULL}, and CHOIR will look for the adjacency matri(ces) generated by
#' function \code{buildParentTree}.
#' @param reduction An optional matrix of dimensionality reduction cell
#' embeddings provided by the user for subsequent clustering steps. By default,
#' this parameter is set to \code{NULL}, and CHOIR will look for the
#' dimensionality reductions generated by function \code{buildParentTree()}.
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
#'
#'
#'@return Returns the object with the following added data stored under the
#' provided key: \describe{
#'   \item{clusters}{Final clusters, full hierarchical cluster tree, and
#'   stepwise cluster results for each progressive pruning step}
#'   \item{parameters}{Record of parameter values used}
#'   \item{records}{Metadata for decision points during hierarchical tree
#'   construction, all recorded permutation test comparisons, and feature
#'   importance scores from all comparisons}
#'   }
#'
#' @export
combineTrees <- function(object,
                         subtree_list,
                         key = "CHOIR",
                         alpha = NULL,
                         p_adjust = NULL,
                         feature_set = NULL,
                         exclude_features = NULL,
                         n_iterations = NULL,
                         n_trees = NULL,
                         use_variance = NULL,
                         min_accuracy = NULL,
                         min_connections = NULL,
                         max_repeat_errors = NULL,
                         distance_awareness = NULL,
                         collect_all_metrics = NULL,
                         sample_max = NULL,
                         downsampling_rate = NULL,
                         min_reads = NULL,
                         normalization_method = NULL,
                         batch_correction_method = NULL,
                         batch_labels = NULL,
                         max_n_batch = NULL,
                         use_assay = NULL,
                         use_slot = NULL,
                         ArchR_matrix = NULL,
                         countsplit = NULL,
                         countsplit_suffix = NULL,
                         input_matrix = NULL,
                         nn_matrix = NULL,
                         reduction = NULL,
                         n_cores = NULL,
                         random_seed = NULL,
                         verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Retrieve/check parameter input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "combineTrees")
  .validInput(subtree_list, "subtree_list", object)
  .validInput(key, "key", list("combineTrees", object))
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }
  # Random seed reproducibility
  if (n_cores > 1) {
    RNGkind("L'Ecuyer-CMRG")
  }

  # Get number of subtrees
  n_subtrees <- length(subtree_list)

  # If certain parameter inputs are not newly provided..
  # Retrieve parameter values from buildParentTree and subtree pruneTree records
  # If things don't match up, throw an error
  parameters <- list(
    alpha = alpha,
    p_adjust = p_adjust,
    feature_set = feature_set,
    exclude_features = exclude_features,
    n_iterations = n_iterations,
    n_trees = n_trees,
    use_variance = use_variance,
    min_accuracy = min_accuracy,
    min_connections = min_connections,
    max_repeat_errors = max_repeat_errors,
    distance_awareness = distance_awareness,
    collect_all_metrics = collect_all_metrics,
    sample_max = sample_max,
    downsampling_rate = downsampling_rate,
    min_reads = min_reads,
    normalization_method = normalization_method,
    batch_correction_method = batch_correction_method,
    batch_labels = batch_labels,
    max_n_batch = max_n_batch,
    use_assay = use_assay,
    use_slot = use_slot,
    ArchR_matrix = ArchR_matrix,
    countsplit = countsplit,
    countsplit_suffix = countsplit_suffix,
    random_seed = random_seed
  )
  parameters_to_check <- names(parameters)[sapply(parameters, is.null)]

  if (length(parameters_to_check) > 0) {
    # Retrieve parameter values from buildParentTree
    buildParentTree_parameters <- .retrieveData(object,
                                                key,
                                                "parameters",
                                                "buildParentTree_parameters")
    buildParentTree_parameters <- buildParentTree_parameters[intersect(names(buildParentTree_parameters),
                                                                       parameters_to_check)]
    # Retrieve parameter values from buildParentTree and subtree pruneTree records
    subtree_pruneTree_parameters <- subtree_list[[1]]$parameters$pruneTree_parameters
    subtree_pruneTree_parameters <- subtree_pruneTree_parameters[intersect(names(subtree_pruneTree_parameters),
                                                                           parameters_to_check)]

    if (identical(buildParentTree_parameters[intersect(names(subtree_pruneTree_parameters),
                                                       names(buildParentTree_parameters))],
                  subtree_pruneTree_parameters[intersect(names(subtree_pruneTree_parameters),
                                                         names(buildParentTree_parameters))])) {
      for (s in 2:n_subtrees) {
        subtree_pruneTree_parameters_s <- subtree_list[[s]]$parameters$pruneTree_parameters
        subtree_pruneTree_parameters_s <- subtree_pruneTree_parameters_s[intersect(names(subtree_pruneTree_parameters_s),
                                                                                   parameters_to_check)]
        if (!identical(subtree_pruneTree_parameters, subtree_pruneTree_parameters_s)) {
          stop(paste0("Parameter records for all provided subtrees must be identical. Otherwise, please supply input to relevant parameters: ",
                      paste(parameters_to_check, collapse = ", ")))
        }
      }
    } else {
      stop(paste0("Parameter records for buildParentTree and for all provided subtrees must be identical. Otherwise, please supply input to relevant parameters: ",
                  paste(names(buildParentTree_parameters), collapse = ", ")))
    }

    # Set parameters
    if ("alpha" %in% parameters_to_check) alpha <- subtree_pruneTree_parameters$alpha
    if ("p_adjust" %in% parameters_to_check) p_adjust <- subtree_pruneTree_parameters$p_adjust
    if ("feature_set" %in% parameters_to_check) feature_set <- subtree_pruneTree_parameters$feature_set
    if ("exclude_features" %in% parameters_to_check) exclude_features <- subtree_pruneTree_parameters$exclude_features
    if ("n_iterations" %in% parameters_to_check) n_iterations <- subtree_pruneTree_parameters$n_iterations
    if ("n_trees" %in% parameters_to_check) n_trees <- subtree_pruneTree_parameters$n_trees
    if ("use_variance" %in% parameters_to_check) use_variance <- subtree_pruneTree_parameters$use_variance
    if ("min_accuracy" %in% parameters_to_check) min_accuracy <- subtree_pruneTree_parameters$min_accuracy
    if ("min_connections" %in% parameters_to_check) min_connections <- subtree_pruneTree_parameters$min_connections
    if ("max_repeat_errors" %in% parameters_to_check) max_repeat_errors <- subtree_pruneTree_parameters$max_repeat_errors
    if ("distance_awareness" %in% parameters_to_check) distance_awareness <- subtree_pruneTree_parameters$distance_awareness
    if ("collect_all_metrics" %in% parameters_to_check) collect_all_metrics <- subtree_pruneTree_parameters$collect_all_metrics
    if ("sample_max" %in% parameters_to_check) sample_max <- subtree_pruneTree_parameters$sample_max
    if ("downsampling_rate" %in% parameters_to_check) downsampling_rate <- subtree_pruneTree_parameters$downsampling_rate
    if ("min_reads" %in% parameters_to_check) min_reads <- subtree_pruneTree_parameters$min_reads
    if ("normalization_method" %in% parameters_to_check) normalization_method <- subtree_pruneTree_parameters$normalization_method
    if ("batch_correction_method" %in% parameters_to_check) batch_correction_method <- subtree_pruneTree_parameters$batch_correction_method
    if ("batch_labels" %in% parameters_to_check) batch_labels <- subtree_pruneTree_parameters$batch_labels
    if ("max_n_batch" %in% parameters_to_check) max_n_batch <- subtree_pruneTree_parameters$max_n_batch
    if ("use_assay" %in% parameters_to_check) use_assay <- buildParentTree_parameters$use_assay
    if ("use_slot" %in% parameters_to_check) use_slot <- buildParentTree_parameters$use_slot
    if ("ArchR_matrix" %in% parameters_to_check) ArchR_matrix <- buildParentTree_parameters$ArchR_matrix
    if ("countsplit" %in% parameters_to_check) countsplit <- subtree_pruneTree_parameters$countsplit
    if ("countsplit_suffix" %in% parameters_to_check) countsplit_suffix <- subtree_pruneTree_parameters$countsplit_suffix
    if ("random_seed" %in% parameters_to_check) random_seed <- subtree_pruneTree_parameters$random_seed
  }

  # Verify parameter validity
  .validInput(alpha, "alpha")
  .validInput(p_adjust, "p_adjust")
  .validInput(feature_set, "feature_set", list("combineTrees", normalization_method))
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
  .validInput(max_n_batch, "max_n_batch")
  .validInput(countsplit, "countsplit")
  .validInput(countsplit_suffix, "countsplit_suffix", countsplit)
  .validInput(use_assay, "use_assay", list(object, countsplit, countsplit_suffix))
  .validInput(use_slot, "use_slot", list(object, use_assay, countsplit, countsplit_suffix))
  .validInput(ArchR_matrix, "ArchR_matrix", list(object, countsplit, countsplit_suffix))
  .validInput(random_seed, "random_seed")

  # Extract cell IDs
  cell_IDs <- .getCellIDs(object, use_assay)

  .validInput(input_matrix, "input_matrix", cell_IDs)
  .validInput(nn_matrix, "nn_matrix", cell_IDs)
  .validInput(reduction, "reduction", list("combineTrees", object))

  # User supplied data should not be mixed with data from prior CHOIR functions
  if (distance_awareness < Inf) {
    if (any(c(!is.null(input_matrix), !is.null(nn_matrix), !is.null(reduction))) &
        !all(c(!is.null(input_matrix), !is.null(nn_matrix), !is.null(reduction)))) {
      warning("If user-supplied data is provided for any of 'input_matrix', 'nn_matrix', or 'reduction', it should be provided for all three.")
    }
  } else {
    if (any(c(!is.null(input_matrix), !is.null(nn_matrix))) &
        !all(c(!is.null(input_matrix), !is.null(nn_matrix)))) {
      warning("If user-supplied data is provided for either 'input_matrix' or 'nn_matrix', it should be provided for both.")
    }
  }
  # User supplied data cannot be used with countsplitting
  if (!is.null(input_matrix) & countsplit == TRUE) {
    warning("Countsplitting is not currently compatible with user-supplied data for 'input_matrix'. Parameter 'count_split' set to FALSE.")
    countsplit <- FALSE
    countsplit_suffix <- NULL
  }

  # ---------------------------------------------------------------------------
  # Prepare object
  # ---------------------------------------------------------------------------

  # Track progress
  if (p_adjust == "none") {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 1/5) Preparing object..")
  } else {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 1/5) Preparing object and calculating compiled adjusted significance threshold..")
  }

  if (p_adjust == "bonferroni") {
    # Get total number of comparisons
    n_total_comparisons <- 0
    for (s in 1:n_subtrees) {
      subtree_s <- subtree_list[[s]]
      # Check if this subtree was run or not
      # If so, collect number of comparisons
      if ("records" %in% names(subtree_s)) {
        n_comparisons <- length(subtree_s$records$comparison_records$percentile_accuracy[!is.na(subtree_s$records$comparison_records$percentile_accuracy)])
        n_total_comparisons <- n_total_comparisons + n_comparisons
      }
    }
    # Set new adjusted alpha
    adjusted_alpha <- alpha/n_total_comparisons
  } else if (p_adjust == "fdr") {
    stop("FDR multiple comparison correction is not yet supported within this function. Please set parameter p_adjust to 'bonferroni' or 'none' to proceed.")
  } else if (p_adjust == "none") {
    adjusted_alpha <- alpha
  }

  # Extract nearest neighbor matrix/matrices
  if (!is.null(nn_matrix)) {
    nn_matrix_provided <- TRUE
  } else if (!is.null(buildParentTree_parameters)) {
    nn_matrix_provided <- FALSE
    nn_matrix <- .retrieveData(object, key, "graph", "P0_graph_nn")
  } else {
    stop("No nearest neighbor adjacency matrix provided.")
  }

  # Retrieve parent tree
  cluster_tree <- .retrieveData(object,
                                key,
                                "clusters",
                                "P0_tree")
  n_levels_parent_tree <- ncol(cluster_tree)

  # ---------------------------------------------------------------------------
  # Merge clusters at new adjusted alpha and compile results
  # ---------------------------------------------------------------------------

  if (p_adjust != "none") {
    # Track progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/5) Merge clusters within each subtree that do not pass new adjusted significance threshold..")

    # Collect condensed comparison records
    all_metrics <- c('comparison', 'cluster1_size', 'cluster2_size', 'sample_size',
                     'mean_accuracy', 'var_accuracy', 'mean_errors',
                     'mean_permuted_accuracy', 'var_permuted_accuracy',
                     'percentile_accuracy', 'percentile_variance',
                     'n_repeat_errors1', 'n_repeat_errors2',
                     'mean_repeat_errors1', 'mean_repeat_errors2',
                     'mean_modified_accuracy', 'var_modified_accuracy',
                     'percentile_modified_accuracy', 'percentile_modified_variance',
                     'batches_used', 'batch_mean_accuracies',
                     'connectivity', 'root_distance', 'subtree_distance', 'time',
                     'decision')
    selected_metrics <- all_metrics[c(1:11,
                                      `if`(collect_all_metrics == TRUE | max_repeat_errors > 0, 12:15, NULL),
                                      `if`(max_repeat_errors > 0, 16:19, NULL),
                                      `if`(batch_correction_method == "Harmony", 20:21, NULL),
                                      `if`(collect_all_metrics == TRUE | min_connections > 0, 22, NULL),
                                      `if`(distance_awareness < Inf, 23:24, NULL),
                                      25:26)]
    comparison_records <- data.frame(matrix(ncol = length(selected_metrics), nrow = 0))
    colnames(comparison_records) <- selected_metrics

    # Collect all subtree cluster labels
    all_cluster_ids <- data.frame(Cell_ID = NULL, subtree_cluster = NULL, merged_subtree_cluster = NULL, parent_cluster = NULL)

    for (s in 1:n_subtrees) {
      subtree_s <- subtree_list[[s]]
      # Check if this subtree has more than 1 cluster
      if (dplyr::n_distinct(subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]][,paste0("CHOIR_clusters_", alpha)]) > 1) {
        # If so, check for necessary merges
        correction_check <- subtree_s$records$comparison_records %>%
          dplyr::mutate(correct = ifelse(decision == "split" &
                                           (percentile_accuracy >= adjusted_alpha |
                                              percentile_variance >= adjusted_alpha), TRUE,
                                         ifelse(decision == "split: repeat error" &
                                                  (percentile_modified_accuracy >= adjusted_alpha |
                                                     percentile_modified_variance >= adjusted_alpha), TRUE, FALSE))) %>%
          dplyr::filter(correct == TRUE) %>%
          tidyr::separate_wider_delim(comparison, delim = " vs. ", names = c("cluster1", "cluster2")) %>%
          dplyr::filter(cluster1 %in% unique(subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$Record_cluster_label),
                        cluster2 %in% unique(subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$Record_cluster_label))

        compiled_cluster_labels <- unique(subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$Record_cluster_label)
        child_IDs <- subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$Record_cluster_label
        subtree_cluster_labels <- paste0("P", s, "_", "L", n_levels_parent_tree + 2, "_", as.numeric(as.factor(child_IDs)))

        if (nrow(correction_check) > 0) {
          # Check whether comparisons should remain unmerged, from bottom up
          for (comp in 1:nrow(correction_check)) {
            # Check whether these clusters appear elsewhere in the corrections
            if (sum(c(correction_check$cluster1, correction_check$cluster2) == correction_check$cluster1[comp]) == 1 &
                sum(c(correction_check$cluster1, correction_check$cluster2) == correction_check$cluster2[comp]) == 1) {
              # If they do not, merge these two clusters
              merge_pair <- correction_check[comp,]
            } else {
              # If these clusters do appear elsewhere in the corrections that are among the current cluster IDs,
              # merge the cluster pair with the closest distance
              if (distance_awareness < Inf) {
                subtree_distances <- dplyr::filter(correction_check,
                                                   cluster1 %in% c(correction_check$cluster1[comp],
                                                                   correction_check$cluster2[comp]) |
                                                     cluster2 %in% c(correction_check$cluster1[comp],
                                                                     correction_check$cluster2[comp]))$subtree_distance
                if (length(subtree_distances) > 0 & any(!is.na(subtree_distances))) {
                  min_subtree_distance <- min(subtree_distances, na.rm = TRUE)
                } else {
                  min_subtree_distance <- NA
                }
                if (!is.na(min_subtree_distance)) {
                  # Use subtree distance if available
                  merge_pair <- correction_check %>%
                    dplyr::filter(cluster1 %in% c(correction_check$cluster1[comp],
                                                  correction_check$cluster2[comp]) |
                                    cluster2 %in% c(correction_check$cluster1[comp],
                                                    correction_check$cluster2[comp]),
                                  subtree_distance == min_subtree_distance)
                } else {
                  # Otherwise use root distance
                  root_distances <- dplyr::filter(correction_check,
                                                  cluster1 %in% c(correction_check$cluster1[comp],
                                                                  correction_check$cluster2[comp]) |
                                                    cluster2 %in% c(correction_check$cluster1[comp],
                                                                    correction_check$cluster2[comp]))$root_distance
                  min_root_distance <- min(root_distances, na.rm = TRUE)
                  merge_pair <- correction_check %>%
                    dplyr::filter(cluster1 %in% c(correction_check$cluster1[comp],
                                                  correction_check$cluster2[comp]) |
                                    cluster2 %in% c(correction_check$cluster1[comp],
                                                    correction_check$cluster2[comp]),
                                  root_distance == min_root_distance)
                }
              } else {
                # Alternately, use accuracy scores
                min_accuracy <- min(dplyr::filter(correction_check,
                                                  cluster1 %in% c(correction_check$cluster1[comp],
                                                                  correction_check$cluster2[comp]) |
                                                    cluster2 %in% c(correction_check$cluster1[comp],
                                                                    correction_check$cluster2[comp]))$mean_accuracy,
                                    na.rm = TRUE)
                merge_pair <- correction_check %>%
                  dplyr::filter(cluster1 %in% c(correction_check$cluster1[comp],
                                                correction_check$cluster2[comp]) |
                                  cluster2 %in% c(correction_check$cluster1[comp],
                                                  correction_check$cluster2[comp]),
                                mean_accuracy == min_accuracy)
              }
            }
            new_labels_list <- .getNewLabels(merge_groups = list(c(merge_pair$cluster1[1],
                                                                   merge_pair$cluster2[1])),
                                             level = 0,
                                             compiled_labels = compiled_cluster_labels)
            compiled_cluster_labels <- new_labels_list[["compiled_cluster_labels"]]
            merged_label <- new_labels_list[["merge_group_labels"]][[1]]

            child_IDs[child_IDs %in% c(merge_pair$cluster1[1], merge_pair$cluster2[1])] <- merged_label
            # Update records
            subtree_s$records$comparison_records$decision[which(subtree_s$records$comparison_records$comparison ==
                                                                  paste0(merge_pair$cluster1[1], " vs. ",
                                                                         merge_pair$cluster2[1]))] <- "merge: adjustment"
          }
        }
        # Filter records
        condensed_records <- subtree_s$records$comparison_records %>% data.frame() %>%
          tidyr::separate_wider_delim(comparison, delim = " vs. ", names = c("cluster1", "cluster2"),
                                      cols_remove = FALSE) %>%
          dplyr::filter(cluster1 %in% child_IDs,
                        cluster2 %in% child_IDs)

        # Provide new cluster labels
        merged_subtree_cluster_labels <- paste0("P", s, "_", "L", n_levels_parent_tree + 1, "_", as.numeric(as.factor(child_IDs)))
        all_cluster_ids <- rbind(all_cluster_ids, data.frame(CellID = subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$CellID,
                                                             subtree_cluster = subtree_cluster_labels,
                                                             merged_subtree_cluster = merged_subtree_cluster_labels,
                                                             parent_cluster = s))
        for (r in 1:nrow(condensed_records)) {
          cluster1_new_label <- merged_subtree_cluster_labels[child_IDs == condensed_records$cluster1[r]][1]
          cluster2_new_label <- merged_subtree_cluster_labels[child_IDs == condensed_records$cluster2[r]][1]
          condensed_records$cluster1[r] <- cluster1_new_label
          condensed_records$cluster2[r] <- cluster2_new_label
          condensed_records$comparison[r] <- paste0(condensed_records$cluster1[r], " vs. ", condensed_records$cluster2[r])
        }
        comparison_records <- rbind(comparison_records, condensed_records[,-c(1,2)])
      } else {
        # Provide new cluster labels
        child_IDs <- subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$Record_cluster_label
        subtree_cluster_labels <- paste0("P", s, "_", "L", n_levels_parent_tree + 2, "_", as.numeric(as.factor(child_IDs)))
        merged_subtree_cluster_labels <- paste0("P", s, "_", "L", n_levels_parent_tree + 1, "_", as.numeric(as.factor(child_IDs)))
        all_cluster_ids <- rbind(all_cluster_ids, data.frame(CellID = subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$CellID,
                                                             subtree_cluster = subtree_cluster_labels,
                                                             merged_subtree_cluster = merged_subtree_cluster_labels,
                                                             parent_cluster = s))
      }
    }
  } else {
    # Track progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/5) Standardize cluster labels and compile records..")

    # Provide new cluster labels
    child_IDs <- subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$Record_cluster_label
    subtree_cluster_labels <- paste0("P", s, "_", "L", n_levels_parent_tree + 2, "_", as.numeric(as.factor(child_IDs)))
    merged_subtree_cluster_labels <- paste0("P", s, "_", "L", n_levels_parent_tree + 1, "_", as.numeric(as.factor(child_IDs)))
    all_cluster_ids <- rbind(all_cluster_ids, data.frame(CellID = subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]]$CellID,
                                                         subtree_cluster = subtree_cluster_labels,
                                                         merged_subtree_cluster = merged_subtree_cluster_labels,
                                                         parent_cluster = s))
    # If more than 1 cluster
    if (dplyr::n_distinct(subtree_s$clusters[paste0("CHOIR_clusters_", alpha)][[1]][,paste0("CHOIR_clusters_", alpha)]) > 1) {
      # Filter records
      condensed_records <- subtree_s$records$comparison_records %>% data.frame() %>%
        tidyr::separate_wider_delim(comparison, delim = " vs. ", names = c("cluster1", "cluster2"),
                                    cols_remove = FALSE) %>%
        dplyr::filter(cluster1 %in% child_IDs,
                      cluster2 %in% child_IDs)
      for (r in 1:nrow(condensed_records)) {
        cluster1_new_label <- merged_subtree_cluster_labels[child_IDs == condensed_records$cluster1[r]][1]
        cluster2_new_label <- merged_subtree_cluster_labels[child_IDs == condensed_records$cluster2[r]][1]
        condensed_records$cluster1[r] <- cluster1_new_label
        condensed_records$cluster2[r] <- cluster2_new_label
        condensed_records$comparison[r] <- paste0(condensed_records$cluster1[r], " vs. ", condensed_records$cluster2[r])
      }
      comparison_records <- rbind(comparison_records, condensed_records[,-c(1,2)])
    }
  }

  # Create new clustering tree
  # Track progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 3/5) Create compiled clustering tree and identify remaining comparisons..")

  # Add new level to parent tree with current set of clusters and number this new level
  rownames(all_cluster_ids) <- all_cluster_ids$CellID
  all_cluster_ids <- all_cluster_ids[cell_IDs,]
  cluster_tree$L_new <- all_cluster_ids$merged_subtree_cluster
  colnames(cluster_tree)[ncol(cluster_tree)] <- paste0("L", ncol(cluster_tree))

  # Record of stepwise cluster IDs
  stepwise_cluster_IDs <- data.frame(CellID = cell_IDs,
                                     subtree_clusters = all_cluster_ids$subtree_cluster,
                                     merged_subtree_clusters = all_cluster_ids$merged_subtree_cluster)
  colnames(stepwise_cluster_IDs)[2] <- paste0("stepwise_cluster_ID_", alpha, "_L", (n_levels_parent_tree + 2))
  colnames(stepwise_cluster_IDs)[3] <- paste0("stepwise_cluster_ID_", alpha, "_L", (n_levels_parent_tree + 1))

  # Identify remaining comparisons

  # Get cluster distances
  child_IDs <- cluster_tree[,ncol(cluster_tree)]

  P0_reduction <- .retrieveData(object,
                                key,
                                "reduction",
                                "P0_reduction")
  root_distances <- .getCentroidDistance(P0_reduction,
                                         child_IDs)

  # Set distance records
  distance_records <- data.frame(cluster_name = unique(child_IDs),
                                 min_root_distance = NA,
                                 min_subtree_distance = NA,
                                 max_pval = NA,
                                 parent_cluster = stringr::str_extract(unique(child_IDs), "P\\d*"))

  min_root_distance_list <- parallel::mclapply(1:nrow(distance_records), FUN = function(i) {
    cluster_i <- distance_records$cluster_name[i]
    parent_i <- distance_records$parent_cluster[i]
    set <- distance_records$cluster_name[distance_records$parent_cluster == parent_i]
    if (length(set) > 1) {
      set <- set[set != cluster_i]
      min(root_distances[cluster_i, set])
    } else {
      NA
    }
  },
  mc.cores = n_cores,
  mc.set.seed = TRUE)
  distance_records$min_root_distance <- unlist(min_root_distance_list)
  distance_records$min_subtree_distance <- NA
  distance_records$max_pval <- NA
  distance_records <- distance_records %>%
    dplyr::select(-parent_cluster) %>%
    dplyr::filter(!is.na(min_root_distance))

  # Create list of permitted comparisons
  permitted_comparisons <- .getPermittedComparisons(object = object,
                                                    key = key,
                                                    type = "combineTrees",
                                                    clusters = child_IDs,
                                                    cell_IDs = cell_IDs,
                                                    nn_matrices = list("P0" = nn_matrix),
                                                    min_connections = min_connections,
                                                    collect_all_metrics = collect_all_metrics,
                                                    distance_awareness = distance_awareness,
                                                    root_distances = root_distances,
                                                    n_cores = n_cores)
  if (distance_awareness < Inf) {
    # Filter permitted comparisons against distance records
    permitted_comparisons <- permitted_comparisons %>%
      dplyr::left_join(distance_records, by = c("cluster1" = "cluster_name")) %>%
      dplyr::rename(min_root_distance1 = min_root_distance) %>%
      dplyr::left_join(distance_records, by = c("cluster2" = "cluster_name")) %>%
      dplyr::rename(min_root_distance2 = min_root_distance) %>%
      dplyr::filter(is.na(min_root_distance1) |
                      is.na(min_root_distance2) |
                      (root_distance > distance_awareness * min_root_distance1 & root_distance > distance_awareness * min_root_distance2)) %>%
      dplyr::select(cluster1, cluster2, tree_name, connectivity, root_distance, subtree_distance)
  }

  # -------------------------------------------------------------------------
  # Prepare compiled tree for pruning
  # -------------------------------------------------------------------------

  # Track progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 4/5) Prepare compiled tree for pruning..")

  # ---------------------------------------------------------------------------
  # Set values for countsplitting
  # ---------------------------------------------------------------------------

  if (countsplit == TRUE) {
    if (is.null(countsplit_suffix)) {
      countsplit_suffix <- c("_1", "_2")
    }
  } else {
    countsplit_suffix <- c("", "")
  }
  if (!is.null(input_matrix)) {
    n_modalities <- 1
  } else {
    if (methods::is(object, "ArchRProject")) {
      n_modalities <- max(length(ArchR_matrix), 1)
    } else {
      n_modalities <- max(length(use_assay), 1)
    }
  }
  # Set new values
  if (methods::is(object, "Seurat")) {
    # Seurat object
    if (is.null(use_assay)) {
      if ("Assay5" %in% methods::is(object[[Seurat::DefaultAssay(object)]])) {
        seurat_version <- "v5"
      } else if ("Assay" %in% methods::is(object[[Seurat::DefaultAssay(object)]])) {
        seurat_version <- "v4"
      } else {
        stop("Assay '", Seurat::DefaultAssay(object),
             "' provided for parameter 'use_assay' is not of class 'Assay' or 'Assay5', please supply valid input!")
      }
    } else {
      if ("Assay5" %in% methods::is(object[[use_assay[1]]])) {
        seurat_version <- "v5"
      } else if ("Assay" %in% methods::is(object[[use_assay[1]]])) {
        seurat_version <- "v4"
      } else {
        stop("Assay '", use_assay[1],
             "' provided for parameter 'use_assay' is not of class 'Assay' or 'Assay5', please supply valid input!")
      }
    }
    # Set values of 'use_assay' and 'use_slot' if necessary
    if (is.null(use_assay)) {
      use_assay <- Seurat::DefaultAssay(object)
    }
    if (is.null(use_slot)) {
      if (use_assay %in% c("RNA", "sketch")) {
        use_slot <- "data"
      } else if (use_assay == "SCT" | use_assay == "integrated") {
        use_slot <- "scale.data"
      } else {
        stop("When using a non-standard assay in a Seurat object, please supply a valid input for the slot parameter.")
      }
    }
    # Set new values
    use_assay_build <- use_assay
    use_assay_prune <- use_assay
    use_slot_build <- paste0(use_slot, countsplit_suffix[1])
    use_slot_prune <- paste0(use_slot, countsplit_suffix[2])
    ArchR_matrix_build <- NULL
    ArchR_matrix_prune <- NULL
    countsplit_text <- paste0("\n - Assay: ", use_assay,
                              "\n - ",
                              ifelse(seurat_version == "v5", "Layer", "Slot"),
                              ifelse(n_modalities == 1, " ", "s "),
                              "used to build tree: ",
                              paste(use_slot_build, collapse = ", "),
                              "\n - ",
                              ifelse(seurat_version == "v5", "Layer", "Slot"),
                              ifelse(n_modalities == 1, " ", "s "),
                              "used to prune tree: ",
                              paste(use_slot_prune, collapse = ", "))
  } else if (methods::is(object, "SingleCellExperiment")) {
    # SingleCellExperiment object
    # Set value of 'use_assay' if necessary
    if (is.null(use_assay)) {
      use_assay <- "logcounts"
    }
    # Set new values
    use_assay_build <- paste0(use_assay, countsplit_suffix[1])
    use_assay_prune <- paste0(use_assay, countsplit_suffix[2])
    use_slot_build <- NULL
    use_slot_prune <- NULL
    ArchR_matrix_build <- NULL
    ArchR_matrix_prune <- NULL
    countsplit_text <- paste0("\n - Assay",
                              ifelse(n_modalities == 1, " ", "s "),
                              "used to build tree: ",
                              paste(use_assay_build, collapse = ", "),
                              "\n - Assay",
                              ifelse(n_modalities == 1, " ", "s "),
                              "used to prune tree: ",
                              paste(use_assay_prune, collapse = ", "))
  } else if (methods::is(object, "ArchRProject")) {
    # ArchR object
    # Set value of 'ArchR_matrix' if necessary
    if (is.null(ArchR_matrix)) {
      ArchR_matrix <- "GeneScoreMatrix"
    }
    if (countsplit == TRUE) {
      warning("Count splitting has not been tested thoroughly outside the context of RNA-seq data.")
    }
    # Set new values
    use_assay_build <- NULL
    use_assay_prune <- NULL
    use_slot_build <- NULL
    use_slot_prune <- NULL
    ArchR_matrix_build <- paste0(ArchR_matrix, countsplit_suffix[1])
    ArchR_matrix_prune <- paste0(ArchR_matrix, countsplit_suffix[2])
    countsplit_text <- paste0("\n - ArchR matri",
                              ifelse(n_modalities == 1, "x ", "ces "),
                              "used to build tree: ",
                              paste(ArchR_matrix_build, collapse = ", "),
                              "\n - ArchR matri",
                              ifelse(n_modalities == 1,  "x ", "ces "),
                              "used to prune tree: ",
                              paste(ArchR_matrix_prune, collapse = ", "))
  }

  # Extract input matrix
  if (!is.null(input_matrix)) {
    input_matrix_provided <- TRUE
    input_matrix <- .getMatrix(use_matrix = input_matrix,
                               exclude_features = exclude_features,
                               verbose = FALSE)
    input_matrix <- BiocGenerics::t(input_matrix)
    n_modalities <- 1
    # Feature names
    features <- colnames(input_matrix)
    # Clean up
    rm(input_matrix)
    # Object type
    if (methods::is(object, "ArchRProject")) {
      object_type <- "ArchRProject"
      .requirePackage("ArchR", installInfo = "Instructions at archrproject.com")
    } else {
      if (methods::is(object, "Seurat")) {
        object_type <- "Seurat"
      } else {
        object_type <- "SingleCellExperiment"
        .requirePackage("SingleCellExperiment", source = "bioc")
      }
    }
    # Need n_modalities to validate some inputs
    .validInput(normalization_method, "normalization_method", list(object, n_modalities, use_assay))
    .validInput(batch_correction_method, "batch_correction_method", n_modalities)
    .validInput(batch_labels, "batch_labels", object)
  } else if (!is.null(buildParentTree_parameters)) {
    input_matrix_provided <- FALSE
    # Number of modalities & object type
    if (methods::is(object, "ArchRProject")) {
      n_modalities <- max(length(ArchR_matrix), 1)
      object_type <- "ArchRProject"
      .requirePackage("ArchR", installInfo = "Instructions at archrproject.com")
    } else {
      n_modalities <- max(length(use_assay), 1)
      if (methods::is(object, "Seurat")) {
        object_type <- "Seurat"
      } else {
        object_type <- "SingleCellExperiment"
        .requirePackage("SingleCellExperiment", source = "bioc")
      }
    }
    # Need n_modalities to validate some inputs
    .validInput(normalization_method, "normalization_method", list(object, n_modalities, use_assay))
    .validInput(batch_correction_method, "batch_correction_method", n_modalities)
    .validInput(batch_labels, "batch_labels", object)

    # Initialize list of feature names
    features <- c()

    # Extract data
    # Use all / var features
    if (feature_set == "all") {
      use_features <- NULL
    } else {
      use_features <- .retrieveData(object, key, "var_features", "P0_var_features")
    }
    # Extract input matrix for random forest comparisons
    if (n_modalities > 1) {
      input_matrix_list <- vector("list", n_modalities)
      for (m in 1:n_modalities) {
        # Match input arguments
        use_assay_prune_m <- .matchArg(use_assay_prune, m)
        use_slot_prune_m <- .matchArg(use_slot_prune, m)
        ArchR_matrix_prune_m <- .matchArg(ArchR_matrix_prune, m)
        normalization_method_m <- .matchArg(normalization_method, m)
        # Features
        if (length(use_features) > 1) {
          use_features_m <- use_features[[m]]
        } else {
          use_features_m <- use_features
        }
        input_matrix_list[[m]] <- .getMatrix(object = object,
                                             use_assay = use_assay_prune_m,
                                             use_slot = use_slot_prune_m,
                                             ArchR_matrix = ArchR_matrix_prune_m,
                                             use_features = use_features_m,
                                             exclude_features = exclude_features,
                                             use_cells = NULL,
                                             verbose = FALSE)
        if (normalization_method_m == "SCTransform") {
          input_matrix_list[[m]] <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(input_matrix_list[[m]]),
                                                                         return.only.var.genes = FALSE,
                                                                         seed.use = random_seed,
                                                                         verbose = FALSE)@assays$SCT@scale.data)
        }
      }
      input_matrix <- do.call(rbind, input_matrix_list)
      input_matrix <- BiocGenerics::t(input_matrix)
      # Clean up
      rm(use_features_m)
    } else {
      input_matrix <- .getMatrix(object = object,
                                 use_assay = use_assay_prune,
                                 use_slot = use_slot_prune,
                                 ArchR_matrix = ArchR_matrix_prune,
                                 use_features = use_features,
                                 exclude_features = exclude_features,
                                 use_cells = NULL,
                                 verbose = FALSE)
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
    # Clean up
    rm(use_features)
  } else {
    stop("No 'input_matrix' supplied. Please supply valid input!")
  }

  # Reduction
  if (!is.null(reduction)) {
    reduction_provided <- TRUE
  } else {
    reduction_provided <- FALSE
  }

  # If 'Harmony' batch correction was used, extract batch labels
  if (batch_correction_method == "Harmony") {
    batches <- .retrieveData(object = object, key = key, type = "cell_metadata", name = batch_labels)
    names(batches) <- cell_IDs
  }

  # Set downsampling rate
  if (downsampling_rate == "auto") {
    downsampling_rate <- min(1, (1/2)^(log10(length(cell_IDs)/5000)))
  }

  # Number of starting clusters
  n_starting_clusters <- dplyr::n_distinct(cluster_tree[, ncol(cluster_tree)])
  # Number of levels in cluster tree
  n_levels <- ncol(cluster_tree)

  # Provided input
  provided_input <- c(
    if (input_matrix_provided) "input_matrix",
    if (nn_matrix_provided) "nn_matrix",
    if (reduction_provided) "reduction"
  )

  # Report object & parameter details
  if (verbose) message("\nInput data:",
                       "\n - Object type: ", object_type,
                       `if`(length(provided_input) > 0, paste0("\n - Provided inputs: ", paste(provided_input, collapse = ", ")), ""),
                       "\n - # of cells: ", length(cell_IDs),
                       "\n - # of batches: ", `if`(batch_correction_method == "none", 1, dplyr::n_distinct(batches)),
                       "\n - # of modalities: ", n_modalities,
                       "\n - # of subtrees: ", n_subtrees,
                       "\n - # of levels: ", n_levels,
                       "\n - # of starting clusters: ", n_starting_clusters,
                       "\n - Countsplitting: ", countsplit,
                       countsplit_text)
  if (verbose) message("\nProceeding with the following parameters:",
                       "\n - Intermediate data stored under key: ", key,
                       "\n - Alpha: ", alpha,
                       "\n - Multiple comparison adjustment: ", p_adjust,
                       "\n - Features to train RF: ", feature_set,
                       "\n - # of excluded features: ", length(exclude_features),
                       "\n - # of permutations: ", n_iterations,
                       "\n - # of RF trees: ", n_trees,
                       "\n - Use variance: ", use_variance,
                       "\n - Minimum accuracy: ", min_accuracy,
                       "\n - Minimum connections: ", min_connections,
                       "\n - Maximum repeated errors: ", max_repeat_errors,
                       "\n - Distance awareness: ", distance_awareness,
                       "\n - All metrics collected: ", collect_all_metrics,
                       "\n - Maximum cells sampled: ", sample_max,
                       "\n - Downsampling rate: ", round(downsampling_rate, 4),
                       "\n - Minimum reads: ", `if`(is.null(min_reads),
                                                    paste0(">0 reads"), paste0(">1 read per ", min_reads, " cells")),
                       "\n - Normalization method: ", normalization_method,
                       "\n - Batch correction method: ", batch_correction_method,
                       `if`(batch_correction_method != 'none', paste0("\n - Metadata column containing batch information: ", batch_labels), ""),
                       `if`(batch_correction_method != 'none', paste0("\n - Maximum # of batches used per permutation test: ", max_n_batch), ""),
                       "\n - # of cores: ", n_cores,
                       "\n - Random seed: ", random_seed,
                       "\n")

  # ---------------------------------------------------------------------------
  # Initialize output data structures
  # ---------------------------------------------------------------------------

  # Data frame to track comparisons to avoid redundancy
  all_metrics <- c('comparison', 'cluster1_size', 'cluster2_size', 'sample_size',
                   'mean_accuracy', 'var_accuracy', 'mean_errors',
                   'mean_permuted_accuracy', 'var_permuted_accuracy',
                   'percentile_accuracy', 'percentile_variance',
                   'n_repeat_errors1', 'n_repeat_errors2',
                   'mean_repeat_errors1', 'mean_repeat_errors2',
                   'mean_modified_accuracy', 'var_modified_accuracy',
                   'percentile_modified_accuracy', 'percentile_modified_variance',
                   'batches_used', 'batch_mean_accuracies', 'batch_var_accuracies',
                   'connectivity', 'root_distance', 'subtree_distance', 'time',
                   'decision')
  selected_metrics <- all_metrics[c(1:11,
                                    `if`(collect_all_metrics == TRUE | max_repeat_errors > 0, 12:15, NULL),
                                    `if`(max_repeat_errors > 0, 16:19, NULL),
                                    `if`(batch_correction_method == "Harmony", 20:22, NULL),
                                    `if`(collect_all_metrics == TRUE | min_connections > 0, 23, NULL),
                                    `if`(methods::is(distance_awareness, "numeric"), 24:25, NULL),
                                    26:27)]
  comparison_records <- data.frame(matrix(ncol = length(selected_metrics), nrow = 0))
  colnames(comparison_records) <- selected_metrics

  # Feature importance records
  feature_importance_records <- data.frame(matrix(ncol = (length(features)+2), nrow = 0))
  colnames(feature_importance_records) <- c('cluster1', 'cluster2', features)

  # Set of all cluster labels in original tree
  compiled_cluster_labels <- unlist(apply(cluster_tree, 2, unique))
  names(compiled_cluster_labels) <- NULL

  # -------------------------------------------------------------------------
  # Iterate through each level of the clustering tree (bottom-up)
  # -------------------------------------------------------------------------

  # Track progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 5/5) Iterating through clustering tree..")

  # Go through tree as usual but deem "split" if not in list of permitted comparisons
  # Update list as you go, removing comparisons as necessary
  # If clusters are merged, add new set of permitted comparisons

  # Get mean cluster size at each level of the tree
  mean_cluster_sizes <- apply(cluster_tree, 2, FUN = function(x) length(cell_IDs)/dplyr::n_distinct(x))
  level_weights <- rep(1, n_levels)
  level_weights[mean_cluster_sizes > 0.1*length(cell_IDs)] <- 5
  level_weights <- (level_weights/sum(level_weights))*95
  names(level_weights) <- paste0("L", seq(0, n_levels - 1))

  # Progress bar
  start_time <- Sys.time()
  hour_start_time <- Sys.time()
  pb <- progress::progress_bar$new(format = "Pruning tree..        [:bar] (:percent) in :elapsedfull",
                                   total = 100, clear = FALSE)
  pb$tick(0)
  percent_done <- 0

  # Starting values
  parent_IDs <- cluster_tree[, n_levels-1]
  child_IDs <- cluster_tree[, n_levels]
  n_current_clusters <- dplyr::n_distinct(child_IDs)
  # Start at bottom 2 levels of tree
  lvl <- n_levels-1
  # Complete?
  complete <- FALSE
  # Progress markers
  progress_markers <- c(10,20,30,40,50,60,70,80,90)

  while (complete == FALSE) {
    # Get all cluster IDs at this level
    unique_parent_IDs <- unique(parent_IDs)
    new_clusters <- c()

    print(paste0("Level ", lvl, ": ", length(unique_parent_IDs), " parent clusters"))

    # For each parent cluster
    for (parent in 1:length(unique_parent_IDs)) {
      # Current parent cluster
      parent_cluster <- unique_parent_IDs[parent]
      parent_inds <- which(parent_IDs == parent_cluster)
      # All child clusters of this parent
      unique_child_IDs <- unique(child_IDs[parent_inds])
      # If there is only one child, move on (this cluster branch does not split at this level of the tree)
      # If there are two or more child clusters, start to make comparisons
      n_child_clusters <- length(unique_child_IDs)
      print(paste0("Parent ", parent, ": ", n_child_clusters, " child clusters"))
      if (n_child_clusters > 1) {
        # Alternately, if this parent and its child clusters have not changed at all since the last level, skip
        if (lvl == n_levels-2) {
          check_identical <- identical(cluster_tree[parent_inds, ncol(cluster_tree)], child_IDs[parent_inds])
        } else if (lvl < n_levels-2) {
          check_identical <- identical(stepwise_cluster_IDs[parent_inds, (ncol(stepwise_cluster_IDs)-1)], child_IDs[parent_inds])
        } else {
          check_identical <- FALSE
        }
        if (dplyr::n_distinct(cluster_tree[parent_inds, lvl+1]) == 1 & # (Check that previous level was originally 1 cluster)
            check_identical == TRUE) {
          # Progress
          if (lvl >= 0) {
            tick_amount <- 0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
            pb$tick(tick_amount)
            if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                           (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
              hour_start_time <- Sys.time()
              pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                                round(difftime(Sys.time(), start_time, units = "min"), 2),
                                " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
            }
            percent_done <- percent_done + tick_amount
          } else {
            if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
              hour_start_time <- Sys.time()
              pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                " : Running additional comparisons, ",
                                round(difftime(Sys.time(), start_time, units = "min"), 2),
                                " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
            }
          }
        } else {
          # Create matrix for comparison results
          result_matrix <- matrix(NA, n_child_clusters, n_child_clusters)
          colnames(result_matrix) <- unique_child_IDs
          rownames(result_matrix) <- unique_child_IDs
          # For each child cluster
          for (child1 in 1:(n_child_clusters-1)) {
            # Name of child cluster 1
            child1_name <- unique_child_IDs[child1]
            # Check if child cluster 1 is part of any permitted comparisons
            if (child1_name %in% unique(c(permitted_comparisons$cluster1, permitted_comparisons$cluster2))) {
              # Get names of cells belonging to child1 cluster
              child1_cells <- cell_IDs[which(child_IDs == child1_name)]
              # If there is only one cell in child1, simply fill in the result_matrix with all "merge"
              if (length(child1_cells) == 1) {
                result_matrix[child1_name, ] <- "merge"
                result_matrix[, child1_name] <- "merge"
                result_matrix[child1_name, child1_name] <- NA
              } else {
                # Compare the child cluster pairwise to each subsequent sibling cluster
                for (child2 in (child1+1):n_child_clusters) {
                  child2_name <- unique_child_IDs[child2]
                  # Get names of cells belonging to child1 cluster
                  child2_cells <- cell_IDs[which(child_IDs == child2_name)]
                  # If there is only 1 cell in child2, automatically merge
                  if (length(child2_cells) == 1) {
                    result_matrix[child1_name, child2_name] <- "merge"
                    result_matrix[child2_name, child1_name] <- "merge"
                  } else {
                    # Check whether clusters are among permitted comparisons
                    if (nrow(dplyr::filter(permitted_comparisons, (cluster1 == child1_name & cluster2 == child2_name) |
                                           (cluster2 == child2_name & cluster1 == child1_name))) == 0) {
                      # Add result to matrix
                      result_matrix[child1_name, child2_name] <- "split"
                      result_matrix[child2_name, child1_name] <- "split"
                    } else {
                      comparison_start_time <- Sys.time()
                      # Check whether clusters have been previously compared
                      previous_comparison <- .checkComparisonRecords(cluster1_name = child1_name,
                                                                     cluster1_cells = child1_cells,
                                                                     cluster2_name = child2_name,
                                                                     cluster2_cells = child2_cells,
                                                                     comparison_records = comparison_records)
                      # If they have, fill in result matrix
                      if (previous_comparison[["previously_compared"]] == TRUE) {
                        # Add result to matrix
                        result_matrix[child1_name, child2_name] <- previous_comparison[["result"]]
                        result_matrix[child2_name, child1_name] <- previous_comparison[["result"]]
                      } else {
                        # Extract info from permitted comparison dataframe
                        P0_distance <- dplyr::filter(permitted_comparisons, (cluster1 == child1_name & cluster2 == child2_name) |
                                                       (cluster2 == child2_name & cluster1 == child1_name))$root_distance[1]
                        adjacent <- dplyr::filter(permitted_comparisons, (cluster1 == child1_name & cluster2 == child2_name) |
                                                    (cluster2 == child2_name & cluster1 == child1_name))$connectivity[1]
                        # Check whether distance between current clusters is greater than previous comparisons
                        distance_conflict <- FALSE
                        if (distance_awareness < Inf & !is.null(distance_records)) {
                          if (child1_name %in% distance_records$cluster_name & child2_name %in% distance_records$cluster_name) {
                            previous_P0_distance <- max(dplyr::filter(distance_records,
                                                                      cluster_name == child1_name |
                                                                        cluster_name == child2_name)$min_root_distance, na.rm = TRUE)
                            if (P0_distance > (previous_P0_distance*distance_awareness)) {
                              distance_conflict <- TRUE
                            }
                          }
                        }
                        if (distance_conflict == TRUE) {
                          # Add result to matrix
                          result_matrix[child1_name, child2_name] <- "split"
                          result_matrix[child2_name, child1_name] <- "split"
                          # Add result to comparison records
                          current_comparison <- matrix(rep(NA, ncol(comparison_records)), nrow = 1, ncol = ncol(comparison_records))
                          colnames(current_comparison) <- colnames(comparison_records)
                          current_comparison <- data.frame(current_comparison)
                          current_comparison$comparison <- paste0(child1_name, " vs. ", child2_name)
                          current_comparison$cluster1_size <- length(child1_cells)
                          current_comparison$cluster2_size <- length(child2_cells)
                          current_comparison$root_distance <- P0_distance
                          current_comparison$subtree_distance <- NA
                          current_comparison$connectivity <- adjacent
                          current_comparison$decision <- "split: distance"
                          current_comparison$time <- round(difftime(Sys.time(), comparison_start_time, units = "secs"), 2)
                          comparison_records <- rbind(comparison_records, current_comparison)
                        } else {
                          comparison_output <- .runPermutationTest(cluster1_name = child1_name,
                                                                   cluster1_cells = child1_cells,
                                                                   cluster1_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                                batches[child1_cells],
                                                                                                NULL),
                                                                   cluster2_name = child2_name,
                                                                   cluster2_cells = child2_cells,
                                                                   cluster2_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                                batches[child2_cells],
                                                                                                NULL),
                                                                   alpha = ifelse(p_adjust != "none", adjusted_alpha, alpha),
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
                                                                   max_n_batch = max_n_batch,
                                                                   input_matrix = input_matrix,
                                                                   nn_matrix = nn_matrix,
                                                                   comparison_records = comparison_records,
                                                                   feature_importance_records = feature_importance_records,
                                                                   P0_distance = P0_distance,
                                                                   P_i_distance = NA,
                                                                   adjacent = adjacent,
                                                                   comparison_start_time = comparison_start_time,
                                                                   n_cores = n_cores,
                                                                   random_seed = random_seed)
                          # Add result to matrix
                          result_matrix[child1_name, child2_name] <- comparison_output[["result"]]
                          result_matrix[child2_name, child1_name] <- comparison_output[["result"]]
                          # If split, add distances to distance_records
                          if (comparison_output[["result"]] == "split" & distance_awareness < Inf) {
                            rownames(distance_records) <- distance_records$cluster_name
                            distance_records <- .addDistance(cluster1_name = child1_name,
                                                             cluster2_name = child2_name,
                                                             P0_distance = P0_distance,
                                                             P_i_distance = NA,
                                                             max_p = comparison_output[["max_p"]],
                                                             distance_records = distance_records)
                          }
                          # Update records
                          comparison_records <- comparison_output[["comparison_records"]]
                          feature_importance_records <- comparison_output[["feature_importance_records"]]
                        }
                      }
                    }
                  }
                  if (lvl >= 0) {
                    tick_amount <- (1/(n_child_clusters - child1))*0.9*(1/(n_child_clusters-1))*0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
                    pb$tick(tick_amount)
                    if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                                   (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
                      hour_start_time <- Sys.time()
                      pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                        " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                                        round(difftime(Sys.time(), start_time, units = "min"), 2),
                                        " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                    }
                    percent_done <- percent_done + tick_amount
                  } else {
                    if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
                      hour_start_time <- Sys.time()
                      pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                        " : Running additional comparisons, ",
                                        round(difftime(Sys.time(), start_time, units = "min"), 2),
                                        " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                    }
                  }
                }
              }
              # Progress
              if (lvl >= 0) {
                tick_amount <- 0.1*(1/(n_child_clusters-1))*0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
                pb$tick(tick_amount)
                if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                               (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
                  hour_start_time <- Sys.time()
                  pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                    " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                                    round(difftime(Sys.time(), start_time, units = "min"), 2),
                                    " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                }
                percent_done <- percent_done + tick_amount
              } else {
                if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
                  hour_start_time <- Sys.time()
                  pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                    " : Running additional comparisons, ",
                                    round(difftime(Sys.time(), start_time, units = "min"), 2),
                                    " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                }
              }
            } else {
              # All split
              result_matrix[child1_name, ] <- "split"
              result_matrix[, child1_name] <- "split"
              result_matrix[child1_name, child1_name] <- NA
              # Progress
              if (lvl >= 0) {
                tick_amount <- 0.1*(1/(n_child_clusters-1))*0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
                pb$tick(tick_amount)
                if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                               (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
                  hour_start_time <- Sys.time()
                  pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                    " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                                    round(difftime(Sys.time(), start_time, units = "min"), 2),
                                    " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                }
                percent_done <- percent_done + tick_amount
              } else {
                if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
                  hour_start_time <- Sys.time()
                  pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                    " : Running additional comparisons, ",
                                    round(difftime(Sys.time(), start_time, units = "min"), 2),
                                    " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                }
              }
            }
          }

          print("Identifying merges")

          # Identify which clusters will merge & update child IDs
          # If all clusters will merge
          if (sum(result_matrix == "merge", na.rm = TRUE) == n_child_clusters*(n_child_clusters - 1)) {
            # Update child IDs
            child_IDs[parent_inds] <- parent_IDs[parent_inds]
            new_clusters <- c(new_clusters, unique(parent_IDs[parent_inds]))
          } else {
            # For each child
            for (child in 1:n_child_clusters) {
              child_name <- unique_child_IDs[child]
              # Count comparisons
              n_merge <- sum(result_matrix[child_name,] == "merge", na.rm = TRUE)
              # If there is more than one cluster that is slated to merge
              if (n_merge > 1) {
                # Identify the set of clusters this cluster is slated to merge with
                partner_clusters <- colnames(result_matrix[, colnames(result_matrix) != child_name])[result_matrix[child_name, colnames(result_matrix) != child_name] == "merge"]
                # Check whether these partner clusters are also slated to merge with each other
                for (partner1 in 1:(length(partner_clusters)-1)) {
                  partner1_name <- partner_clusters[partner1]
                  for (partner2 in (partner1 + 1):length(partner_clusters)) {
                    partner2_name <- partner_clusters[partner2]
                    # If not, assess whether to bridge the clusters or not
                    if (result_matrix[partner1_name, partner2_name] == "split") {
                      # Get cell IDs belonging to each cluster
                      child_cells <- cell_IDs[which(child_IDs == child_name)]
                      partner1_cells <- cell_IDs[which(child_IDs == partner1_name)]
                      partner2_cells <- cell_IDs[which(child_IDs == partner2_name)]
                      # Compare child + partner1 vs. partner2
                      child_partner1_name <- paste0(child_name, ".", partner1_name)
                      child_partner1_cells <- c(child_cells, partner1_cells)
                      # If there is only 1 cell in partner 2, automatically merge
                      if (length(partner2_cells) == 1) {
                        comparison1_output <- list("result" = "merge")
                        comparison1_proceed <- FALSE
                      } else {
                        # Check whether clusters have been previously compared
                        previous_comparison <- .checkComparisonRecords(cluster1_name = child_partner1_name,
                                                                       cluster1_cells = child_partner1_cells,
                                                                       cluster2_name = partner2_name,
                                                                       cluster2_cells = partner2_cells,
                                                                       comparison_records = comparison_records,
                                                                       type = "bridge")
                        # If they have, fill in result matrix
                        if (previous_comparison[["previously_compared"]] == TRUE) {
                          comparison1_output <- list("result" = previous_comparison[["result"]])
                          comparison1_proceed <- FALSE
                        } else {
                          comparison1_proceed <- TRUE
                        }
                      }
                      # Compare child + partner1 vs. partner2
                      child_partner2_name <- paste0(child_name, ".", partner2_name)
                      child_partner2_cells <- c(child_cells, partner2_cells)
                      # If there is only 1 cell in partner 1, automatically merge
                      if (length(partner1_cells) == 1) {
                        comparison2_output <- list("result" = "merge")
                        comparison2_proceed <- FALSE
                      } else {
                        # Check whether clusters have been previously compared
                        previous_comparison <- .checkComparisonRecords(cluster1_name = child_partner2_name,
                                                                       cluster1_cells = child_partner2_cells,
                                                                       cluster2_name = partner1_name,
                                                                       cluster2_cells = partner1_cells,
                                                                       comparison_records = comparison_records,
                                                                       type = "bridge")
                        # If they have, fill in result matrix
                        if (previous_comparison[["previously_compared"]] == TRUE) {
                          comparison2_output <- list("result" = previous_comparison[["result"]])
                          comparison2_proceed <- FALSE
                        } else {
                          comparison2_proceed <- TRUE
                        }
                      }
                      if (comparison1_proceed == TRUE | comparison2_proceed == TRUE) {
                        if (comparison1_proceed == TRUE) {
                          comparison_start_time <- Sys.time()
                          child_partner1_comparison_names <- c(paste0(child_name, " vs. ", partner1_name),
                                                               paste0(partner1_name, " vs. ", child_name))
                          # Get distance between clusters
                          distance_check <- .checkDistance(object = object,
                                                           key = key,
                                                           cluster1_name = child_partner1_name,
                                                           cluster1_cells = child_partner1_cells,
                                                           cluster2_name = partner2_name,
                                                           cluster2_cells = partner2_cells,
                                                           distance_awareness = distance_awareness,
                                                           use_input_matrix = "P0",
                                                           reduction = NULL,
                                                           distance_records = NULL)
                          # Run comparison 1
                          comparison1_output <- .runPermutationTest(cluster1_name = child_partner1_name,
                                                                    cluster1_cells = child_partner1_cells,
                                                                    cluster1_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                                 batches[child_partner1_cells],
                                                                                                 NULL),
                                                                    cluster2_name = partner2_name,
                                                                    cluster2_cells = partner2_cells,
                                                                    cluster2_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                                 batches[partner2_cells],
                                                                                                 NULL),
                                                                    alpha = ifelse(p_adjust != "none", adjusted_alpha, alpha),
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
                                                                    max_n_batch = max_n_batch,
                                                                    input_matrix = input_matrix,
                                                                    nn_matrix = nn_matrix,
                                                                    comparison_records = comparison_records,
                                                                    feature_importance_records = feature_importance_records,
                                                                    P0_distance = distance_check[["P0_distance"]],
                                                                    P_i_distance = distance_check[["P_i_distance"]],
                                                                    comparison_start_time = comparison_start_time,
                                                                    n_cores = n_cores,
                                                                    random_seed = random_seed)
                          # If split, add distances to distance_records
                          if (comparison1_output[["result"]] == "split" & methods::is(distance_awareness, "numeric")) {
                            rownames(distance_records) <- distance_records$cluster_name
                            distance_records <- .addDistance(cluster1_name = child_partner1_name,
                                                             cluster2_name = partner2_name,
                                                             P0_distance = distance_check[["P0_distance"]],
                                                             P_i_distance = distance_check[["P_i_distance"]],
                                                             max_p = comparison1_output[["max_p"]],
                                                             distance_records = distance_records)
                          }
                          # Update records
                          comparison_records <- comparison1_output[["comparison_records"]]
                          feature_importance_records <- comparison1_output[["feature_importance_records"]]
                          comparison1 <- comparison_records %>% dplyr::filter(comparison == paste0(child_partner1_name, " vs. ", partner2_name))
                        }
                        if (comparison2_proceed == TRUE) {
                          comparison_start_time <- Sys.time()
                          child_partner2_comparison_names <- c(paste0(child_name, " vs. ", partner2_name),
                                                               paste0(partner2_name, " vs. ", child_name))
                          # Get distance between clusters
                          distance_check <- .checkDistance(object = object,
                                                           key = key,
                                                           cluster1_name = child_partner2_name,
                                                           cluster1_cells = child_partner2_cells,
                                                           cluster2_name = partner1_name,
                                                           cluster2_cells = partner1_cells,
                                                           distance_awareness = distance_awareness,
                                                           use_input_matrix = "P0",
                                                           reduction = NULL,
                                                           distance_records = NULL)
                          # Run comparison 2
                          comparison2_output <- .runPermutationTest(cluster1_name = child_partner2_name,
                                                                    cluster1_cells = child_partner2_cells,
                                                                    cluster1_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                                 batches[child_partner2_cells],
                                                                                                 NULL),
                                                                    cluster2_name = partner1_name,
                                                                    cluster2_cells = partner1_cells,
                                                                    cluster2_cell_batches = `if`(batch_correction_method == "Harmony",
                                                                                                 batches[partner1_cells],
                                                                                                 NULL),
                                                                    alpha = ifelse(p_adjust != "none", adjusted_alpha, alpha),
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
                                                                    max_n_batch = max_n_batch,
                                                                    input_matrix = input_matrix,
                                                                    nn_matrix = nn_matrix,
                                                                    comparison_records = comparison_records,
                                                                    feature_importance_records = feature_importance_records,
                                                                    P0_distance = distance_check[["P0_distance"]],
                                                                    P_i_distance = distance_check[["P_i_distance"]],
                                                                    comparison_start_time = comparison_start_time,
                                                                    n_cores = n_cores,
                                                                    random_seed = random_seed)
                          # If split, add distances to distance_records
                          if (comparison2_output[["result"]] == "split" & methods::is(distance_awareness, "numeric")) {
                            rownames(distance_records) <- distance_records$cluster_name
                            distance_records <- .addDistance(cluster1_name = child_partner2_name,
                                                             cluster2_name = partner1_name,
                                                             P0_distance = distance_check[["P0_distance"]],
                                                             P_i_distance = distance_check[["P_i_distance"]],
                                                             max_p = comparison2_output[["max_p"]],
                                                             distance_records = distance_records)
                          }
                          # Update records
                          comparison_records <- comparison2_output[["comparison_records"]]
                          feature_importance_records <- comparison2_output[["feature_importance_records"]]
                          comparison2 <- comparison_records %>% dplyr::filter(comparison == paste0(child_partner2_name, " vs. ", partner1_name))
                        }
                      }
                      # Update results matrix
                      if (comparison1_output[["result"]] == "merge" & comparison2_output[["result"]] == "split") {
                        result_matrix[child_name, partner1_name] <- "split"
                        result_matrix[partner1_name, child_name] <- "split"
                      } else if (comparison1_output[["result"]] == "split" & comparison2_output[["result"]] == "merge") {
                        result_matrix[child_name, partner2_name] <- "split"
                        result_matrix[partner2_name, child_name] <- "split"
                      } else if (comparison1_output[["result"]] == "split" & comparison2_output[["result"]] == "split") {
                        # If they both result in a "split", merge child cluster with the partner w/ the lowest distance
                        if (methods::is(distance_awareness, "numeric")) {
                          child_partner1_distance <- dplyr::filter(comparison_records,
                                                                   comparison %in% child_partner1_comparison_names)$root_distance
                          child_partner2_distance <- dplyr::filter(comparison_records,
                                                                   comparison %in% child_partner2_comparison_names)$root_distance
                          # If distance has not been calculated
                          if (length(child_partner1_distance) == 0) {
                            distance_check <- .checkDistance(object = object,
                                                             key = key,
                                                             cluster1_name = child_name,
                                                             cluster1_cells = child_cells,
                                                             cluster2_name = partner1_name,
                                                             cluster2_cells = partner1_cells,
                                                             distance_awareness = distance_awareness,
                                                             use_input_matrix = "P0",
                                                             reduction = NULL,
                                                             distance_records = NULL)
                            child_partner1_distance <- distance_check[["P0_distance"]]
                          }
                          if (length(child_partner2_distance) == 0) {
                            distance_check <- .checkDistance(object = object,
                                                             key = key,
                                                             cluster1_name = child_name,
                                                             cluster1_cells = child_cells,
                                                             cluster2_name = partner2_name,
                                                             cluster2_cells = partner2_cells,
                                                             distance_awareness = distance_awareness,
                                                             use_input_matrix = "P0",
                                                             reduction = NULL,
                                                             distance_records = NULL)
                            child_partner2_distance <- distance_check[["P0_distance"]]
                          }
                          # Compare
                          if (child_partner1_distance > child_partner2_distance) {
                            result_matrix[child_name, partner1_name] <- "split"
                            result_matrix[partner1_name, child_name] <- "split"
                          } else {
                            result_matrix[child_name, partner2_name] <- "split"
                            result_matrix[partner2_name, child_name] <- "split"
                          }
                        } else {
                          # Alternately, use accuracy scores
                          child_partner1_mean_accuracy <- dplyr::filter(comparison_records,
                                                                        comparison %in% child_partner1_comparison_names)$mean_accuracy
                          child_partner2_mean_accuracy <- dplyr::filter(comparison_records,
                                                                        comparison %in% child_partner2_comparison_names)$mean_accuracy
                          if (child_partner1_mean_accuracy > child_partner2_mean_accuracy) {
                            result_matrix[child_name, partner1_name] <- "split"
                            result_matrix[partner1_name, child_name] <- "split"
                          } else {
                            result_matrix[child_name, partner2_name] <- "split"
                            result_matrix[partner2_name, child_name] <- "split"
                          }
                        }
                      } # Else if both result in "merge", do nothing, allow bridge
                    }
                  }
                }
              }
            }
            # Loop through child clusters again to find merge groups
            merge_group_list <- vector(mode = "list", length = n_child_clusters)
            # Change this if all children merge
            all_merge <- FALSE
            # For each child
            for (child in 1:n_child_clusters) {
              child_name <- unique_child_IDs[child]
              # Stop if it becomes established that all children will merge
              if (all_merge == FALSE) {
                # count cells
                n_merge <- sum(result_matrix[child_name,] == "merge", na.rm = TRUE)
                n_splits <- sum(result_matrix[child_name,] == "split", na.rm = TRUE)
                # if all merge
                if (n_merge == (n_child_clusters-1)) {
                  all_merge <- TRUE
                } else if (n_splits == (n_child_clusters-1)) { # If all splits
                  merge_group_list[[child]] <- child_name
                } else {
                  merge_names <- rownames(result_matrix)[result_matrix[child_name,] == "merge"]
                  merge_names <- merge_names[!is.na(merge_names)]
                  merge_names <- c(child_name, merge_names)
                  # Check if this group overlaps with any that are already identified
                  if (child > 1) {
                    for (i in 1:(child-1)) {
                      overlap <- length(intersect(merge_names, merge_group_list[[i]]))
                      if (overlap > 0) {
                        # include all members & delete old group
                        merge_names <- unique(c(merge_group_list[[i]], merge_names))
                        merge_group_list[i] <- list(NULL)
                      }
                    }
                  }
                  merge_group_list[[child]] <- merge_names
                }
              }
            }
            # Convert merge groups to new cluster names
            new_labels_list <- .getNewLabels(merge_groups = merge_group_list,
                                             level = lvl,
                                             compiled_labels = compiled_cluster_labels)
            compiled_cluster_labels <- new_labels_list[["compiled_cluster_labels"]]
            merge_group_labels <- new_labels_list[["merge_group_labels"]]

            # Update child_IDs
            if (all_merge == TRUE) {
              child_IDs[parent_inds] <- parent_IDs[parent_inds]
              new_clusters <- c(new_clusters, unique(parent_IDs[parent_inds]))
            } else {
              # Make key
              merge_group_key <- data.frame(old = unique_child_IDs,
                                            new = unique_child_IDs)
              rownames(merge_group_key) <- merge_group_key$old
              for (m_g in 1:length(merge_group_list)) {
                if (!is.null(merge_group_list[[m_g]])) {
                  merge_group_key[merge_group_list[[m_g]],"new"] <- merge_group_labels[[m_g]]
                }
              }

              new_cluster_labels <- child_IDs[parent_inds]
              for (child in 1:n_child_clusters) {
                new_cluster_labels[new_cluster_labels == unique_child_IDs[child]] <- merge_group_key[unique_child_IDs[child], "new"][1]
              }
              new_clusters <- c(new_clusters, unique(new_cluster_labels)[!(unique(new_cluster_labels) %in% unique(child_IDs[parent_inds]))])
              child_IDs[parent_inds] <- new_cluster_labels
            }
          }
        }
      } else {
        # Progress
        if (lvl >= 0) {
          tick_amount <- 0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
          pb$tick(tick_amount)
          if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                         (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
            hour_start_time <- Sys.time()
            pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                              " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                              round(difftime(Sys.time(), start_time, units = "min"), 2),
                              " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
          }
          percent_done <- percent_done + tick_amount
        } else {
          if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
            hour_start_time <- Sys.time()
            pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                              " : Running additional comparisons, ",
                              round(difftime(Sys.time(), start_time, units = "min"), 2),
                              " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
          }
        }
      }
      # Progress
      if (lvl >= 0) {
        tick_amount <- 0.1*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
        pb$tick(tick_amount)
        if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                       (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                            round(difftime(Sys.time(), start_time, units = "min"), 2),
                            " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
        }
        percent_done <- percent_done + tick_amount
      } else {
        if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : Running additional comparisons, ",
                            round(difftime(Sys.time(), start_time, units = "min"), 2),
                            " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
        }
      }
    }

    # If there are new clusters, edit list of permitted comparisons
    if (length(new_clusters) > 0 & dplyr::n_distinct(child_IDs) > 1) {
      permitted_comparisons <- .getPermittedComparisons(object = object,
                                                        key = key,
                                                        clusters = child_IDs,
                                                        cell_IDs = cell_IDs,
                                                        nn_matrices = list("P0" = nn_matrix),
                                                        min_connections = min_connections,
                                                        collect_all_metrics = collect_all_metrics,
                                                        distance_awareness = distance_awareness,
                                                        new_clusters = new_clusters,
                                                        permitted_comparisons = permitted_comparisons,
                                                        n_cores = n_cores)
      new_clusters <- c()
    }

    # Update parent IDs
    if (lvl > 1) {
      # Move up clustering tree
      parent_IDs <- cluster_tree[, lvl-1]
    } else {
      # Continue until all clusters have been compared
      parent_IDs <- rep(paste0("P0_L", paste(rep(0, (((-1)*lvl) + 2)), collapse = ""), "_1"), length(cell_IDs))
    }
    # Record stepwise changes
    stepwise_child_IDs_df <- data.frame(stepwise_cluster_IDs = child_IDs)
    if (lvl >= 0) {
      colnames(stepwise_child_IDs_df) <- paste0("stepwise_cluster_ID_", alpha, "_L", lvl)
    } else {
      colnames(stepwise_child_IDs_df) <- paste0("stepwise_cluster_ID_", alpha, "_L", paste(rep(0, abs(lvl) + 1), collapse = ""))
    }
    stepwise_cluster_IDs <- cbind(stepwise_cluster_IDs, stepwise_child_IDs_df)
    new_clusters <- c()
    # Progress
    if (lvl >= 0) {
      tick_amount <- 0.1*level_weights[paste0("L", lvl)]
      pb$tick(tick_amount)
      if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                     (difftime(Sys.time(), hour_start_time, units = "hours") >= 1))) {
        hour_start_time <- Sys.time()
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                          " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                          round(difftime(Sys.time(), start_time, units = "min"), 2),
                          " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
      }
      percent_done <- percent_done + tick_amount
    } else {
      if (verbose & (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5)) {
        hour_start_time <- Sys.time()
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                          " : Running additional comparisons, ",
                          round(difftime(Sys.time(), start_time, units = "min"), 2),
                          " min elapsed. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
      }
    }

    # Check for completion if beyond root of clustering tree
    if (lvl < 1) {
      if ((dplyr::n_distinct(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)]) ==
           dplyr::n_distinct(paste(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)], stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)-1])))) {
        complete <- TRUE
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Completed: all clusters compared."))
        pb$tick(5)
      } else if (dplyr::n_distinct(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)]) == 1) {
        complete <- TRUE
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Completed: only one cluster remaining."))
        pb$tick(5)
      } else {
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Additional comparisons necessary."))
      }
    }
    # Increment level
    lvl <- lvl - 1
  }

  # Finalize cluster IDs
  n_final_clusters <- dplyr::n_distinct(child_IDs)
  final_clusters <- data.frame(CellID = cell_IDs,
                               Record_cluster_label = child_IDs)
  # Convert to numeric cluster labels based on cluster size
  cluster_key <- final_clusters %>%
    dplyr::group_by(Record_cluster_label) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::arrange(-n)
  cluster_key$CHOIR_ID <- 1:n_final_clusters
  final_clusters <- merge(final_clusters, cluster_key, by = "Record_cluster_label", all.x = TRUE) %>%
    dplyr::select(CellID, CHOIR_ID, Record_cluster_label)
  colnames(final_clusters)[2] <- paste0("CHOIR_clusters_", alpha)

  # Rearrange back to original order
  rownames(final_clusters) <- final_clusters$CellID
  final_clusters <- final_clusters[cell_IDs, ]

  # Add to object
  object <- .storeData(object, key, "final_clusters",
                       data.frame(CHOIR_IDs = final_clusters[, paste0("CHOIR_clusters_", alpha)]),
                       paste0("CHOIR_clusters_", alpha))
  object <- .storeData(object, key, "clusters", final_clusters, paste0("CHOIR_clusters_", alpha))
  object <- .storeData(object, key, "clusters", stepwise_cluster_IDs, paste0("stepwise_clusters_", alpha))

  # -------------------------------------------------------------------------
  # Report results & warnings
  # -------------------------------------------------------------------------

  if (verbose) {
    message("")
    if (max_repeat_errors > 0 & sum(comparison_records$decision == "split: repeat error") > 0) {
      message(" - Repeatedly misassigned cells affected the result of ", sum(comparison_records$decision == "split: repeat error"),
              " comparisons. Set 'max_repeat_errors' parameter to 0 if you would like to disable this setting.")
    }
    if (methods::is(distance_awareness, "numeric") & sum(comparison_records$decision == "split: distance") > 0) {
      message(" - In ", sum(comparison_records$decision == "split: distance"),
              " comparisons, clusters were split due to cluster distance threshold.")
    }
    if (min_connections > 0 & sum(comparison_records$decision == "split: min connections") > 0) {
      message(" - In ", sum(comparison_records$decision == "split: min connections"),
              " comparisons, clusters were split due to the minimum number of nearest neighbor connections.")
    }
    if (sum(comparison_records$decision == "merge: min accuracy") > 0) {
      message(" - In ", sum(comparison_records$decision == "merge: min accuracy"),
              " comparisons, clusters were merged due to the minimum accuracy threshold.")
    }
    if (p_adjust != "none" & sum(comparison_records$decision == "merge: adjustment") > 0) {
      message(" - In ", sum(comparison_records$decision == "merge: adjustment"),
              " comparisons, clusters were merged after correction for multiple comparisons.")
      message(" - Final adjusted significance threshold = ", ifelse(adjusted_alpha > 0.0001, round(adjusted_alpha, 5), signif(adjusted_alpha, digits=5)), ".")
    }
    if (sum(comparison_records$decision == "merge: batch-dependent") > 0) {
      message(" - In ", sum(comparison_records$decision == "merge: batch-dependent"),
              " comparisons, clusters were merged due to batch-dependence.")
    }
    message("\n", format(Sys.time(), "%Y-%m-%d %X"), " : Identified ", n_final_clusters, " clusters.")
  }

  # Record parameters used and add to original object
  parameter_list <- list("alpha" = alpha,
                         "p_adjust" = p_adjust,
                         "adjusted_alpha" = adjusted_alpha,
                         "feature_set" = feature_set,
                         "exclude_features" = exclude_features,
                         "n_iterations" = n_iterations,
                         "n_trees" = n_trees,
                         "use_variance" = use_variance,
                         "min_accuracy" = min_accuracy,
                         "min_connections" = min_connections,
                         "max_repeat_errors" = max_repeat_errors,
                         "distance_awareness" = distance_awareness,
                         "collect_all_metrics" = collect_all_metrics,
                         "sample_max" = sample_max,
                         "downsampling_rate" = downsampling_rate,
                         "normalization_method" = normalization_method,
                         "batch_correction_method" = batch_correction_method,
                         "batch_labels" = batch_labels,
                         "max_n_batch" = max_n_batch,
                         "input_matrix_provided" = input_matrix_provided,
                         "nn_matrix_provided" = nn_matrix_provided,
                         "reduction_provided" = reduction_provided,
                         "random_seed" = random_seed)

  object <- .storeData(object, key, "parameters", parameter_list, "combineTrees_parameters")

  # Add records to object
  object <- .storeData(object, key, "records", comparison_records, "comparison_records")
  if (collect_all_metrics == TRUE) {
    object <- .storeData(object, key, "records", feature_importance_records, "feature_importance_records")
  }

  # Return object with new additions
  return(object)
}
