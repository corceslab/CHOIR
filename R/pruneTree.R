#' Prune clustering tree using random forest classifiers
#'
#' To identify a final set of clusters, this function will move iteratively from
#' the bottom up to prune the provided hierarchical clustering tree using a
#' framework of random forest classifiers and permutation tests.
#'
#' If \code{CHOIR::buildTree()} was run prior to this function, most parameters
#' will be retrieved from the object. Alternately, parameter values can be
#' supplied. For multi-modal data, optionally supply parameter inputs as
#' vectors/lists that sequentially specify the value for each modality.
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{ArchRProject}. For multi-omic data, we recommend using
#' \code{ArchRProject} objects.
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
#' @param distance_approx A Boolean value indicating whether or not to use
#' approximate distance calculations. Defaults to \code{TRUE}, which will use
#' centroid-based distances. Setting distance approximation to \code{FALSE} will
#' substantially increase the computational time and memory required,
#' particularly for large datasets. Using approximated distances (\code{TRUE})
#' rather than absolute distances (\code{FALSE}) is unlikely to have a
#' meaningful effect on the distance thresholds imposed by CHOIR.
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
#' \code{FALSE}. Setting this parameter to \code{FALSE} or increasing the input
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
#' @param cluster_params A list of additional parameters to be passed to
#' \code{Seurat} function \code{FindClusters} for clustering at each level of
#' the tree. By default, when the \code{Seurat::FindClusters} parameter
#' \code{group.singletons} is set to \code{TRUE}, singletons are grouped into
#' the nearest cluster. Alternately, if \code{group.singletons} is set to
#' \code{FALSE}, CHOIR will relabel clusters such that each singleton
#' constitutes its own cluster.
#' @param use_assay For \code{Seurat} or \code{SingleCellExperiment} objects, a
#' character string or vector indicating the assay(s) to use in the provided
#' object. The default value, \code{NULL}, will choose the current active assay
#' for \code{Seurat} objects and the \code{logcounts} assay for
#' \code{SingleCellExperiment} objects.
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
#' @param cluster_tree An optional dataframe containing the cluster IDs of each
#' cell across the levels of a hierarchical clustering tree. Default = \code{NULL}
#' will use the hierarchical clustering tree generation by function
#' \code{buildTree()}.
#' @param input_matrix An optional matrix containing the feature x cell data
#' provided by the user, on which to train the random forest classifiers. By
#' default, this parameter is set to \code{NULL}, and CHOIR will look for the
#' feature x cell matri(ces) indicated by function \code{buildTree}.
#' @param nn_matrix An optional matrix containing the nearest neighbor adjacency
#' of the cells, provided by the user. By default, this parameter is set to
#' \code{NULL}, and CHOIR will look for the adjacency matri(ces) generated by
#' function \code{buildTree}.
#' @param snn_matrix An optional matrix containing the shared nearest neighbor
#' adjacency of the cells, provided by the user. By default, this parameter is
#' set to \code{NULL}, and CHOIR will look for the adjacency matri(ces)
#' generated by function \code{buildTree}.
#' @param dist_matrix An optional distance matrix of cell to cell distances
#' (based on dimensionality reduction cell embeddings), provided by the user. By
#' default, this parameter is set to \code{NULL}, and CHOIR will look for the
#' distance matri(ces) generated by function \code{buildTree}.
#' @param reduction An optional matrix of dimensionality reduction cell
#' embeddings provided by the user for subsequent clustering steps. By default,
#' this parameter is set to \code{NULL}, and CHOIR will look for the
#' dimensionality reductions generated by function \code{buildTree}.
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
#' @return Returns the object with the following added data stored under the
#' provided key: \describe{
#'   \item{clusters}{Final clusters and stepwise cluster results for each
#'   progressive pruning step}
#'   \item{parameters}{Record of parameter values used}
#'   \item{records}{Metadata for all recorded permutation test comparisons and
#'   feature importance scores from all comparisons}
#'   }
#'
#' @export
#'
pruneTree <- function(object,
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
                      distance_approx = NULL,
                      distance_awareness = 2,
                      collect_all_metrics = FALSE,
                      sample_max = NULL,
                      downsampling_rate = NULL,
                      min_reads = NULL,
                      normalization_method = NULL,
                      batch_correction_method = NULL,
                      batch_labels = NULL,
                      max_n_batch = NULL,
                      cluster_params = NULL,
                      use_assay = NULL,
                      countsplit = NULL,
                      countsplit_suffix = NULL,
                      cluster_tree = NULL,
                      input_matrix = NULL,
                      nn_matrix = NULL,
                      snn_matrix = NULL,
                      dist_matrix = NULL,
                      reduction = NULL,
                      n_cores = NULL,
                      random_seed = NULL,
                      verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Retrieve/check parameter input validity
  # ---------------------------------------------------------------------------

  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 1/2) Checking inputs and preparing object..")

  .validInput(object, "object", "pruneTree")
  .validInput(key, "key", list("pruneTree", object))
  .validInput(distance_awareness, "distance_awareness")
  .validInput(collect_all_metrics, "collect_all_metrics")
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

  # Retrieve parameter values from buildTree
  buildTree_parameters <- .retrieveData(object,
                                        key,
                                        "parameters",
                                        "buildTree_parameters")
  default_parameters <- list("alpha" = 0.05,
                             "p_adjust" = "bonferroni",
                             "feature_set" = "var",
                             "exclude_features" = NULL,
                             "n_iterations" = 100,
                             "n_trees" = 50,
                             "use_variance" = TRUE,
                             "min_accuracy" = 0.5,
                             "min_connections" = 1,
                             "max_repeat_errors" = 20,
                             "distance_approx" = TRUE,
                             "sample_max" = Inf,
                             "downsampling_rate" = "auto",
                             "min_reads" = NULL,
                             "normalization_method" = "none",
                             "batch_correction_method" = "none",
                             "batch_labels" = NULL,
                             "max_n_batch" = Inf,
                             "cluster_params" = list(algorithm = 1,
                                                     group.singletons = TRUE),
                             "use_assay"  = NULL,
                             "countsplit" = FALSE,
                             "countsplit_suffix" = NULL,
                             "countsplit_text" = "",
                             "random_seed" = 1)

  # For any parameters set to NULL, use parameters from buildTree or defaults
  alpha <- .retrieveParam(alpha, "alpha", buildTree_parameters, default_parameters)
  p_adjust <- .retrieveParam(p_adjust, "p_adjust", buildTree_parameters, default_parameters)
  feature_set <- .retrieveParam(feature_set, "feature_set", buildTree_parameters, default_parameters)
  exclude_features <- .retrieveParam(exclude_features, "exclude_features", buildTree_parameters, default_parameters)
  n_iterations <- .retrieveParam(n_iterations, "n_iterations", buildTree_parameters, default_parameters)
  n_trees <- .retrieveParam(n_trees, "n_trees", buildTree_parameters, default_parameters)
  use_variance <- .retrieveParam(use_variance, "use_variance", buildTree_parameters, default_parameters)
  min_accuracy <- .retrieveParam(min_accuracy, "min_accuracy", buildTree_parameters, default_parameters)
  min_connections <- .retrieveParam(min_connections, "min_connections", buildTree_parameters, default_parameters)
  max_repeat_errors <- .retrieveParam(max_repeat_errors, "max_repeat_errors", buildTree_parameters, default_parameters)
  distance_approx <- .retrieveParam(distance_approx, "distance_approx", buildTree_parameters, default_parameters)
  sample_max <- .retrieveParam(sample_max, "sample_max", buildTree_parameters, default_parameters)
  downsampling_rate <- .retrieveParam(downsampling_rate, "downsampling_rate", buildTree_parameters, default_parameters)
  min_reads <- .retrieveParam(min_reads, "min_reads", buildTree_parameters, default_parameters)
  normalization_method <- .retrieveParam(normalization_method, "normalization_method", buildTree_parameters, default_parameters)
  batch_correction_method <- .retrieveParam(batch_correction_method, "batch_correction_method", buildTree_parameters, default_parameters)
  batch_labels <- .retrieveParam(batch_labels, "batch_labels", buildTree_parameters, default_parameters)
  max_n_batch <- .retrieveParam(max_n_batch, "max_n_batch", buildTree_parameters, default_parameters)
  cluster_params <- .retrieveParam(cluster_params, "cluster_params", buildTree_parameters, default_parameters)
  use_assay <- .retrieveParam(use_assay, "use_assay", buildTree_parameters, default_parameters)
  countsplit <- .retrieveParam(countsplit, "countsplit", buildTree_parameters, default_parameters)
  countsplit_suffix <- .retrieveParam(countsplit_suffix, "countsplit_suffix", buildTree_parameters, default_parameters)
  countsplit_text <- .retrieveParam(NULL, "countsplit_text", buildTree_parameters, default_parameters)
  random_seed <- .retrieveParam(random_seed, "random_seed", buildTree_parameters, default_parameters)
  # Verify parameter validity
  .validInput(alpha, "alpha")
  .validInput(p_adjust, "p_adjust")
  .validInput(feature_set, "feature_set", list("pruneTree", normalization_method))
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
  .validInput(cluster_params, "cluster_params")
  .validInput(use_assay, "use_assay", object)
  .validInput(countsplit, "countsplit")
  .validInput(countsplit_suffix, "countsplit_suffix", countsplit)
  .validInput(random_seed, "random_seed")

  # Add additional parameters if not provided
  if (!any(names(cluster_params) == "verbose")) {
    cluster_params$verbose <- FALSE
  }
  if (!any(names(cluster_params) == "algorithm")) {
    cluster_params$algorithm <- 1
  }
  if (!any(names(cluster_params) == "group.singletons")) {
    cluster_params$group.singletons <- TRUE
  }

  # Extract cell IDs
  cell_IDs <- .getCellIDs(object, use_assay)

  .validInput(cluster_tree, "cluster_tree", object)
  .validInput(input_matrix, "input_matrix", cell_IDs)
  .validInput(nn_matrix, "nn_matrix", cell_IDs)
  .validInput(snn_matrix, "snn_matrix", cell_IDs)
  .validInput(dist_matrix, "dist_matrix", object)
  .validInput(reduction, "reduction", list("pruneTree", object))

  # User supplied data should not be mixed with data from buildTree()
  if (methods::is(distance_awareness, "numeric") & distance_approx == FALSE) {
    if (any(c(!is.null(cluster_tree), !is.null(input_matrix), !is.null(nn_matrix), !is.null(snn_matrix), !is.null(dist_matrix))) &
        !all(c(!is.null(cluster_tree), !is.null(input_matrix), !is.null(nn_matrix), !is.null(snn_matrix), !is.null(dist_matrix)))) {
      warning("If user-supplied data is provided for any of 'cluster_tree', 'input_matrix', 'nn_matrix', 'snn_matrix', or 'dist_matrix', it should be provided for all five.")
    }
  } else if (methods::is(distance_awareness, "numeric") & distance_approx == TRUE) {
    if (any(c(!is.null(cluster_tree), !is.null(input_matrix), !is.null(nn_matrix), !is.null(snn_matrix), !is.null(reduction))) &
        !all(c(!is.null(cluster_tree), !is.null(input_matrix), !is.null(nn_matrix), !is.null(snn_matrix), !is.null(reduction)))) {
      warning("If user-supplied data is provided for any of 'cluster_tree', 'input_matrix', 'nn_matrix', 'snn_matrix', or 'reduction', it should be provided for all five.")
    }
  } else {
    if (any(c(!is.null(cluster_tree), !is.null(input_matrix), !is.null(nn_matrix), !is.null(snn_matrix))) &
        !all(c(!is.null(cluster_tree), !is.null(input_matrix), !is.null(nn_matrix), !is.null(snn_matrix)))) {
      warning("If user-supplied data is provided for any of 'cluster_tree', 'input_matrix', 'nn_matrix', or 'snn_matrix', it should be provided for all four.")
    }
  }
  # User supplied data cannot be used with countsplitting
  if (!is.null(input_matrix) & countsplit == TRUE) {
    warning("Countsplitting is not currently compatible with user-supplied data for 'input_matrix'. Parameter 'count_split' set to FALSE.")
    countsplit <- FALSE
    countsplit_suffix <- NULL
  }

  # Retrieve subtree data
  if (!is.null(buildTree_parameters)) {
    # Subtree names
    subtree_names <- buildTree_parameters[["subtree_names"]]
    subtree_sizes <- buildTree_parameters[["subtree_sizes"]]
    subtree_names_filtered <- subtree_names[subtree_sizes > 3]
    n_subtrees <- length(subtree_names)
    n_subtrees_filtered <- length(subtree_names_filtered)
  } else {
    # If there are no records of parameters used in buildTree function
    warning("No record of parameters used in buildTree() function.")
    # Subtree names
    subtree_names <- "P0"
    subtree_sizes <- length(cell_IDs)
    subtree_names_filtered <- subtree_names[subtree_sizes > 3]
    n_subtrees <- length(subtree_names)
    n_subtrees_filtered <- length(subtree_names_filtered)
  }

  # ---------------------------------------------------------------------------
  # Prepare object
  # ---------------------------------------------------------------------------

  # Extract cluster tree
  if (is.null(cluster_tree)) {
    cluster_tree <- .retrieveData(object, key, "clusters", "full_tree")
    .validInput(cluster_tree, "cluster_tree", object)
    cluster_tree_provided <- FALSE
  } else {
    cluster_tree_provided <- TRUE
  }
  # Check whether provided clustering tree is ordered correctly & strictly hierarchical
  cluster_tree <- .checkHierarchy(cluster_tree)
  # Check whether provided clustering tree has appropriately named clusters
  # If not, rename them
  cluster_tree <- .checkClusterLabels(cluster_tree)
  rownames(cluster_tree) <- cell_IDs
  # Number of starting clusters
  n_starting_clusters <- dplyr::n_distinct(cluster_tree[, ncol(cluster_tree)])
  # Number of levels in cluster tree
  n_levels <- ncol(cluster_tree)
  # Will change based on multiple comparison correction, if applicable
  adjusted_alpha <- alpha

  # Extract input matrix/matrices
  if (!is.null(input_matrix)) {
    input_matrix_provided <- TRUE
    input_matrix <- .getMatrix(use_matrix = input_matrix,
                               exclude_features = exclude_features,
                               verbose = FALSE)
    input_matrix <- BiocGenerics::t(input_matrix)
    input_matrices <- list("P0" = input_matrix)
    n_input_matrices <- 1
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
    .validInput(distance_approx, "distance_approx", list(length(cell_IDs), object_type, n_modalities))
  } else if (!is.null(buildTree_parameters)) {
    input_matrix_provided <- FALSE
    # Parameters for extracting matrices
    use_assay <- buildTree_parameters[["use_assay"]]
    use_slot <- buildTree_parameters[["use_slot"]]
    ArchR_matrix <- buildTree_parameters[["ArchR_matrix"]]
    use_assay_build <- buildTree_parameters[["use_assay_build"]]
    use_assay_prune <- buildTree_parameters[["use_assay_prune"]]
    use_slot_build <- buildTree_parameters[["use_slot_build"]]
    use_slot_prune <- buildTree_parameters[["use_slot_prune"]]
    ArchR_matrix_build <- buildTree_parameters[["ArchR_matrix_build"]]
    ArchR_matrix_prune <- buildTree_parameters[["ArchR_matrix_prune"]]
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
    .validInput(distance_approx, "distance_approx", list(length(cell_IDs), object_type, n_modalities))

    # If reduction was not recalculated for buildTree, there will only be one input matrix
    if (buildTree_parameters[["subtree_reductions"]] == FALSE) {
      n_input_matrices <- 1
    } else {
      n_input_matrices <- n_subtrees_filtered
    }
    # Initialize list of feature names
    features <- c()
    # Initialize list of input matrices
    input_matrices <- vector(mode = "list", n_subtrees_filtered)

    # Extract data
    for (subtree in 1:n_input_matrices) {
      subtree_name <- subtree_names_filtered[subtree]
      use_cells_subtree <- .retrieveData(object, key, "cell_IDs", paste0(subtree_name, "_cell_IDs"))
      # Use all / var features
      if (feature_set == "all") {
        use_features_subtree <- NULL
      } else {
        use_features_subtree <- .retrieveData(object, key, "var_features", paste0(subtree_name, "_var_features"))
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
          if (length(use_features_subtree) > 1) {
            use_features_subtree_m <- use_features_subtree[[m]]
          } else {
            use_features_subtree_m <- use_features_subtree
          }
          input_matrix_list[[m]] <- .getMatrix(object = object,
                                               use_assay = use_assay_prune_m,
                                               use_slot = use_slot_prune_m,
                                               ArchR_matrix = ArchR_matrix_prune_m,
                                               use_features = use_features_subtree_m,
                                               exclude_features = exclude_features,
                                               use_cells = use_cells_subtree,
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
        rm(use_features_subtree_m)
      } else {
        input_matrix <- .getMatrix(object = object,
                                   use_assay = use_assay_prune,
                                   use_slot = use_slot_prune,
                                   ArchR_matrix = ArchR_matrix_prune,
                                   use_features = use_features_subtree,
                                   exclude_features = exclude_features,
                                   use_cells = use_cells_subtree,
                                   verbose = FALSE)
        if (normalization_method == "SCTransform") {
          input_matrix <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(input_matrix),
                                                               return.only.var.genes = FALSE,
                                                               seed.use = random_seed,
                                                               verbose = FALSE)@assays$SCT@scale.data)
        }
        input_matrix <- BiocGenerics::t(input_matrix)
      }
      # Add to list of input_matrices
      input_matrices[[subtree]] <- input_matrix
      # Add features to vector
      features <- unique(c(features, colnames(input_matrix)))
    }
    names(input_matrices) <- subtree_names_filtered
    # Clean up
    rm(use_features_subtree)
    rm(input_matrix)
  } else {
    stop("No 'input_matrix' supplied. Please supply valid input!")
  }

  # Extract nearest neighbor matrix/matrices
  if (!is.null(nn_matrix)) {
    nn_matrix_provided <- TRUE
    nn_matrices <- list("P0" = nn_matrix)
    # Clean up
    rm(nn_matrix)
  } else if (!is.null(buildTree_parameters)) {
    nn_matrix_provided <- FALSE
    # Initialize list of nearest neighbor adjacency matrices
    nn_matrices <- vector(mode = "list", n_subtrees_filtered)
    # For each subtree
    for (subtree in 1:n_subtrees_filtered) {
      subtree_name <- subtree_names_filtered[subtree]
      nn_matrices[[subtree]] <- .retrieveData(object, key, "graph", paste0(subtree_name, "_graph_nn"))
    }
    names(nn_matrices) <- subtree_names_filtered
  } else {
    stop("No nearest neighbor adjacency matrix provided.")
  }

  # Check on shared nearest neighbor matrix/matrices
  if (!is.null(snn_matrix)) {
    snn_matrix_provided <- TRUE
  } else if (!is.null(buildTree_parameters)) {
    snn_matrix_provided <- FALSE
  } else {
    stop("No shared nearest neighbor adjacency matrix provided.")
  }

  # Reduction
  if (!is.null(reduction)) {
    reduction_provided <- TRUE
  } else {
    reduction_provided <- FALSE
  }

  # Distance matrix
  if (!is.null(dist_matrix)) {
    dist_matrix_provided <- TRUE
  } else {
    dist_matrix_provided <- FALSE
  }

  # If 'Harmony' batch correction was used, extract batch labels
  if (batch_correction_method == "Harmony") {
    batches <- .retrieveData(object = object, key = key, type = "cell_metadata", name = batch_labels)
    if ("Rle" %in% methods::is(batches)) {
      batches <- methods::as(batches, "character")
    }
    names(batches) <- cell_IDs
  }

  # Set downsampling rate
  if (downsampling_rate == "auto") {
    downsampling_rate <- min(1, (1/2)^(log10(length(cell_IDs)/5000)))
    if (batch_correction_method == "none") {
      downsampling_rate <- downsampling_rate*0.5
    }
  }

  # Provided input
  provided_input <- c(
    if (cluster_tree_provided) "cluster_tree",
    if (input_matrix_provided) "input_matrix",
    if (nn_matrix_provided) "nn_matrix",
    if (snn_matrix_provided) "snn_matrix",
    if (dist_matrix_provided) "dist_matrix",
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
                       "\n - Distance approximation: ", distance_approx,
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
                       "\n - Clustering parameters provided: ", `if`(length(cluster_params) == 0, "No",
                                                                     paste0("\n     - ", paste0(paste0(names(cluster_params), ": ",
                                                                                                       cluster_params),
                                                                                                collapse = "\n     - "))),
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

  # Record of distances
  if (methods::is(distance_awareness, "numeric")) {
    distance_records <- data.frame(cluster_name = NULL,
                                   min_root_distance = NULL,
                                   min_subtree_distance = NULL,
                                   max_pval = NULL)
  } else {
    distance_records <- NULL
  }

  # Record of stepwise cluster IDs
  stepwise_cluster_IDs <- data.frame(CellID = cell_IDs)

  # Record underclustering checks
  checked_for_underclustering <- c()
  results_of_underclustering_check <- c()
  underclustering_buffer <- FALSE

  # Set of all cluster labels in original tree
  compiled_cluster_labels <- unlist(apply(cluster_tree, 2, unique))
  names(compiled_cluster_labels) <- NULL

  # -------------------------------------------------------------------------
  # Iterate through each level of the cluster tree (bottom-up)
  # -------------------------------------------------------------------------

  # Track progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/2) Iterating through clustering tree..")

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

      if (n_child_clusters > 1) {
        # Create matrix for comparison results
        result_matrix <- matrix(NA, n_child_clusters, n_child_clusters)
        colnames(result_matrix) <- unique_child_IDs
        rownames(result_matrix) <- unique_child_IDs
        # For each child cluster
        for (child1 in 1:(n_child_clusters-1)) {
          # Name of child cluster 1
          child1_name <- unique_child_IDs[child1]
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
              comparison_start_time <- Sys.time()
              child2_name <- unique_child_IDs[child2]
              # Get names of cells belonging to child1 cluster
              child2_cells <- cell_IDs[which(child_IDs == child2_name)]
              # If there is only 1 cell in child2, automatically merge
              if (length(child2_cells) == 1) {
                result_matrix[child1_name, child2_name] <- "merge"
                result_matrix[child2_name, child1_name] <- "merge"
              } else {
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
                  # Determine which input & nn matrices to use
                  roots <- unlist(stringr::str_extract_all(c(child1_name, child2_name), "P\\d*"))
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
                  # Check whether distance between current clusters is greater than previous comparisons
                  distance_check <- .checkDistance(object = object,
                                                   key = key,
                                                   cluster1_name = child1_name,
                                                   cluster1_cells = child1_cells,
                                                   cluster2_name = child2_name,
                                                   cluster2_cells = child2_cells,
                                                   distance_awareness = distance_awareness,
                                                   distance_approx = distance_approx,
                                                   use_input_matrix = use_input_matrix,
                                                   dist_matrix = dist_matrix,
                                                   reduction = reduction,
                                                   distance_records = distance_records)
                  if (distance_check[["distance_conflict"]] == TRUE) {
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
                    current_comparison$root_distance <- distance_check[["P0_distance"]]
                    current_comparison$subtree_distance <- distance_check[["P_i_distance"]]
                    current_comparison$decision <- "split: distance"
                    current_comparison$time <- round(difftime(Sys.time(), comparison_start_time, units = "secs"), 2)
                    comparison_records <- rbind(comparison_records, current_comparison)
                  } else if (distance_check[["distance_conflict"]] == FALSE) {
                    # Run comparison
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
                                                             input_matrix = input_matrices[[use_input_matrix]],
                                                             nn_matrix = nn_matrices[[use_nn_matrix]],
                                                             comparison_records = comparison_records,
                                                             feature_importance_records = feature_importance_records,
                                                             P0_distance = distance_check[["P0_distance"]],
                                                             P_i_distance = distance_check[["P_i_distance"]],
                                                             comparison_start_time = comparison_start_time,
                                                             n_cores = n_cores,
                                                             random_seed = random_seed)
                    # Add result to matrix
                    result_matrix[child1_name, child2_name] <- comparison_output[["result"]]
                    result_matrix[child2_name, child1_name] <- comparison_output[["result"]]
                    # If split, add distances to distance_records
                    if (comparison_output[["result"]] == "split" & methods::is(distance_awareness, "numeric")) {
                      distance_records <- .addDistance(cluster1_name = child1_name,
                                                       cluster2_name = child2_name,
                                                       P0_distance = distance_check[["P0_distance"]],
                                                       P_i_distance = distance_check[["P_i_distance"]],
                                                       max_p = comparison_output[["max_p"]],
                                                       distance_records = distance_records)
                    }
                    # Update records
                    comparison_records <- comparison_output[["comparison_records"]]
                    feature_importance_records <- comparison_output[["feature_importance_records"]]
                  }
                }
              }
              if (lvl >= 0) {
                tick_amount <- (1/(n_child_clusters - child1))*0.9*(1/(n_child_clusters-1))*0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
                pb$tick(tick_amount)
                if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                               (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
                  hour_start_time <- Sys.time()
                  pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                    " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                                    round(difftime(Sys.time(), start_time, units = "min"), 2),
                                    " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
                }
                percent_done <- percent_done + tick_amount
              }
            }
          }
          if (lvl >= 0) {
            tick_amount <- 0.1*(1/(n_child_clusters-1))*0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
            pb$tick(tick_amount)
            if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                           (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
              hour_start_time <- Sys.time()
              pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                                " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                                round(difftime(Sys.time(), start_time, units = "min"), 2),
                                " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
            }
            percent_done <- percent_done + tick_amount
          }
        }
        # Identify which clusters will merge & update child IDs
        # If all clusters will merge
        if (sum(result_matrix == "merge", na.rm = TRUE) == n_child_clusters*(n_child_clusters - 1)) {
          # Update child IDs
          child_IDs[parent_inds] <- parent_IDs[parent_inds]
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
                      # Determine which input & nn matrices to use
                      roots <- unlist(stringr::str_extract_all(c(child_name, partner1_name, partner2_name), "P\\d*"))
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
                                                         distance_approx = distance_approx,
                                                         use_input_matrix = use_input_matrix,
                                                         dist_matrix = dist_matrix,
                                                         reduction = reduction,
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
                                                                  input_matrix = input_matrices[[use_input_matrix]],
                                                                  nn_matrix = nn_matrices[[use_nn_matrix]],
                                                                  comparison_records = comparison_records,
                                                                  feature_importance_records = feature_importance_records,
                                                                  P0_distance = distance_check[["P0_distance"]],
                                                                  P_i_distance = distance_check[["P_i_distance"]],
                                                                  comparison_start_time = comparison_start_time,
                                                                  n_cores = n_cores,
                                                                  random_seed = random_seed)
                        # If split, add distances to distance_records
                        if (comparison1_output[["result"]] == "split" & methods::is(distance_awareness, "numeric")) {
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
                                                         distance_approx = distance_approx,
                                                         use_input_matrix = use_input_matrix,
                                                         dist_matrix = dist_matrix,
                                                         reduction = reduction,
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
                                                                  input_matrix = input_matrices[[use_input_matrix]],
                                                                  nn_matrix = nn_matrices[[use_nn_matrix]],
                                                                  comparison_records = comparison_records,
                                                                  feature_importance_records = feature_importance_records,
                                                                  P0_distance = distance_check[["P0_distance"]],
                                                                  P_i_distance = distance_check[["P_i_distance"]],
                                                                  comparison_start_time = comparison_start_time,
                                                                  n_cores = n_cores,
                                                                  random_seed = random_seed)
                        # If split, add distances to distance_records
                        if (comparison2_output[["result"]] == "split" & methods::is(distance_awareness, "numeric")) {
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
                        if (use_input_matrix == "P0") {
                          child_partner1_distance <- dplyr::filter(comparison_records,
                                                                   comparison %in% child_partner1_comparison_names)$root_distance
                          child_partner2_distance <- dplyr::filter(comparison_records,
                                                                   comparison %in% child_partner2_comparison_names)$root_distance
                        } else {
                          child_partner1_distance <- dplyr::filter(comparison_records,
                                                                   comparison %in% child_partner1_comparison_names)$subtree_distance
                          child_partner2_distance <- dplyr::filter(comparison_records,
                                                                   comparison %in% child_partner2_comparison_names)$subtree_distance
                        }
                        # If distance has not been calculated
                        if (length(child_partner1_distance) == 0) {
                          distance_check <- .checkDistance(object = object,
                                                           key = key,
                                                           cluster1_name = child_name,
                                                           cluster1_cells = child_cells,
                                                           cluster2_name = partner1_name,
                                                           cluster2_cells = partner1_cells,
                                                           distance_awareness = distance_awareness,
                                                           distance_approx = distance_approx,
                                                           use_input_matrix = use_input_matrix,
                                                           dist_matrix = dist_matrix,
                                                           reduction = reduction,
                                                           distance_records = NULL)
                          if (use_input_matrix == "P0") {
                            child_partner1_distance <- distance_check[["P0_distance"]]
                          } else {
                            child_partner1_distance <- distance_check[["P_i_distance"]]
                          }
                        }
                        if (length(child_partner2_distance) == 0) {
                          distance_check <- .checkDistance(object = object,
                                                           key = key,
                                                           cluster1_name = child_name,
                                                           cluster1_cells = child_cells,
                                                           cluster2_name = partner2_name,
                                                           cluster2_cells = partner2_cells,
                                                           distance_awareness = distance_awareness,
                                                           distance_approx = distance_approx,
                                                           use_input_matrix = use_input_matrix,
                                                           dist_matrix = dist_matrix,
                                                           reduction = reduction,
                                                           distance_records = NULL)
                          if (use_input_matrix == "P0") {
                            child_partner2_distance <- distance_check[["P0_distance"]]
                          } else {
                            child_partner2_distance <- distance_check[["P_i_distance"]]
                          }
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
            child_IDs[parent_inds] <- new_cluster_labels
          }
        }
      } else {
        # No need to update child IDs if there is only one child cluster
        # (Cluster will retain same ID unless it merges, to retain ability to use records)
        # Progress
        if (lvl >= 0) {
          tick_amount <- 0.9*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
          pb$tick(tick_amount)
          if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                         (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
            hour_start_time <- Sys.time()
            pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                              " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                              round(difftime(Sys.time(), start_time, units = "min"), 2),
                              " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
          }
          percent_done <- percent_done + tick_amount
        }
      }
      # Progress
      if (lvl >= 0) {
        tick_amount <- 0.1*(1/length(unique_parent_IDs))*(0.9*level_weights[paste0("L", lvl)])
        pb$tick(tick_amount)
        if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                       (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                            round(difftime(Sys.time(), start_time, units = "min"), 2),
                            " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
        }
        percent_done <- percent_done + tick_amount
      }
    }
    # Check multiple comparison adjustment
    # Starting after at least 10 total comparisons & at least 3 "split" calls
    # Or if we're at the top of the tree with at least 1 "split" call
    if (p_adjust != "none" &
        ((lvl <= 1 & sum(comparison_records$decision %in% c("split", "split: repeat error")) >= 1) |
         (length(comparison_records$percentile_accuracy[!is.na(comparison_records$percentile_accuracy)]) >=
          ifelse(p_adjust == "fdr", 100, 50) &
          sum(comparison_records$decision %in% c("split", "split: repeat error")) >=
          ifelse(p_adjust == "fdr", 5, 1)))) {

      # Set a new alpha threshold
      p_accuracy <- comparison_records$percentile_accuracy[!is.na(comparison_records$percentile_accuracy)]
      p_variance <- comparison_records$percentile_variance[!is.na(comparison_records$percentile_variance)]
      if (p_adjust == "fdr") {
        # Based on FDR
        p_accuracy_adj <- p_accuracy[stats::p.adjust(p_accuracy, method = "fdr") < alpha]
        if (length(p_accuracy_adj) > 0) {
          p_accuracy_max <- (max(which(sort(p_accuracy) == max(p_accuracy_adj)))/length(p_accuracy))*alpha
        } else {
          p_accuracy_max <- adjusted_alpha
        }
        p_variance_adj <- p_variance[stats::p.adjust(p_variance, method = "fdr") < alpha]
        if (length(p_variance_adj) > 0) {
          p_variance_max <- (max(which(sort(p_variance) == max(p_variance_adj)))/length(p_variance))*alpha
        } else {
          p_variance_max <- adjusted_alpha
        }
        adjusted_alpha <- min(adjusted_alpha, max(p_variance_max, p_variance_max, na.rm = TRUE), na.rm = TRUE)
      } else if (p_adjust == "bonferroni") {
        adjusted_alpha <- alpha/length(p_accuracy)
      }
      # Remove distance records in which the comparison p-value was above the new adjusted alpha value
      distance_records <- dplyr::filter(distance_records, max_pval < adjusted_alpha)

      # Check for cases where p-values are now below alpha threshold
      # Where both clusters are still among the current cluster IDs
      if (use_variance == TRUE) {
        correction_check <- comparison_records %>%
          dplyr::mutate(correct = ifelse(decision == "split" &
                                           (percentile_accuracy >= adjusted_alpha |
                                              percentile_variance >= adjusted_alpha), TRUE,
                                         ifelse(decision == "split: repeat error" &
                                                  (percentile_modified_accuracy >= adjusted_alpha |
                                                     percentile_modified_variance >= adjusted_alpha), TRUE, FALSE))) %>%
          dplyr::filter(correct == TRUE) %>%
          tidyr::separate_wider_delim(comparison, delim = " vs. ", names = c("cluster1", "cluster2")) %>%
          dplyr::filter(cluster1 %in% unique(child_IDs),
                        cluster2 %in% unique(child_IDs))
      } else {
        correction_check <- comparison_records %>%
          dplyr::mutate(correct = ifelse(decision == "split" &
                                           percentile_accuracy >= adjusted_alpha, TRUE,
                                         ifelse(decision == "split: repeat error" &
                                                  percentile_modified_accuracy >= adjusted_alpha, TRUE, FALSE))) %>%
          dplyr::filter(correct == TRUE) %>%
          tidyr::separate_wider_delim(comparison, delim = " vs. ", names = c("cluster1", "cluster2")) %>%
          dplyr::filter(cluster1 %in% unique(child_IDs),
                        cluster2 %in% unique(child_IDs))
      }


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
            if (methods::is(distance_awareness, "numeric")) {
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
                                           level = lvl,
                                           compiled_labels = compiled_cluster_labels)
          compiled_cluster_labels <- new_labels_list[["compiled_cluster_labels"]]
          merged_label <- new_labels_list[["merge_group_labels"]][[1]]

          child_IDs[child_IDs %in% c(merge_pair$cluster1[1], merge_pair$cluster2[1])] <- merged_label
          # Update records
          comparison_records$decision[which(comparison_records$comparison ==
                                              paste0(merge_pair$cluster1[1], " vs. ",
                                                     merge_pair$cluster2[1]))] <- "merge: adjustment"
          n_current_clusters <- n_current_clusters - 1
        }
      }
    }
    # Update parent IDs
    if (lvl > 1) {
      # Move up clustering tree
      parent_IDs <- cluster_tree[, lvl-1]
    } else {
      # Continue until all clusters have been compared
      parent_IDs <- rep(paste0("P0_L", paste(rep(0, (((-1)*lvl) + 1)), collapse = ""), "_1"), length(cell_IDs))
    }
    # Record stepwise changes
    stepwise_child_IDs_df <- data.frame(stepwise_cluster_IDs = child_IDs)
    if (lvl >= 0) {
      colnames(stepwise_child_IDs_df) <- paste0("stepwise_cluster_ID_", alpha, "_L", lvl)
    } else {
      colnames(stepwise_child_IDs_df) <- paste0("stepwise_cluster_ID_", alpha, "_L", paste(rep(0, abs(lvl) + 1), collapse = ""))
    }
    stepwise_cluster_IDs <- cbind(stepwise_cluster_IDs, stepwise_child_IDs_df)
    # Progress
    if (lvl >= 0) {
      tick_amount <- 0.1*level_weights[paste0("L", lvl)]
      pb$tick(tick_amount)
      if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                     (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
        hour_start_time <- Sys.time()
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                          " : ", round((percent_done + tick_amount)), "% (", n_levels - lvl, "/", n_levels ," levels) in ",
                          round(difftime(Sys.time(), start_time, units = "min"), 2),
                          " min. ", dplyr::n_distinct(child_IDs), " clusters remaining."))
      }
      percent_done <- percent_done + tick_amount
    }
    # Check for completion if beyond root of clustering tree
    if (lvl < 1) {
      if (underclustering_buffer == FALSE &
          (dplyr::n_distinct(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)]) ==
           dplyr::n_distinct(paste(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)], stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)-1])))) {
        # Check if any clusters from bottom of tree are still present
        check_df <- data.frame(original = cluster_tree[, n_levels],
                            current = stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)])
        check <- check_df %>%
          dplyr::group_by(current, original) %>%
          dplyr::filter(!(current %in% checked_for_underclustering)) %>%
          dplyr::filter(!(original %in% checked_for_underclustering)) %>%
          dplyr::summarise(cells = dplyr::n()) %>%
          dplyr::group_by(current) %>%
          dplyr::summarise(n = dplyr::n()) %>%
          dplyr::filter(n == 1)

        clusters_to_check <- unique(child_IDs)[unique(child_IDs) %in% results_of_underclustering_check]
        clusters_to_check <- clusters_to_check[!(clusters_to_check %in% checked_for_underclustering)]
        clusters_to_check <- c(clusters_to_check, unique(check$current))
        checked_for_underclustering <- unique(c(checked_for_underclustering,
                                         clusters_to_check,
                                         unique(dplyr::filter(check_df, current %in% check$current)$original)))

        if (length(clusters_to_check) > 0) {
          # Replace parent IDs with current child IDs
          parent_IDs <- child_IDs
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Checking for underclustering in ",
                            length(clusters_to_check), " clusters."))
          for (u_clust in 1:length(clusters_to_check)) {
            current_cell_inds <- which(child_IDs == clusters_to_check[u_clust])
            current_cell_IDs <- cell_IDs[current_cell_inds]
            use_input_matrix <- unlist(stringr::str_extract_all(clusters_to_check[u_clust], "P\\d*"))
            if ("subtree_reductions" %in% names(buildTree_parameters)) {
              subtree_reductions <- buildTree_parameters[["subtree_reductions"]]
              if (subtree_reductions == FALSE) {
                use_input_matrix <- "P0"
              }
            }
            # Build subtree
            subtree_list <- .getTree(snn_matrix = `if`(!is.null(snn_matrix),
                                                       snn_matrix[current_cell_IDs, current_cell_IDs],
                                                       .retrieveData(object, key, "graph", paste0(use_input_matrix, "_graph_snn"))[current_cell_IDs, current_cell_IDs]),
                                     nn_matrix = nn_matrices[[use_input_matrix]][current_cell_IDs, current_cell_IDs],
                                     dist_matrix = `if`(distance_approx == FALSE,
                                                        `if`(dist_matrix_provided == TRUE, dist_matrix[current_cell_IDs, current_cell_IDs],
                                                             .retrieveData(object, key, "reduction", paste0(use_input_matrix, "_reduction_dist"))[current_cell_IDs, ]), NULL),
                                     reduction = `if`(distance_approx == TRUE,
                                                      `if`(reduction_provided == TRUE, reduction[current_cell_IDs, ],
                                                           .retrieveData(object, key, "reduction", paste0(use_input_matrix, "_reduction"))[current_cell_IDs, ]), NULL),
                                     input_matrix = input_matrices[[use_input_matrix]][current_cell_IDs, ],
                                     distance_approx = distance_approx,
                                     tree_type = "subtree",
                                     cluster_params = cluster_params,
                                     alpha = ifelse(p_adjust != "none", adjusted_alpha, alpha),
                                     exclude_features = exclude_features,
                                     n_iterations = n_iterations,
                                     n_trees = n_trees,
                                     use_variance = use_variance,
                                     min_accuracy = min_accuracy,
                                     min_connections = min_connections,
                                     max_repeat_errors = max_repeat_errors,
                                     sample_max = sample_max,
                                     downsampling_rate = downsampling_rate,
                                     batch_correction_method = batch_correction_method,
                                     batches = `if`(batch_correction_method == "Harmony",
                                                    batches[current_cell_IDs],
                                                    NULL),
                                     max_n_batch = max_n_batch,
                                     tree_records = data.frame(tree_type = NULL,
                                                               tree_name = NULL,
                                                               num_cells = NULL,
                                                               resolution = NULL,
                                                               num_clusters = NULL,
                                                               silhouette = NULL,
                                                               neighbors_distance = NULL,
                                                               neighbors_mean_accuracy = NULL,
                                                               neighbors_var_accuracy = NULL,
                                                               neighbors_percentile_accuracy = NULL,
                                                               neighbors_percentile_variance = NULL,
                                                               neighbors_decision = NULL,
                                                               stop_branching_reason = NULL),
                                     n_cores = n_cores,
                                     random_seed = random_seed)
            # Add bottom subcluster IDs to child IDs
            # Convert to new cluster names
            subcluster_labels <- as.numeric(as.factor(subtree_list[["cluster_tree"]][, ncol(subtree_list[["cluster_tree"]])]))
            if (dplyr::n_distinct(subcluster_labels) > 1) {
              underclustering_buffer <- TRUE
              unique_cluster_labels <- unique(c(cluster_tree[, n_levels], child_IDs))
              unique_cluster_labels <- unique_cluster_labels[grepl(paste0(use_input_matrix, "_"), unique_cluster_labels)]
              max_P_label <- max(as.numeric(sub("P\\d*_L\\d*_", "", unique_cluster_labels)))
              new_cluster_IDs <- paste0(use_input_matrix, "_", "L",
                                        ifelse(lvl >= 0, lvl, paste(rep(0, abs(lvl) + 1), collapse = "")),
                                        "_", subcluster_labels + max_P_label)
              child_IDs[current_cell_inds] <- new_cluster_IDs
              results_of_underclustering_check <- c(results_of_underclustering_check, unique(new_cluster_IDs))
            }
          }
        } else {
          complete <- TRUE
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Completed: all clusters compared."))
          pb$tick(5)
        }
      } else if (dplyr::n_distinct(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)]) == 1) {
        complete <- TRUE
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Completed: only one cluster remaining."))
        pb$tick(5)
      } else {
        pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"), " : Additional comparisons necessary. ",
                          dplyr::n_distinct(stepwise_cluster_IDs[, ncol(stepwise_cluster_IDs)]), " clusters remaining."))
        underclustering_buffer <- FALSE
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

  # Pull out accuracy scores for comparisons between final clusters
  final_cluster_mean_accuracies <- matrix(NA, nrow = n_final_clusters, ncol = n_final_clusters)
  final_cluster_distances <- matrix(NA, nrow = n_final_clusters, ncol = n_final_clusters)
  for (i in 1:(n_final_clusters-1)) {
    cluster_i_name <- cluster_key[cluster_key$CHOIR_ID == i, "Record_cluster_label"]
    for (j in (i+1):n_final_clusters) {
      cluster_j_name <- cluster_key[cluster_key$CHOIR_ID == j, "Record_cluster_label"]

      comparison_ij <- dplyr::filter(comparison_records,
                                  comparison == paste0(cluster_i_name, " vs. ", cluster_j_name) |
                                    comparison == paste0(cluster_j_name, " vs. ", cluster_i_name))
      if (nrow(comparison_ij) == 1) {
        final_cluster_mean_accuracies[i, j] <- comparison_ij$mean_accuracy
        final_cluster_mean_accuracies[j, i] <- comparison_ij$mean_accuracy
        final_cluster_distances[i, j] <- comparison_ij$root_distance
        final_cluster_distances[j, i] <- comparison_ij$root_distance
      }
    }
  }

  # Add to object
  object <- .storeData(object, key, "final_clusters",
                       data.frame(CHOIR_IDs = final_clusters[, paste0("CHOIR_clusters_", alpha)]),
                       paste0("CHOIR_clusters_", alpha))
  object <- .storeData(object, key, "clusters", final_clusters, paste0("CHOIR_clusters_", alpha))
  object <- .storeData(object, key, "clusters", stepwise_cluster_IDs, paste0("stepwise_clusters_", alpha))
  object <- .storeData(object, key, "records", final_cluster_mean_accuracies, paste0("CHOIR_clusters_", alpha, "_accuracies"))
  object <- .storeData(object, key, "records", final_cluster_distances, paste0("CHOIR_clusters_", alpha, "_distances"))

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
  parameter_list <- list("subtree_names" = names(input_matrices),
                         "alpha" = alpha,
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
                         "distance_approx" = distance_approx,
                         "distance_awareness" = distance_awareness,
                         "collect_all_metrics" = collect_all_metrics,
                         "sample_max" = sample_max,
                         "downsampling_rate" = downsampling_rate,
                         "min_reads" = min_reads,
                         "normalization_method" = normalization_method,
                         "batch_correction_method" = batch_correction_method,
                         "batch_labels" = batch_labels,
                         "max_n_batch" = max_n_batch,
                         "countsplit" = countsplit,
                         "countsplit_suffix" = countsplit_suffix,
                         "countsplit_text" = countsplit_text,
                         "cluster_tree_provided" = cluster_tree_provided,
                         "input_matrix_provided" = input_matrix_provided,
                         "nn_matrix_provided" = nn_matrix_provided,
                         "snn_matrix_provided" = snn_matrix_provided,
                         "dist_matrix_provided" = dist_matrix_provided,
                         "reduction_provided" = reduction_provided,
                         "random_seed" = random_seed)

  object <- .storeData(object, key, "parameters", parameter_list, "pruneTree_parameters")

  # Add records to object
  object <- .storeData(object, key, "records", comparison_records, "comparison_records")
  if (collect_all_metrics == TRUE) {
    object <- .storeData(object, key, "records", feature_importance_records, "feature_importance_records")
  }

  # Return object with new additions
  return(object)
}
