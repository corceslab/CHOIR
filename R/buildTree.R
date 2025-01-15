#' Build full hierarchical clustering tree
#'
#' This function performs the first step of the CHOIR algorithm. It constructs
#' a hierarchical clustering tree starting from a single cluster encompassing
#' all cells. First, a root tree is constructed, from which subtrees are
#' subsequently generated. Each branch is subdivided until all cells are
#' demonstrably overclustered.
#'
#' For multi-modal data, optionally supply parameter inputs as vectors/lists
#' that sequentially specify the value for each modality.
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
#' @param max_clusters Indicates the extent to which the hierarchical clustering
#' tree will be expanded. Defaults to “auto”, which will expand the tree until
#' cases of underclustering have been eliminated in all branches. Alternatively,
#' supply a numerical value indicating the maximum number of clusters to which
#' to the tree should be expanded. Using the default value for this parameter is
#' highly recommended to avoid instances of underclustering. Setting a numerical
#' value in this parameter hampers the ability of CHOIR to ensure that
#' underclustering has not occurred.
#' @param min_cluster_depth A numerical value indicating the maximum cluster
#' size at the bottom of the clustering tree, prior to pruning branches.
#' Defaults to 2000. Increasing this parameter can cause a computational
#' bottleneck when generating the initial clustering tree for some large
#' datasets; therefore, the default value is recommended. However, changing this
#' value is unlikely to have meaningful effects on the final cluster results.
#' @param normalization_method A character string or vector indicating which
#' normalization method to use. In general, input data should be supplied to
#' CHOIR after normalization, except when the user wishes to use
#' \code{Seurat SCTransform} normalization. Permitted values are “none” or
#' “SCTransform”. Defaults to “none”. Because CHOIR has not been tested
#' thoroughly with \code{SCTransform} normalization, we do not recommend this
#' approach at this time. For multi-omic datasets, provide a vector with a value
#' corresponding to each provided value of \code{use_assay} or
#' \code{ArchR_matrix} in the same order.
#' @param subtree_reductions A Boolean value indicating whether to generate a
#' new dimensionality reduction and set of highly variable features for each
#' subtree. Defaults to \code{TRUE}, which enables CHOIR to compare similar
#' clusters using a more nuanced set of features. Setting this parameter to
#' \code{FALSE} may decrease computational time, but may result in
#' instances of underclustering.
#' @param reduction_method A character string or vector indicating which
#' dimensionality reduction method to use. Permitted values are “PCA” for
#' principal component analysis, “LSI” for latent semantic indexing, and
#' “IterativeLSI” for iterative latent semantic indexing. These three methods
#' implement the \code{Seurat} function \code{RunPCA}, the \code{Signac}
#' function \code{RunSVD}, and the \code{ArchR} function \code{addIterativeLSI},
#' respectively. The default value, \code{NULL}, will select a method based on
#' the input data type, specifically “IterativeLSI” for \code{ArchR} objects,
#' “LSI” for \code{Seurat} or \code{SingleCellExperiment} objects when parameter
#' \code{atac} is \code{TRUE}, and “PCA” in all other cases. For multi-omic
#' datasets, provide a vector with a value corresponding to each provided value
#' of \code{use_assay} or \code{ArchR_matrix} in the same order.
#' @param reduction_params A list of additional parameters to be passed to the
#' selected dimensionality reduction method. By default, CHOIR will use the
#' default parameter settings of the dimensionality reduction method indicated
#' by the input to parameter reduction_method. Input to this parameter is passed
#' to each downstream dimensionality reduction method and will overwrite or
#' augment those defaults. Altering the performance of the dimensionality
#' reduction in CHOIR will affect downstream clustering results, but not in ways
#' that are easily predictable.
#' @param n_var_features A numerical value indicating how many variable features
#' to identify. Defaults to 2000 features for most data inputs, or 25000
#' features for ATAC-seq data. Increasing the number of features may increase
#' the computational time and memory required. If the provided value is either
#' substantially higher or lower, instances of underclustering may occur. For
#' multi-omic datasets, provide a vector with a value corresponding to each
#' provided value of \code{use_assay} or \code{ArchR_matrix} in the same order.
#' @param batch_correction_method A character string indicating which batch
#' correction method to use. Permitted values are “Harmony” and “none”. Defaults
#' to “none”. Batch correction should only be used when the different batches
#' are not expected to also have unique cell types or cell states. Using batch
#' correction would ensure that clusters do not originate from a single batch,
#' thereby making the final cluster calls more conservative.
#' @param batch_correction_params A list of additional parameters to be passed
#' to the selected batch correction method for each iteration. Only applicable
#' when \code{batch_correction_method} is “Harmony”.
#' @param batch_labels A character string that, if applying batch correction,
#' specifies the name of the column in the input object metadata containing the
#' batch labels. Defaults to \code{NULL}.
#' @param neighbor_params A list of additional parameters to be passed to
#' \code{Seurat} function \code{FindNeighbors} (or, in the case of multi-modal
#' data for \code{Seurat} or \code{SingleCellExperiment} objects, \code{Seurat}
#' function \code{FindMultiModalNeighbors}).
#' @param cluster_params A list of additional parameters to be passed to
#' \code{Seurat} function \code{FindClusters} for clustering at each level of
#' the tree. By default, when the \code{Seurat::FindClusters} parameter
#' \code{group.singletons} is set to \code{TRUE}, CHOIR relabels clusters such
#' that each singleton constitutes its own cluster.
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
#' provide a vector with a value corresponding to each modality. When
#' "GeneScoreMatrix" is provided, the "GeneScoreMatrix" will be used as input
#' to the random forest classifiers, but the "TileMatrix" will be used for the
#' initial dimensionality reduction(s).
#' @param ArchR_depthcol For \code{ArchR} objects, a character string or vector
#' indicating which column to use for correlation with sequencing depth. The
#' default value, \code{NULL}, will use the “nFrags” column for ATAC-seq data or
#' the “Gex_nUMI” for RNA-seq data. For multi-omic datasets, provide a vector
#' with a value corresponding to each provided value of \code{ArchR_matrix} in
#' the same order.
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
#' @param reduction An optional matrix of dimensionality reduction cell
#' embeddings provided by the user for subsequent clustering steps. By default,
#' this parameter is set to \code{NULL}, and the dimensionality reduction(s)
#' will be calculated using the method specified by the \code{reduction_method}
#' parameter.
#' @param var_features An optional character vector of names of variable
#' features to be used for subsequent clustering steps. By default, this
#' parameter is set to \code{NULL}, and variable features will be calculated as
#' part of running CHOIR. Input to this parameter is required when a
#' dimensionality reduction is supplied to parameter \code{reduction}. For
#' multi-omic datasets, concatenate feature names for all modalities.
#' @param atac A Boolean value or vector indicating whether the provided data is
#' ATAC-seq data. For multi-omic datasets, provide a vector with a value
#' corresponding to each provided value of \code{use_assay} or
#' \code{ArchR_matrix} in the same order. Defaults to \code{FALSE}.
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
#'   \item{cell_IDs}{Cell IDs belonging to each subtree}
#'   \item{clusters}{Full hierarchical clustering tree}
#'   \item{graph}{All calculated nearest neighbor and shared nearest neighbor
#'   adjacency matrices}
#'   \item{parameters}{Record of parameter values used}
#'   \item{records}{Metadata for decision points during hierarchical tree
#'   construction}
#'   \item{reduction}{Cell embeddings for all calculated dimensionality
#'   reductions}
#'   \item{var_features}{Variable features for all calculated dimensionality
#'   reductions}
#'   }
#'
#' @export
#'
buildTree <- function(object,
                      key = "CHOIR",
                      alpha = 0.05,
                      p_adjust = "bonferroni",
                      feature_set = "var",
                      exclude_features = NULL,
                      n_iterations = 100,
                      n_trees = 50,
                      use_variance = TRUE,
                      min_accuracy = 0.5,
                      min_connections = 1,
                      max_repeat_errors = 20,
                      distance_approx = TRUE,
                      sample_max = Inf,
                      downsampling_rate = "auto",
                      min_reads = NULL,
                      max_clusters = "auto",
                      min_cluster_depth = 2000,
                      normalization_method = "none",
                      subtree_reductions = TRUE,
                      reduction_method = NULL,
                      reduction_params = list(),
                      n_var_features = NULL,
                      batch_correction_method = "none",
                      batch_correction_params = list(),
                      batch_labels = NULL,
                      neighbor_params = list(),
                      cluster_params = list(algorithm = 1,
                                            group.singletons = TRUE),
                      use_assay = NULL,
                      use_slot = NULL,
                      ArchR_matrix = NULL,
                      ArchR_depthcol = NULL,
                      countsplit = FALSE,
                      countsplit_suffix = NULL,
                      reduction = NULL,
                      var_features = NULL,
                      atac = FALSE,
                      n_cores = NULL,
                      random_seed = 1,
                      verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(max_clusters, "max_clusters")

  # Progress
  if (verbose  & max_clusters == "auto") message(format(Sys.time(), "%Y-%m-%d %X"),
                                                 " : (Step 1/7) Checking inputs and preparing object..")
  if (verbose  & max_clusters != "auto") message(format(Sys.time(), "%Y-%m-%d %X"),
                                                 " : (Step 1/5) Checking inputs and preparing object..")


  .validInput(object, "object", "buildTree")
  .validInput(key, "key", list("buildTree", object))
  .validInput(alpha, "alpha")
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
  .validInput(max_clusters, "max_clusters")
  .validInput(min_cluster_depth, "min_cluster_depth")
  .validInput(countsplit, "countsplit")
  .validInput(countsplit_suffix, "countsplit_suffix", countsplit)
  .validInput(use_assay, "use_assay", list(object, countsplit, countsplit_suffix))
  .validInput(use_slot, "use_slot", list(object, use_assay, countsplit, countsplit_suffix))
  .validInput(ArchR_matrix, "ArchR_matrix", list(object, countsplit, countsplit_suffix))

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
    } else if (methods::is(object, "SingleCellExperiment")) {
      object_type <- "SingleCellExperiment"
      .requirePackage("SingleCellExperiment", source = "bioc")
    }
  }

  # Get cell IDs
  cell_IDs <- .getCellIDs(object, use_assay = use_assay)

  .validInput(distance_approx, "distance_approx", list(length(cell_IDs), object, n_modalities))
  .validInput(ArchR_depthcol, "ArchR_depthcol", list(object, n_modalities))
  .validInput(reduction, "reduction", list("buildTree", object))
  .validInput(var_features, "var_features", reduction)
  .validInput(atac, "atac", n_modalities)
  .validInput(normalization_method, "normalization_method", list(object, n_modalities, use_assay))
  .validInput(feature_set, "feature_set", list("buildTree", normalization_method, atac))
  .validInput(subtree_reductions, "subtree_reductions", list(reduction, max_clusters))
  .validInput(reduction_method, "reduction_method", list(object, n_modalities))
  .validInput(reduction_params, "reduction_params", list(object, ArchR_matrix, use_assay, reduction_method))
  .validInput(n_var_features, "n_var_features", n_modalities)
  .validInput(batch_correction_method, "batch_correction_method", n_modalities)
  .validInput(batch_correction_params, "batch_correction_params", list(object, ArchR_matrix, use_assay, batch_correction_method))
  .validInput(batch_labels, "batch_labels", object)
  .validInput(neighbor_params, "neighbor_params", list(object, n_modalities))
  .validInput(cluster_params, "cluster_params")
  .validInput(n_cores, "n_cores")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # Add additional parameters if not provided
  if (!any(names(neighbor_params) == "verbose")) {
    neighbor_params$verbose <- FALSE
  }
  if (!any(names(cluster_params) == "verbose")) {
    cluster_params$verbose <- FALSE
  }
  if (!any(names(cluster_params) == "algorithm")) {
    cluster_params$algorithm <- 1
  }
  if (!any(names(cluster_params) == "group.singletons")) {
    cluster_params$group.singletons <- TRUE
  }

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }
  # Set downsampling rate
  if (downsampling_rate == "auto") {
    downsampling_rate <- min(1, (1/2)^(log10(length(cell_IDs)/5000)))
    if (batch_correction_method == "none") {
      downsampling_rate <- downsampling_rate*0.5
    }
  }
  # Random seed reproducibility
  if (n_cores > 1) {
    RNGkind("L'Ecuyer-CMRG")
  }

  # Check that required packages are loaded
  # Seurat needs to be loaded to use Seurat:::FindModalityWeights() and Seurat:::MultiModalNN()
  if (n_modalities >= 2 & !methods::is(object, "ArchRProject")) {
    .requirePackage("Seurat", source = "cran")
  }

  # If batch_correction_method is Harmony, make sure batch IDs is a character vector
  if (batch_correction_method == "Harmony") {
    batches <- as.character(.retrieveData(object = object, key = key, type = "cell_metadata", name = batch_labels))
    if (!("character" %in% methods::is(batches))) {
      warning("Metadata column '", batch_labels, "' converted to type character.")
    }
    batches <- as.character(batches)
    .storeData(object = object, key = key, type = "cell_metadata", name = batch_labels, input_data = batches)
    names(batches) <- cell_IDs
  } else {
    batches <- NULL
  }

  # ---------------------------------------------------------------------------
  # Set values for countsplitting
  # ---------------------------------------------------------------------------

  if (countsplit == TRUE) {
    if (is.null(countsplit_suffix)) {
      countsplit_suffix <- c("_1", "_2")
    } else {
      countsplit_suffix <- c("", "")
    }
  }
  # Set new values
  if (methods::is(object, "Seurat")) {
    # Seurat object
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
                              paste(ArchR_matrix_prune, collapse = ", "),
                              "\n - ArchR depth column",
                              ifelse(n_modalities == 1, ": ", "s: "),
                              paste(ArchR_depthcol, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Report object & parameter details
  # ---------------------------------------------------------------------------

  if (verbose) message("\nInput data:",
                       "\n - Object type: ", ifelse(object_type == "Seurat", paste0(object_type, " (", seurat_version, ")"), object_type),
                       `if`(!is.null(reduction) | !is.null(var_features), "\n - Provided inputs: ", ""),
                       `if`(!is.null(reduction), "reduction", ""),
                       `if`(!is.null(reduction) & !is.null(var_features), ", ", ""),
                       `if`(!is.null(reduction), "var_features", ""),
                       "\n - # of cells: ", length(cell_IDs),
                       "\n - # of batches: ", `if`(batch_correction_method == "none", 1, dplyr::n_distinct(batches)),
                       "\n - # of modalities: ", n_modalities,
                       "\n - ATAC data: ", paste(atac, collapse = ", "),
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
                       "\n - Maximum cells sampled: ", sample_max,
                       "\n - Downsampling rate: ", round(downsampling_rate, 4),
                       "\n - Minimum reads: ", `if`(is.null(min_reads),
                                                    paste0(">0 reads"), paste0(">1 read per ", min_reads, " cells")),
                       "\n - Maximum clusters: ", max_clusters,
                       "\n - Minimum cluster depth: ", min_cluster_depth,
                       "\n - Normalization method: ", normalization_method,
                       "\n - Subtree dimensionality reductions: ", subtree_reductions,
                       "\n - Dimensionality reduction method: ", `if`(is.null(reduction_method), "Default", reduction_method),
                       "\n - Dimensionality reduction parameters provided: ", `if`(length(reduction_params) == 0, "No",
                                                                                   paste0("\n     - ", paste0(paste0(names(reduction_params), ": ",
                                                                                                                     reduction_params),
                                                                                                              collapse = "\n     - "))),
                       "\n - # of variable features: ", `if`(!is.null(var_features), length(var_features),
                                                             `if`(is.null(n_var_features), "Default", n_var_features)),
                       "\n - Batch correction method: ", batch_correction_method,
                       "\n - Batch correction parameters provided: ", `if`(length(batch_correction_params) == 0, "No",
                                                                           paste0("\n     - ", paste0(paste0(names(batch_correction_params), ": ",
                                                                                                             batch_correction_params),
                                                                                                      collapse = "\n     - "))),
                       `if`(batch_correction_method != 'none', paste0("\n - Metadata column containing batch information: ", batch_labels), ""),
                       "\n - Nearest neighbor parameters provided: ", `if`(length(neighbor_params) == 0, "No",
                                                                           paste0("\n     - ", paste0(paste0(names(neighbor_params), ": ",
                                                                                                             neighbor_params),
                                                                                                      collapse = "\n     - "))),
                       "\n - Clustering parameters provided: ", `if`(length(cluster_params) == 0, "No",
                                                                     paste0("\n     - ", paste0(paste0(names(cluster_params), ": ",
                                                                                                       cluster_params),
                                                                                                collapse = "\n     - "))),
                       "\n - # of cores: ", n_cores,
                       "\n - Random seed: ", random_seed,
                       "\n")

  # ---------------------------------------------------------------------------
  # Step 2: Initial dimensionality reduction
  # ---------------------------------------------------------------------------

  # Run dimensionality reduction if not supplied by user
  if (is.null(reduction)) {
    if (verbose  & max_clusters == "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/7) Running initial dimensionality reduction..")
    if (verbose  & max_clusters != "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/5) Running initial dimensionality reduction..")
    P0_dim_reduction <- .runDimReduction(object = object,
                                         normalization_method = normalization_method,
                                         reduction_method = reduction_method,
                                         reduction_params = reduction_params,
                                         n_var_features = n_var_features,
                                         batch_correction_method = batch_correction_method,
                                         batch_correction_params = batch_correction_params,
                                         batch_labels = batch_labels,
                                         use_assay = use_assay_build,
                                         use_slot = use_slot_build,
                                         ArchR_matrix = ArchR_matrix_build,
                                         ArchR_depthcol = ArchR_depthcol,
                                         atac = atac,
                                         return_full = methods::is(object, "ArchRProject"),
                                         n_cores = n_cores,
                                         random_seed = random_seed,
                                         verbose = verbose)
  } else {
    if (verbose  & max_clusters == "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/7) Setting initial dimensionality reduction..")
    if (verbose  & max_clusters != "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/5) Setting initial dimensionality reduction..")
    P0_dim_reduction <- list("reduction_coords" = reduction,
                             "var_features" = var_features)
  }
  # Store output in object
  object <- .storeData(object, key, "reduction", P0_dim_reduction[["reduction_coords"]], "P0_reduction")
  object <- .storeData(object, key, "var_features", P0_dim_reduction[["var_features"]], "P0_var_features")
  object <- .storeData(object, key, "cell_IDs", cell_IDs, "P0_cell_IDs")
  # Store full reduction
  if (methods::is(object, "ArchRProject")) {
    object@reducedDims$CHOIR_P0_reduction <- P0_dim_reduction[["full_reduction"]]
    # Clean up
    P0_dim_reduction[["full_reduction"]] <- NULL
  } else {
    object <- .storeData(object = object,
                         key = key,
                         type = "full_reduction",
                         input_data = P0_dim_reduction[["reduction_coords"]],
                         name = "CHOIR_P0_reduction",
                         reduction_method = reduction_method,
                         use_assay = use_assay_build,
                         atac = atac)
  }
  # Clean up
  P0_dim_reduction[["P0_cell_IDs"]] <- NULL

  # ---------------------------------------------------------------------------
  # Step 2: Find nearest neighbors & calculate distance matrix for dimensionality reduction
  # ---------------------------------------------------------------------------
  if (verbose & max_clusters == "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 3/7) Generating initial nearest neighbors graph..")
  if (verbose & max_clusters != "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 3/5) Generating initial nearest neighbors graph..")

  # 1 vs. multiple dimensionality reductions
  if (n_modalities < 2 | methods::is(object, "ArchRProject")) {
    # Number & names of cells
    n_cells <- nrow(P0_dim_reduction[["reduction_coords"]])
    cell_IDs <- rownames(P0_dim_reduction[["reduction_coords"]])

    # Find neighbors
    P0_nearest_neighbors <- do.call(Seurat::FindNeighbors, c(list("object" = P0_dim_reduction[["reduction_coords"]]),
                                                             neighbor_params))
    # Dimensionality reduction distance matrix
    if (distance_approx == FALSE) {
      P0_reduction_dist <- stats::dist(P0_dim_reduction[["reduction_coords"]])
      object <- .storeData(object, key, "reduction", P0_reduction_dist, "P0_reduction_dist")
    }
  } else {
    # Number & names of cells
    n_cells <- nrow(P0_dim_reduction[["reduction_coords"]][[1]])
    cell_IDs <- rownames(P0_dim_reduction[["reduction_coords"]][[1]])
    # Prep for finding neighbors
    tmp <- matrix(stats::rnorm(nrow(P0_dim_reduction[["reduction_coords"]][[1]]) * 3, 10), ncol = nrow(P0_dim_reduction[["reduction_coords"]][[1]]), nrow = 3)
    colnames(tmp) <- rownames(P0_dim_reduction[["reduction_coords"]][[1]])
    rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
    tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = 'tmp')
    dim_list = vector("list", length = n_modalities)
    for (i in 1:n_modalities) {
      tmp_seurat[[paste0("DR_", i)]] <- Seurat::CreateDimReducObject(embeddings = P0_dim_reduction[["reduction_coords"]][[i]],
                                                                     key = paste0("DR_", i, "_"), assay = 'tmp')
      dim_list[[i]] <- 1:ncol(P0_dim_reduction[["reduction_coords"]][[i]])
    }
    # Find neighbors
    P0_nearest_neighbors <- do.call(Seurat::FindMultiModalNeighbors, c(list("object" = tmp_seurat,
                                                                            "reduction.list" = list(paste0("DR_", seq(1, n_modalities))),
                                                                            "dim.list" = dim_list,
                                                                            "knn.graph.name" = "nn",
                                                                            "snn.graph.name" = "snn"),
                                                                       neighbor_params))@graphs
    # Dimensionality reduction distance matrix
    if (distance_approx == FALSE) {
      P0_reduction_dist <- .getMultiModalDistance(tmp_seurat,
                                                  reduction_list = list(paste0("DR_", seq(1, n_modalities))),
                                                  dim_list = dim_list)
      object <- .storeData(object, key, "reduction", P0_reduction_dist, "P0_reduction_dist")
    }
    # Clean up
    rm(tmp)
    rm(tmp_seurat)
  }
  # Store output in object
  object <- .storeData(object, key, "graph", P0_nearest_neighbors[["nn"]], "P0_graph_nn")
  object <- .storeData(object, key, "graph", P0_nearest_neighbors[["snn"]], "P0_graph_snn")

  # ---------------------------------------------------------------------------
  # Step 3: Identify starting clustering resolution
  # ---------------------------------------------------------------------------
  if (verbose & max_clusters == "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 4/7) Identify starting clustering resolution..")
  if (verbose & max_clusters != "auto") message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 4/5) Identify starting clustering resolution..")

  P0_starting_resolution <- .getStartingResolution(snn_matrix = P0_nearest_neighbors[["snn"]],
                                                   cluster_params = cluster_params,
                                                   random_seed = random_seed,
                                                   verbose = verbose)

  # ---------------------------------------------------------------------------
  # Step 4: Build root tree
  # ---------------------------------------------------------------------------

  # Create dataframe to store tree generation records
  tree_records <- data.frame(tree_type = NULL,
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
                             stop_branching_reason = NULL)

  if (max_clusters == "auto") {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 5/7) Building root clustering tree..")
    P0_tree_list <- .getTree(snn_matrix = P0_nearest_neighbors[["snn"]],
                             dist_matrix = `if`(distance_approx == FALSE, P0_reduction_dist, NULL),
                             reduction = `if`(distance_approx == TRUE, P0_dim_reduction[["reduction_coords"]], NULL),
                             distance_approx = distance_approx,
                             tree_type = "silhouette",
                             cluster_params = cluster_params,
                             starting_resolution = P0_starting_resolution[["starting_resolution"]],
                             res0_clusters = P0_starting_resolution[["res0_clusters"]],
                             decimal_places = P0_starting_resolution[["decimal_places"]],
                             tree_records = tree_records,
                             n_cores = n_cores,
                             random_seed = random_seed)
    P0_tree <- P0_tree_list[["cluster_tree"]]
    if (verbose) message("\n                      ", P0_tree_list[["stop_reason"]])

    # Store records
    tree_records <- P0_tree_list[["tree_records"]]

    # Clean up
    rm(P0_tree_list)

    # Optimize tree
    if (ncol(P0_tree) >= 2) {
      P0_tree <- .optimizeTree(cluster_tree = P0_tree,
                               dist_matrix = `if`(distance_approx == FALSE, P0_reduction_dist, NULL),
                               reduction = `if`(distance_approx == TRUE, P0_dim_reduction[["reduction_coords"]], NULL),
                               distance_approx = distance_approx)
      P0_tree$CellID <- rownames(P0_tree)
      P0_tree <- P0_tree[match(cell_IDs, P0_tree$CellID), ]
      P0_tree <- dplyr::select(P0_tree, -CellID)
    }

    # Reassign column names
    colnames(P0_tree) <- paste0("L", seq(1, ncol(P0_tree)))
    # Rename cluster IDs with parent cluster & level indicators
    for (col in 1:ncol(P0_tree)) {
      P0_tree[,col] <- paste0("P0_", colnames(P0_tree)[col], "_", P0_tree[,col])
    }
  } else {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 5/5) Building clustering tree..")
    full_tree_list <- .getTree(snn_matrix = P0_nearest_neighbors[["snn"]],
                               tree_type = "full",
                               max_clusters = max_clusters,
                               cluster_params = cluster_params,
                               starting_resolution = P0_starting_resolution[["starting_resolution"]],
                               res0_clusters = P0_starting_resolution[["res0_clusters"]],
                               decimal_places = P0_starting_resolution[["decimal_places"]],
                               tree_records = tree_records,
                               n_cores = n_cores,
                               random_seed = random_seed)
    full_tree <- full_tree_list[["cluster_tree"]]
    if (verbose) message("\n                      ", full_tree_list[["stop_reason"]])

    # Store records
    tree_records <- full_tree_list[["tree_records"]]

    # Clean up
    rm(full_tree_list)

    # Optimize tree
    if (ncol(full_tree) >= 2) {
      full_tree <- .optimizeTree(cluster_tree = full_tree,
                                 dist_matrix = `if`(distance_approx == FALSE, P0_reduction_dist, NULL),
                                 reduction = `if`(distance_approx == TRUE, P0_dim_reduction[["reduction_coords"]], NULL),
                                 distance_approx = distance_approx)
      full_tree$CellID <- rownames(full_tree)
      full_tree <- full_tree[match(cell_IDs, full_tree$CellID), ]
      full_tree <- dplyr::select(full_tree, -CellID)
    }

    # Reassign column names
    colnames(full_tree) <- paste0("L", seq(1, ncol(full_tree)))
    # Rename cluster IDs with parent cluster & level indicators
    for (col in 1:ncol(full_tree)) {
      full_tree[,col] <- paste0("P0_", colnames(full_tree)[col], "_", full_tree[,col])
    }
  }
  # Clean up
  rm(P0_nearest_neighbors)
  rm(P0_starting_resolution)
  if (distance_approx == FALSE) {
    rm(P0_reduction_dist)
  }

  # ---------------------------------------------------------------------------
  # Step 5: Subcluster tree
  # ---------------------------------------------------------------------------
  if (max_clusters == "auto") {

    P0_clusters <- unique(P0_tree[, ncol(P0_tree)])
    P0_clusters_n <- length(P0_clusters)

    if (verbose) message("                      Identified ", P0_clusters_n, " clusters in root tree.")
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 6/7) Subclustering root tree..")

    # Initiate list of subtrees
    subtree_list <- vector("list", P0_clusters_n)
    max_columns <- 1

    # For each cluster at the lowest level of the tree, generate subclustered tree
    # Progress bar
    pb <- progress::progress_bar$new(format = "Generating subtrees.. [:bar] :percent in :elapsedfull",
                                     total = 100, clear = FALSE)

    subtree_weights <- (table(P0_tree[, ncol(P0_tree)])/length(cell_IDs))*100

    pb$tick(0)
    percent_done <- 0
    start_time <- Sys.time()
    hour_start_time <- Sys.time()

    # Record subtree sizes
    subtree_sizes <- c()
    total_clusters <- P0_clusters_n

    for (i in 1:P0_clusters_n) {
      # ---------------------------------------------------------------------
      # Step 5A: Subset cluster i & obtain dimensionality reduction
      # ---------------------------------------------------------------------
      P_i <- P0_clusters[i]
      subset_i <- P0_tree[, ncol(P0_tree)] == P_i
      cell_IDs_i <- rownames(P0_tree)[subset_i]
      subtree_sizes <- c(subtree_sizes, length(cell_IDs_i))
      if (length(cell_IDs_i) <= 3) {
        # If there are 3 cells or fewer in cluster i, do not proceed
        P_i_tree <- data.frame(L = rep(0, length(cell_IDs_i)))
        colnames(P_i_tree) <- paste0("L", ncol(P0_tree) + 1)
        rownames(P_i_tree) <- cell_IDs_i
        # Progress
        tick_amount <- subtree_weights[P_i]
        if (verbose) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (Subtree ", i,
                            "/", P0_clusters_n, ", ", length(cell_IDs_i), " cells), ",
                            total_clusters, " total clusters."))
        }
        pb$tick(tick_amount)
        percent_done <- percent_done + tick_amount
      } else {
        # Use subset of dimensionality reduction or recalculate
        # If there are less than 20 cells in cluster i, do not generate new dimensionality reduction
        if (subtree_reductions == FALSE | length(cell_IDs_i) < 20) {
          P_i_dim_reduction <- list("reduction_coords" = P0_dim_reduction[["reduction_coords"]][cell_IDs_i,],
                                    "var_features" = P0_dim_reduction[["var_features"]])
        } else if (subtree_reductions == TRUE) {
          P_i_dim_reduction <- .runDimReduction(object = object,
                                                normalization_method = normalization_method,
                                                reduction_method = reduction_method,
                                                reduction_params = reduction_params,
                                                n_var_features = n_var_features,
                                                batch_correction_method = batch_correction_method,
                                                batch_correction_params = batch_correction_params,
                                                batch_labels = batch_labels,
                                                use_assay = use_assay_build,
                                                use_slot = use_slot_build,
                                                ArchR_matrix = ArchR_matrix_build,
                                                ArchR_depthcol = ArchR_depthcol,
                                                atac = atac,
                                                use_cells = cell_IDs_i,
                                                n_cores = n_cores,
                                                random_seed = random_seed,
                                                verbose = FALSE)
        }
        # Store output in original object
        object <- .storeData(object, key, "reduction", P_i_dim_reduction[["reduction_coords"]], paste0("P", i, "_reduction"))
        object <- .storeData(object, key, "var_features", P_i_dim_reduction[["var_features"]], paste0("P", i, "_var_features"))
        object <- .storeData(object, key, "cell_IDs", cell_IDs_i, paste0("P", i, "_cell_IDs"))
        # Extract input matrix for random forest comparisons
        if (n_modalities > 1) {
          input_matrix_list <- vector("list", n_modalities)
          for (m in 1:n_modalities) {
            # Match input arguments
            use_assay_prune_m <- .matchArg(use_assay_prune, m)
            use_slot_prune_m <- .matchArg(use_slot_prune, m)
            ArchR_matrix_prune_m <- .matchArg(ArchR_matrix_prune, m)
            normalization_method_m <- .matchArg(normalization_method, m)
            # Variable features
            if (feature_set == "var") {
              use_features_m <- P_i_dim_reduction[["var_features"]][[m]]
            } else {
              use_features_m <- NULL
            }
            input_matrix_list[[m]] <- .getMatrix(object = object,
                                                 use_assay = use_assay_prune_m,
                                                 use_slot = use_slot_prune_m,
                                                 ArchR_matrix = ArchR_matrix_prune_m,
                                                 use_features = use_features_m,
                                                 use_cells = cell_IDs_i,
                                                 verbose = FALSE)
            if (normalization_method_m == "SCTransform") {
              input_matrix_list[[m]] <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(input_matrix_list[[m]]),
                                                                             return.only.var.genes = FALSE,
                                                                             seed.use = random_seed,
                                                                             verbose = FALSE)@assays$SCT@scale.data)
            }
          }
          input_matrix <- do.call(rbind, input_matrix_list)
        } else {
          input_matrix <- .getMatrix(object = object,
                                     use_assay = use_assay_prune,
                                     use_slot = use_slot_prune,
                                     ArchR_matrix = ArchR_matrix_prune,
                                     use_features = `if`(feature_set == "var",
                                                         P_i_dim_reduction[["var_features"]],
                                                         NULL),
                                     use_cells = cell_IDs_i,
                                     verbose = FALSE)
          if (normalization_method == "SCTransform") {
            input_matrix <- suppressWarnings(Seurat::SCTransform(Seurat::CreateSeuratObject(input_matrix),
                                                                 return.only.var.genes = FALSE,
                                                                 seed.use = random_seed,
                                                                 verbose = FALSE)@assays$SCT@scale.data)
          }
        }
        # Progress
        tick_amount <- subtree_weights[P_i]*0.1
        if (verbose) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (Subtree ", i,
                            "/", P0_clusters_n, ", ", length(cell_IDs_i), " cells), ",
                            total_clusters, " total clusters."))
        }
        pb$tick(tick_amount)
        percent_done <- percent_done + tick_amount

        # ---------------------------------------------------------------------
        # Step 5B: Recalculate nearest neighbors & distance matrix for dimensionality reduction
        # ---------------------------------------------------------------------

        # 1 vs. multiple dimensionality reductions
        if (n_modalities < 2 | methods::is(object, "ArchRProject")) {
          # Find neighbors
          P_i_nearest_neighbors <- do.call(Seurat::FindNeighbors, c(list("object" = P_i_dim_reduction[["reduction_coords"]]),
                                                                    neighbor_params))
          # Dimensionality reduction distance matrix
          if (distance_approx == FALSE) {
            P_i_reduction_dist <- stats::dist(P_i_dim_reduction[["reduction_coords"]])
            object <- .storeData(object, key, "reduction", P_i_reduction_dist, paste0("P", i, "_reduction_dist"))
          }
        } else {
          # Prep for finding neighbors
          tmp <- matrix(stats::rnorm(nrow(P_i_dim_reduction[["reduction_coords"]][[1]]) * 3, 10), ncol = nrow(P_i_dim_reduction[["reduction_coords"]][[1]]), nrow = 3)
          colnames(tmp) <- rownames(P_i_dim_reduction[["reduction_coords"]][[1]])
          rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
          tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = 'tmp')
          dim_list = vector("list", length = n_modalities)
          for (m in 1:n_modalities) {
            tmp_seurat[[paste0("DR_", m)]] <- Seurat::CreateDimReducObject(embeddings = P_i_dim_reduction[["reduction_coords"]][[m]],
                                                                           key = paste0("DR_", m, "_"), assay = 'tmp')
            dim_list[[m]] <- 1:ncol(P0_dim_reduction[["reduction_coords"]][[m]])
          }
          # Find neighbors
          P_i_nearest_neighbors <- do.call(Seurat::FindMultiModalNeighbors, c(list("object" = tmp_seurat,
                                                                                   "reduction.list" = list(paste0("DR_", seq(1, n_modalities))),
                                                                                   "dim.list" = dim_list,
                                                                                   "knn.graph.name" = "nn",
                                                                                   "snn.graph.name" = "snn"),
                                                                              neighbor_params))@graphs
          # Dimensionality reduction distance matrix
          if (distance_approx == FALSE) {
            P_i_reduction_dist <- .getMultiModalDistance(tmp_seurat,
                                                         reduction_list = list(paste0("DR_", seq(1, n_modalities))),
                                                         dim_list = dim_list)
            object <- .storeData(object, key, "reduction", P_i_reduction_dist, paste0("P", i, "_reduction_dist"))
          }
          # Clean up
          rm(tmp)
          rm(tmp_seurat)
        }
        # Store output in original object
        object <- .storeData(object, key, "graph", P_i_nearest_neighbors[["nn"]], paste0("P", i, "_graph_nn"))
        object <- .storeData(object, key, "graph", P_i_nearest_neighbors[["snn"]], paste0("P", i, "_graph_snn"))
        # Progress
        tick_amount <- subtree_weights[P_i]*0.05
        if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                       (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (Subtree ", i,
                            "/", P0_clusters_n, ", ", length(cell_IDs_i), " cells), ",
                            total_clusters, " total clusters."))
        }
        pb$tick(tick_amount)
        percent_done <- percent_done + tick_amount

        # ---------------------------------------------------------------------
        # Step 5C: Subdivide until silhouette score is again maximized
        # ---------------------------------------------------------------------

        if (length(cell_IDs_i) > min_cluster_depth) {
          # Get starting resolution
          P_i_starting_resolution <- .getStartingResolution(snn_matrix = P_i_nearest_neighbors[["snn"]],
                                                            cluster_params = cluster_params,
                                                            random_seed = random_seed,
                                                            verbose = FALSE)
          # Maximize silhouette score
          P_i_tree_list <- .getTree(snn_matrix = P_i_nearest_neighbors[["snn"]],
                                    dist_matrix = `if`(distance_approx == FALSE, P_i_reduction_dist, NULL),
                                    reduction = `if`(distance_approx == TRUE, P_i_dim_reduction[["reduction_coords"]], NULL),
                                    distance_approx = distance_approx,
                                    tree_type = "silhouette",
                                    cluster_params = cluster_params,
                                    starting_resolution = P_i_starting_resolution[["starting_resolution"]],
                                    res0_clusters = P_i_starting_resolution[["res0_clusters"]],
                                    decimal_places = P_i_starting_resolution[["decimal_places"]],
                                    tree_records = tree_records,
                                    tree_id = paste0("P", i),
                                    n_cores = n_cores,
                                    random_seed = random_seed)
          P_i_tree <- P_i_tree_list[["cluster_tree"]]

          # Store records
          tree_records <- P_i_tree_list[["tree_records"]]

          # Clean up
          rm(P_i_tree_list)

          # Optimize tree
          if (ncol(P_i_tree) >= 2) {
            P_i_tree <- .optimizeTree(cluster_tree = P_i_tree,
                                      dist_matrix = `if`(distance_approx == FALSE, P_i_reduction_dist, NULL),
                                      reduction = `if`(distance_approx == TRUE, P_i_dim_reduction[["reduction_coords"]], NULL),
                                      distance_approx = distance_approx)
          }

          # Reassign column names
          colnames(P_i_tree) <- paste0("L", seq(ncol(P0_tree) + 1,
                                                ncol(P0_tree) + ncol(P_i_tree)))
        } else {
          P_i_tree <- data.frame(L = rep(0, length(cell_IDs_i)))
          colnames(P_i_tree) <- paste0("L", ncol(P0_tree) + 1)
          rownames(P_i_tree) <- cell_IDs_i
        }

        # Clusters in current tree
        P_i_clusters <- unique(P_i_tree[,ncol(P_i_tree)])
        P_i_clusters_n <- length(P_i_clusters)

        # Progress
        tick_amount <- subtree_weights[P_i]*0.05
        if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                       (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (Subtree ", i,
                            "/", P0_clusters_n, ", ", length(cell_IDs_i), " cells), ",
                            total_clusters + P_i_clusters_n, " total clusters."))
        }
        pb$tick(tick_amount)
        percent_done <- percent_done + tick_amount

        # ---------------------------------------------------------------------
        # Step 5D: For each current cluster, subdivide until overclustered
        # ---------------------------------------------------------------------

        for (j in 1:P_i_clusters_n) {
          # Subset
          cell_IDs_j <- cell_IDs_i[P_i_tree[, ncol(P_i_tree)] == P_i_clusters[j]]
          if (length(cell_IDs_j) > 3) {
            current_snn <- P_i_nearest_neighbors[["snn"]][cell_IDs_j, cell_IDs_j]
            current_nn <- P_i_nearest_neighbors[["nn"]][cell_IDs_j, cell_IDs_j]
            current_input_matrix <- BiocGenerics::t(input_matrix[, cell_IDs_j])
            # Build subtree
            P_j_tree_list <- .getTree(snn_matrix = current_snn,
                                      nn_matrix = current_nn,
                                      dist_matrix = `if`(distance_approx == FALSE, P_i_reduction_dist[cell_IDs_j, cell_IDs_j], NULL),
                                      reduction = `if`(distance_approx == TRUE, P_i_dim_reduction[["reduction_coords"]][cell_IDs_j,], NULL),
                                      input_matrix = current_input_matrix,
                                      distance_approx = distance_approx,
                                      tree_type = "subtree",
                                      cluster_params = cluster_params,
                                      min_cluster_depth = min_cluster_depth,
                                      alpha = `if`(p_adjust == "bonferroni", alpha/(P0_clusters_n*4), alpha),
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
                                      tree_records = tree_records,
                                      tree_id = paste0("P", i, "_", j),
                                      n_cores = n_cores,
                                      random_seed = random_seed)
            P_j_tree <- P_j_tree_list[["cluster_tree"]]

            # Store records
            tree_records <- P_j_tree_list[["tree_records"]]

            # Optimize tree
            if (ncol(P_j_tree) >= 2) {
              P_j_tree <- .optimizeTree(cluster_tree = P_j_tree,
                                        dist_matrix = `if`(distance_approx == FALSE, P_i_reduction_dist[cell_IDs_j, cell_IDs_j], NULL),
                                        reduction = `if`(distance_approx == TRUE, P_i_dim_reduction[["reduction_coords"]][cell_IDs_j, ], NULL),
                                        distance_approx = distance_approx)
            }

            # Reassign column names
            colnames(P_j_tree) <- paste0("L", seq(ncol(P0_tree) + ncol(P_i_tree) + 1,
                                                  ncol(P0_tree) + ncol(P_i_tree) + ncol(P_j_tree)))

            # Add to main subtree
            for (col in (ncol(P0_tree) + ncol(P_i_tree) + 1):(ncol(P0_tree) + ncol(P_i_tree) + ncol(P_j_tree))) {
              level_name <- paste0("L", col)
              if (!(level_name %in% colnames(P_i_tree))) {
                P_i_tree[, level_name] <- P_i_tree[, paste0("L", col-1)]
              }
              P_i_tree[cell_IDs_j, level_name] <- P_j_tree[cell_IDs_j, level_name] + max(P_i_tree[, level_name]) + 1
            }
          }

          # Progress
          tick_amount <- subtree_weights[P_i]*0.75*(1/P_i_clusters_n)
          if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                         (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
            hour_start_time <- Sys.time()
            pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                              " : ", round((percent_done + tick_amount)), "% (Subtree ", i,
                              "/", P0_clusters_n, ", ", length(cell_IDs_i), " cells), ",
                              total_clusters + dplyr::n_distinct(P_i_tree[,ncol(P_i_tree)]),
                              " total clusters."))
          }
          pb$tick(tick_amount)
          percent_done <- percent_done + tick_amount
        }

        # Cluster IDs to sequential integers
        for (col in 1:ncol(P_i_tree)) {
          P_i_tree[, col] <- as.numeric(as.factor(P_i_tree[, col]))
        }

        # Get max number of columns
        max_columns <- max(max_columns, ncol(P_i_tree))

        # Progress
        tick_amount <- subtree_weights[P_i]*0.05
        if (verbose & ((((percent_done + tick_amount) %/% 10) - (percent_done %/% 10) > 0) |
                       (difftime(Sys.time(), hour_start_time, units = "hours") >= 0.5))) {
          hour_start_time <- Sys.time()
          pb$message(paste0(format(Sys.time(), "%Y-%m-%d %X"),
                            " : ", round((percent_done + tick_amount)), "% (Subtree ", i,
                            "/", P0_clusters_n, ", ", length(cell_IDs_i), " cells), ",
                            total_clusters + dplyr::n_distinct(P_i_tree[,ncol(P_i_tree)]),
                            " total clusters."))
        }
        pb$tick(tick_amount)
        percent_done <- percent_done + tick_amount
      }
      # Add to list
      subtree_list[[i]] <- P_i_tree
      total_clusters <- total_clusters + dplyr::n_distinct(P_i_tree[,ncol(P_i_tree)])
    }
    # Clean up
    rm(P0_dim_reduction)
    rm(P_i_tree)
    if (exists('P_i_tree_list')) {
      rm(P_i_tree_list)
    }
    rm(P_j_tree)
    rm(P_j_tree_list)

    # ---------------------------------------------------------------------------
    # Step 6: Stitch trees together
    # ---------------------------------------------------------------------------
    if (verbose) message("\n", format(Sys.time(), "%Y-%m-%d %X"), " : (Step 7/7) Compiling full clustering tree..")

    # Create shell dataframe for subtree clusters
    subtrees <- data.frame(matrix(rep(NA, max_columns + 1), nrow = 1, ncol = max_columns + 1))
    colnames(subtrees) <- c(paste0("L", seq(ncol(P0_tree) + 1, ncol(P0_tree) + max_columns)),
                            "CellID")

    # Populate cluster IDs up to max_columns
    for (i in 1:P0_clusters_n) {
      subtree_i <- data.frame(subtree_list[[i]])
      if (ncol(subtree_i) < max_columns) {
        for (col in seq(ncol(P0_tree) + ncol(subtree_i) + 1,
                        ncol(P0_tree) + max_columns)) {
          subtree_i[, paste0("L", col)] <- subtree_i[, ncol(subtree_i)]
        }
      }
      # Rename cluster IDs with parent cluster & level indicators
      for (col in 1:ncol(subtree_i)) {
        subtree_i[,col] <- paste0("P", i, "_", colnames(subtree_i)[col], "_", subtree_i[,col])
      }
      subtree_i$CellID <- rownames(subtree_i)
      subtrees <- rbind(subtrees, subtree_i)
    }

    # Merge subtrees with P0 tree
    subtrees <- subtrees[-1,]
    P0_tree$CellID <- rownames(P0_tree)
    P0_tree <- P0_tree[match(cell_IDs, P0_tree$CellID), ]
    subtrees <- subtrees[match(cell_IDs, subtrees$CellID), ]
    full_tree <- cbind(dplyr::select(P0_tree, -CellID), dplyr::select(subtrees, -CellID))

    # Clean up
    rm(P0_tree)
    rm(subtrees)
  }

  if (verbose) message("                      Full tree has ", ncol(full_tree), " levels and ", dplyr::n_distinct(full_tree[,ncol(full_tree)]), " clusters.")

  tree_records <- rbind(tree_records, data.frame(tree_type = "full",
                                                 tree_name = "P0",
                                                 num_cells = length(cell_IDs),
                                                 resolution = NA,
                                                 num_clusters = dplyr::n_distinct(full_tree[,ncol(full_tree)]),
                                                 silhouette = NA,
                                                 neighbors_distance = NA,
                                                 neighbors_mean_accuracy = NA,
                                                 neighbors_var_accuracy = NA,
                                                 neighbors_percentile_accuracy = NA,
                                                 neighbors_percentile_variance = NA,
                                                 neighbors_decision = NA,
                                                 stop_branching_reason = "Full tree compiled."))

  # ---------------------------------------------------------------------------
  # Store data in object
  # ---------------------------------------------------------------------------

  # Add full tree to original object
  object <- .storeData(object, key, "clusters", full_tree, "full_tree")

  # Add tree records to object
  object <- .storeData(object, key, "records", tree_records, "buildTree_records")

  # Record parameters used and add to original object
  parameter_list <- list("subtree_names" = `if`(max_clusters == "auto", c("P0", paste0("P", seq(1, P0_clusters_n))), "P0"),
                         "subtree_sizes" = `if`(max_clusters == "auto", c(n_cells, subtree_sizes), NULL),
                         "alpha" = alpha,
                         "p_adjust" = p_adjust,
                         "feature_set" = feature_set,
                         "exclude_features" = exclude_features,
                         "n_iterations" = n_iterations,
                         "n_trees" = n_trees,
                         "use_variance" = use_variance,
                         "min_accuracy" = min_accuracy,
                         "min_connections" = min_connections,
                         "max_repeat_errors" = max_repeat_errors,
                         "distance_approx"  = distance_approx,
                         "sample_max" = sample_max,
                         "downsampling_rate" = downsampling_rate,
                         "min_reads" = min_reads,
                         "max_clusters" = max_clusters,
                         "min_cluster_depth" = min_cluster_depth,
                         "normalization_method" = normalization_method,
                         "subtree_reductions" = subtree_reductions,
                         "reduction_method" = reduction_method,
                         "reduction_params" = reduction_params,
                         "n_var_features" = n_var_features,
                         "batch_correction_method" = batch_correction_method,
                         "batch_correction_params" = batch_correction_params,
                         "batch_labels" = batch_labels,
                         "neighbor_params" = neighbor_params,
                         "cluster_params" = cluster_params,
                         "use_assay" = use_assay,
                         "use_slot" = use_slot,
                         "ArchR_matrix" = ArchR_matrix,
                         "ArchR_depthcol" = ArchR_depthcol,
                         "countsplit" = countsplit,
                         "countsplit_suffix" = countsplit_suffix,
                         "use_assay_build" = use_assay_build,
                         "use_assay_prune" = use_assay_prune,
                         "use_slot_build" = use_slot_build,
                         "use_slot_prune" = use_slot_prune,
                         "ArchR_matrix_build" = ArchR_matrix_build,
                         "ArchR_matrix_prune" = ArchR_matrix_prune,
                         "countsplit_text" = countsplit_text,
                         "reduction_provided" = !is.null(reduction),
                         "var_features_provided" = !is.null(var_features),
                         "atac" = atac,
                         "random_seed" = random_seed)

  object <- .storeData(object, key, "parameters", parameter_list, "buildTree_parameters")

  # Return object with new additions
  return(object)
}
