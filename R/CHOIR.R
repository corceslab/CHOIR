#' Run CHOIR clustering
#'
#' This function runs CHOIR clustering to identify a set of clusters that
#' represent statistically distinct populations. CHOIR is based on the premise
#' that, if clusters contain biologically different cell types or states, a
#' classifier that considers features present in cells from each cluster should
#' be able to distinguish the clusters with a higher level of accuracy than
#' classifiers trained on randomly permuted cluster labels.
#'
#' First, a hierarchical clustering tree is generated using a top-down approach
#' that proceeds from an initial partition, in which all cells are in the same
#' cluster, to a partition in which all cells are demonstrably overclustered.
#' Second, to identify a final set of clusters, this hierarchical clustering
#' tree is pruned from the bottom up using a framework of random forest
#' classifiers and permutation tests.
#'
#' For multi-modal data, optionally supply parameter inputs as vectors/lists
#' that sequentially specify the value for each modality.
#'
#' @param object An object of class 'Seurat', 'SingleCellExperiment', or
#' 'ArchRProject'.
#' @param key The name under which CHOIR-related data for this run is stored in
#' the object. Defaults to 'CHOIR'.
#' @param alpha A numeric value indicating the significance level used for
#' permutation test comparisons of cluster prediction accuracies. Defaults to
#' 0.05.
#' @param p_adjust A string indicating which multiple comparison adjustment to
#' use. Permitted values are 'bonferroni', 'fdr', and 'none'. Defaults to
#' 'bonferroni'.
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
#' @param distance_approx A boolean value indicating whether or not to use
#' approximate distance calculations. Default = \code{TRUE} will use
#' centroid-based distances.
#' @param distance_awareness A numeric value representing the distance threshold
#' above which a cluster will not merge with another cluster. Specifically,
#' this value is multiplied by the distance between a cluster and its
#' closest distinguishable neighbor to set the threshold. Default = 2 sets
#' this threshold at a 2-fold increase in distance. Alternately, to omit all
#' distance calculations, set to \code{FALSE}.
#' @param collect_all_metrics A boolean value indicating whether to collect and
#' save additional metrics from the random forest classifier comparisons,
#' including feature importances and tree depth. Defaults to \code{FALSE}.
#' @param sample_max A numeric value indicating the maximum number of cells used
#' per cluster to train/test each random forest classifier. Default = \code{Inf}
#' does not cap the number of cells used.
#' @param downsampling_rate A numeric value indicating the proportion of cells
#' used per cluster to train/test each random forest classifier. Default =
#' "auto" sets the downsampling rate according to the dataset size, for
#' efficiency.
#' @param max_clusters Indicates the extent to which the hierarchical clustering
#' tree will be expanded. Default = 'auto' will expand the tree until instances
#' of underclustering have been eliminated in all branches. Alternately, supply
#' a numeric value indicating the maximum number of clusters to expand the
#' tree to.
#' @param min_cluster_depth A numeric value indicating the maximum cluster size
#' at the bottom of the clustering tree, prior to pruning branches. Defaults to
#' 2000.
#' @param normalization_method A character string or vector indicating which
#' normalization method to use. In general, input data should be supplied to
#' CHOIR after normalization, except in cases when the user wishes to use
#' \code{Seurat::SCTransform()} normalization. Permitted values are 'none' or
#' 'SCTransform'. Defaults to 'none'.
#' @param subtree_reductions A boolean value indicating whether to generate a
#' new dimensionality reduction for each subtree. Defaults to \code{TRUE}.
#' @param reduction_method A character string or vector indicating which
#' dimensionality reduction method to use. Permitted values are 'PCA' for
#' principal component analysis, 'LSI' for latent semantic indexing, and
#' 'IterativeLSI' for iterative latent semantic indexing. Default = \code{NULL}
#' will specify a method automatically based on the input data type.
#' @param reduction_params A list of additional parameters to be passed to
#' the selected dimensionality reduction method.
#' @param n_var_features A numeric value indicating how many variable features
#' to identify. Default = \code{NULL} will use 2000 features, or 25000 features
#' for ATAC-seq data.
#' @param batch_correction_method A character string or vector indicating which
#' batch correction method to use. Permitted values are 'Harmony' and
#' 'none'. Defaults to 'none'.
#' @param batch_correction_params A list of additional parameters to be passed
#' to the selected batch correction method for each iteration. Only applicable
#' when \code{batch_correction_method = 'Harmony'}.
#' @param batch_labels If applying batch correction, a character string or
#' vector indicating the name of the column containing the batch labels.
#' Defaults to \code{NULL}.
#' @param neighbor_params A list of additional parameters to be passed to
#' \code{Seurat::FindNeighbors()} (or, in the case of multi-modal data for
#' Seurat or SingleCellExperiment objects,
#' \code{Seurat::FindMultiModalNeighbors()}).
#' @param cluster_params A list of additional parameters to be passed to
#' Seurat::FindClusters() for clustering at each level of the tree. Note that if
#' \code{group.singletons} is set to \code{TRUE}, \code{CHOIR} relabels initial
#' clusters such that each singleton constitutes its own cluster.
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
#' @param ArchR_depthcol For ArchR objects, a character string or vector
#' indicating which column to use for correlation with sequencing depth.
#' Default = \code{NULL} will use the 'nFrags' column for ATAC-seq data or the
#' 'Gex_nUMI' for RNA-seq data.
#' @param countsplit A boolean value indicating whether or not to use
#' countsplit input data (see A. Neufeld \code{countsplit} package), such that
#' one matrix of counts is used for clustering tree generation, and a separate
#' matrix is used for all random forest classifier permutation testing. Defaults
#' to \code{FALSE}.
#' @param countsplit_suffix A character vector indicating the suffixes
#' that distinguish the two countsplit matrices to be used. Suffixes are
#' appended onto input string/vector for \code{use_slot} for Seurat objects,
#' \code{use_assay} for SingleCellExperiment objects, or \code{ArchR_matrix} for
#' ArchR objects. When countsplitting is enabled, default = \code{NULL} uses
#' suffixes "_1" and "_2".
#' @param reduction An optional matrix of dimensionality reduction cell
#' embeddings to be used for subsequent clustering steps. Defaults to
#' \code{NULL}, whereby dimensionality reduction(s) will instead be calculated
#' using method specified by \code{reduction_method}.
#' @param var_features An optional character vector of variable features to be
#' used for subsequent clustering steps. Defaults to \code{NULL}, whereby new
#' sets of variable features will instead be generated.
#' @param atac A boolean value or vector indicating whether the provided data is
#' ATAC-seq data. Defaults to \code{FALSE}. For multi-omic datasets containing
#' ATAC-seq data, it is important to supply this parameter as a vector
#' corresponding to each modality in order.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param random_seed A numeric value indicating the random seed to be used.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Can be set to \code{FALSE} for a
#' cleaner output.
#'
#' @return Returns the object with the following added data stored under the
#' provided key: \describe{
#'   \item{cell_IDs}{Cell IDs belonging to each subtree}
#'   \item{clusters}{Final clusters, full hierarchical cluster tree, and
#'   stepwise cluster results for each progressive pruning step}
#'   \item{graph}{All calculated nearest neighbor and shared nearest neighbor
#'   adjacency matrices}
#'   \item{parameters}{Record of parameter values used}
#'   \item{records}{Metadata for decision points during hierarchical tree
#'   construction, all recorded permutation test comparisons, and feature
#'   importance scores from all comparisons}
#'   \item{reduction}{Cell embeddings for all calculated dimensionality
#'   reductions}
#'   \item{var_features}{Variable features for all calculated dimensionality
#'   reductions}
#'   }
#'
#' @export
#'
CHOIR <- function(object,
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
                  distance_awareness = 2,
                  collect_all_metrics = FALSE,
                  sample_max = Inf,
                  downsampling_rate = "auto",
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

  .validInput(distance_awareness, "distance_awareness")
  .validInput(collect_all_metrics, "collect_all_metrics")

  # ---------------------------------------------------------------------------
  # Part 1: Build tree
  # ---------------------------------------------------------------------------
  if (verbose) {
    message("----------------------------------------")
    message("- CHOIR - Part 1: Build clustering tree")
    message("----------------------------------------")
  }
  object <- buildTree(object = object,
                      key = key,
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
                      sample_max = sample_max,
                      downsampling_rate = downsampling_rate,
                      max_clusters = max_clusters,
                      min_cluster_depth = min_cluster_depth,
                      distance_approx = distance_approx,
                      normalization_method = normalization_method,
                      subtree_reductions = subtree_reductions,
                      reduction_method = reduction_method,
                      reduction_params = reduction_params,
                      n_var_features = n_var_features,
                      batch_correction_method = batch_correction_method,
                      batch_correction_params = batch_correction_params,
                      batch_labels = batch_labels,
                      neighbor_params = neighbor_params,
                      cluster_params = cluster_params,
                      use_assay = use_assay,
                      use_slot = use_slot,
                      ArchR_matrix = ArchR_matrix,
                      ArchR_depthcol = ArchR_depthcol,
                      countsplit = countsplit,
                      countsplit_suffix = countsplit_suffix,
                      reduction = reduction,
                      var_features = var_features,
                      atac = atac,
                      n_cores = n_cores,
                      random_seed = random_seed,
                      verbose = verbose)

  # ---------------------------------------------------------------------------
  # Part 2: Prune clustering tree
  # ---------------------------------------------------------------------------
  # pruneTree function inherits parameter values used in buildTree
  if (verbose) {
    message("")
    message("----------------------------------------")
    message("- CHOIR - Part 2: Prune clustering tree")
    message("----------------------------------------")
  }
  object <- pruneTree(object = object,
                      key = key,
                      distance_awareness = distance_awareness,
                      collect_all_metrics = collect_all_metrics,
                      n_cores = n_cores,
                      verbose = verbose)

  # Output
  return(object)
}
