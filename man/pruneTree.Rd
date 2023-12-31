% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruneTree.R
\name{pruneTree}
\alias{pruneTree}
\title{Prune clustering tree using random forest classifiers}
\usage{
pruneTree(
  object,
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
  normalization_method = NULL,
  batch_correction_method = NULL,
  batch_labels = NULL,
  cluster_params = NULL,
  use_assay = NULL,
  cluster_tree = NULL,
  input_matrix = NULL,
  nn_matrix = NULL,
  dist_matrix = NULL,
  reduction = NULL,
  n_cores = NULL,
  random_seed = NULL,
  verbose = TRUE
)
}
\arguments{
\item{object}{An object of class 'Seurat', 'SingleCellExperiment', or
'ArchRProject'.}

\item{key}{The name under which CHOIR-related data for this run is stored in
the object. Defaults to 'CHOIR'.}

\item{alpha}{A numerical value indicating the significance level used for
permutation test comparisons of cluster distinguishability. Defaults to 0.05.}

\item{p_adjust}{A string indicating which multiple comparison
adjustment to use. Permitted values are 'fdr', 'bonferroni', and 'none'.
Defaults to 'bonferroni'.}

\item{feature_set}{A string indicating whether to train random forest
classifiers on 'all' features or only variable ('var') features. Defaults to 'var'.}

\item{exclude_features}{A character vector indicating features that should be
excluded from input to the random forest classifier. Default = \code{NULL}
will not exclude any features.}

\item{n_iterations}{A numeric value indicating the number of iterations run
for each permutation test comparison. Defaults to 100.}

\item{n_trees}{A numeric value indicating the number of trees in each random
forest. Defaults to 50.}

\item{use_variance}{A boolean value indicating whether to use the variance of
the random forest accuracy scores as part of the permutation test threshold.
Defaults to \code{TRUE}. If only ATAC-seq data is supplied, it is recommended
to set to \code{FALSE}.}

\item{min_accuracy}{A numeric value indicating the minimum accuracy required
of the random forest classifier, below which clusters will be automatically
merged. Defaults to 0.5.}

\item{min_connections}{A numeric value indicating the minimum number of
nearest neighbors between two clusters for them to be considered 'adjacent'.
Non-adjacent clusters will not be merged. Defaults to 1.}

\item{max_repeat_errors}{Used to account for situations in which random
forest classifier errors are concentrated among a few cells that are
repeatedly misassigned. Numeric value indicating the maximum number of such
'repeat errors' that will be taken into account. If set to 0, 'repeat errors'
will not be evaluated. Defaults to 20.}

\item{distance_approx}{A boolean value indicating whether or not to use
approximate distance calculations. Default = TRUE will use centroid-based
distances.}

\item{distance_awareness}{To omit all distance calculations, set to
\code{FALSE}. Otherwise, a numeric value representing the distance threshold
above which a cluster will not merge with another cluster. Specifically,
this value is a multiplier applied to the distance between a cluster and its
closest distinguishable neighbor, giving the threshold. Default = 2 sets
this threshold at a 2-fold increase in distance.}

\item{collect_all_metrics}{A boolean value indicating whether to collect and save
additional metrics from the random forest classifiers, including feature
importances and tree depth. Defaults to \code{FALSE}.}

\item{sample_max}{A numeric value indicating the maximum number of cells used
per cluster to train/test each random forest classifier. Default = \code{Inf}
does not cap the number of cells used.}

\item{downsampling_rate}{A numeric value indicating the proportion of cells used
per cluster to train/test each random forest classifier. Default = "auto" sets
the downsampling rate according to the dataset size, for efficiency.}

\item{normalization_method}{A character string or vector indicating which
normalization method to use. In general, input data should be supplied to
CHOIR after normalization, except in cases when the user wishes to use
\code{Seurat::SCTransform()} normalization. Permitted values are 'none' or
'SCTransform'. Defaults to 'none'.}

\item{batch_correction_method}{A character string or vector indicating which
batch correction method to use. Permitted values are 'Harmony', 'VAE', and
'none'. Defaults to 'none'.}

\item{batch_labels}{If applying batch correction, the name of the column
containing the batch labels. Defaults to \code{NULL}.}

\item{cluster_params}{A list of additional parameters to be passed to
Seurat::FindClusters() for clustering at each level of the tree. Note that if
'group.singletons' is set to TRUE, clusters are relabeled such that each
singleton constitutes its own cluster.}

\item{use_assay}{For Seurat or SingleCellExperiment objects, a character
string or vector indicating the assay(s) to use in the provided object.
Default = \code{NULL} will choose the current active assay for Seurat objects
and the \code{log_counts} assay for SingleCellExperiment objects.}

\item{cluster_tree}{An optional dataframe containing the cluster IDs of each
cell across the levels of a hierarchical clustering tree. Default = \code{NULL}
will use the hierarchical clustering tree generation by function
\code{buildTree()}.}

\item{input_matrix}{An optional matrix containing the feature x cell data on
which to train the random forest classifiers. Default = \code{NULL} will use
the feature x cell matri(ces) indicated by function \code{buildTree()}.}

\item{nn_matrix}{An optional matrix containing the nearest neighbor adjacency
of the cells. Default = \code{NULL} will look for the adjacency matri(ces)
generated by function \code{buildTree()}.}

\item{dist_matrix}{An optional distance matrix of cell to cell distances (based
on dimensionality reduction cell embeddings). Default = \code{NULL} will look
for the distance matri(ces) generated by function \code{buildTree()}.}

\item{reduction}{An optional matrix of dimensionality reduction cell
embeddings to be used for distance calculations. Defaults = \code{NULL} will
look for the dimensionality reductions generated by function \code{buildTree()}.}

\item{n_cores}{A numeric value indicating the number of cores to use for
parallelization. Default = \code{NULL} will use the number of available cores
minus 2.}

\item{random_seed}{A numeric value indicating the random seed to be used.}

\item{verbose}{A boolean value indicating whether to use verbose output
during the execution of this function. Can be set to \code{FALSE} for a
cleaner output.}
}
\value{
Returns the object with the following added data stored under the provided
key: \describe{
\item{final_clusters}{Clustering results for each provided threshold value}
\item{stepwise_clusters}{A dataframe of clustering results for each progressive step through the
clustering tree at each provided threshold value}
\item{comparison_records}{A dataframe of all recorded comparisons}
\item{feature_importances}{If 'collect_all_metrics' is true, a dataframe containing the feature importance scores for each feature across all comparisons}
}
}
\description{
This function will move iteratively up the provided hierarchical clustering
tree, and use permutation tests of random forest classifier prediction
accuracies to identify which clusters should be merged, in order to identify
robust final clusters.
}
