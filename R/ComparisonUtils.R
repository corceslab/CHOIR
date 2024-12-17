# ---------------------------------------------------------------------------
# Helper functions for permutation test comparisons
# ---------------------------------------------------------------------------

# Run a permutation test comparison ---------------------------
#
# This function is a wrapper around .runRF() and decides whether to split or
# merge two clusters based on the prediction accuracies achieved by the random
# forest classifiers on true and permuted cluster labels.
#
# cluster1_name -- A string indicating the name of the first cluster
# cluster1_cells -- A vector indicating the cells belonging to the first cluster
# cluster1_cell_batches -- A vector indicating the batch IDs for cluster 1
# cluster2_name -- A string indicating the name of the second cluster
# cluster2_cells -- A vector indicating the cells belonging to the second cluster
# cluster2_cell_batches -- A vector indicating the batch IDs for cluster 2
# alpha -- A numeric value indicating the significance level used for the permutation test
# n_iterations A numeric value indicating the number of iterations run for each permutation test
# n_trees A numeric value indicating the number of trees in each random forest
# use_variance -- A boolean value indicating whether to use variance in permutation test
# min_accuracy -- A numeric value indicating the minimum accuracy below which clusters will be automatically merged
# min_connections -- A numeric value indicating the minimum number of nearest neighbors between two clusters for them to be considered "adjacent"
# max_repeat_errors -- A numeric value indicating the maximum number of cells that will be considered as repeated errors
# collect_all_metrics -- A boolean value indicating whether to collect and save additional metrics from the random forest classifiers
# sample_max -- A numeric indicating max cells to sample for random forest
# downsampling_rate -- A numeric indicating how much to downsample cells from each cluster for train/test
# min_reads -- A numeric used to filter out features that do not have more than 1 read for this many cells in at least one of the clusters
# input_matrix -- A matrix of sequencing data, on which the random forest classifier will be trained/tested
# nn_matrix -- A nearest neighbors adjacency matrix
# comparison_records -- A dataframe containing the comparison records output of the previous runs of this function
# feature_importance_records -- A dataframe containing the feature importance output of the previous runs of this function
# P0_distance -- Distance between two clusters using root tree dimensionality reduction
# P_i_distance -- Distance between two clusters using subtree dimensionality reduction, if applicable
# comparison_start_time -- Time that comparison was begun
# n_cores -- A numeric value indicating the number of cores to use for parallelization
# random_seed -- A numeric value indicating the random seed used
#
# Returns a list containing the following 4 elements:
# comparison_result -- A string, either "merge" or "split", indicating the result of the comparison.
# comparison_records -- A dataframe of all recorded comparisons
# feature_importances -- If 'collect_all_metrics' is TRUE, a dataframe containing the feature importance scores for the comparison
# max_p -- Max p-value for comparison (for filtering if adjusted alpha threshold decreases)
.runPermutationTest <- function(cluster1_name,
                                cluster1_cells,
                                cluster1_cell_batches = NULL,
                                cluster2_name,
                                cluster2_cells,
                                cluster2_cell_batches = NULL,
                                alpha,
                                n_iterations,
                                n_trees,
                                use_variance,
                                min_accuracy,
                                min_connections,
                                max_repeat_errors,
                                collect_all_metrics,
                                sample_max,
                                downsampling_rate,
                                min_reads,
                                input_matrix,
                                nn_matrix,
                                comparison_records = NULL,
                                feature_importance_records = NULL,
                                P0_distance = NA,
                                P_i_distance = NA,
                                comparison_start_time = NULL,
                                n_cores,
                                random_seed) {

  # Set default
  max_p <- 1

  # Check number of adjacent cells
  if ((min_connections > 0 | collect_all_metrics == TRUE) & !is.null(nn_matrix)) {
    adjacent <- sum(as.matrix(nn_matrix[cluster1_cells, cluster2_cells])) +
      sum(as.matrix(nn_matrix[cluster2_cells, cluster1_cells]))
    # If clusters don't have enough adjacent cells,
    # don't compare & don't allow to merge
    if (adjacent < min_connections) {
      proceed <- FALSE
      comparison_result <- "split"
      # Add to comparison records
      if (!is.null(comparison_records)) {
        current_comparison <- matrix(rep(NA, ncol(comparison_records)), nrow = 1, ncol = ncol(comparison_records))
        colnames(current_comparison) <- colnames(comparison_records)
        current_comparison <- data.frame(current_comparison)
        current_comparison$comparison <- paste0(cluster1_name, " vs. ", cluster2_name)
        current_comparison$cluster1_size <- length(cluster1_cells)
        current_comparison$cluster2_size <- length(cluster2_cells)
        current_comparison$connectivity <- adjacent
        current_comparison$decision <- "split: min connections"
        if ("root_distance" %in% colnames(comparison_records)) {
          current_comparison$root_distance <- P0_distance
          current_comparison$subtree_distance <- P_i_distance
        }
      }
    } else {
      proceed <- TRUE
    }
  } else {
    adjacent <- NA
    proceed <- TRUE
  }

  # If relevant, make a list of batches to check
  if (!is.null(cluster1_cell_batches) & proceed == TRUE) {
    # Only use batches which fulfill these conditions in both clusters:
    # (1) Contain at least 2 cells
    # (2) Make up at least 20% as many cells as the most prevalent batch in the cluster
    cluster1_batch_freq <- table(cluster1_cell_batches)
    cluster1_batches <- names(cluster1_batch_freq)[cluster1_batch_freq >= max(2, 0.20*max(cluster1_batch_freq))]
    cluster2_batch_freq <- table(cluster2_cell_batches)
    cluster2_batches <- names(cluster2_batch_freq)[cluster2_batch_freq >= max(2, 0.20*max(cluster2_batch_freq))]
    batches <- intersect(cluster1_batches, cluster2_batches)
    if (length(batches) == 0) {
      comparison_result <- "merge"
      proceed <- FALSE
      # Add to comparison records
      if (!is.null(comparison_records)) {
        current_comparison <- matrix(rep(NA, ncol(comparison_records)), nrow = 1, ncol = ncol(comparison_records))
        colnames(current_comparison) <- colnames(comparison_records)
        current_comparison <- data.frame(current_comparison)
        current_comparison$comparison <- paste0(cluster1_name, " vs. ", cluster2_name)
        current_comparison$cluster1_size <- length(cluster1_cells)
        current_comparison$cluster2_size <- length(cluster2_cells)
        current_comparison$connectivity <- adjacent
        current_comparison$decision <- "merge: batch-dependent"
        if ("root_distance" %in% colnames(comparison_records)) {
          current_comparison$root_distance <- P0_distance
          current_comparison$subtree_distance <- P_i_distance
        }
      }
    } else {
      # Make a vector of repeating batches to evaluate
      use_batch <- c(rep(batches, round(n_iterations/2) %/% length(batches)),
                     batches[1:(round(n_iterations/2) %% length(batches))])
      # Set balanced sample length
      n_sampled <- min(sample_max,
                       max(1, min(round(downsampling_rate*(min(length(cluster1_cells), length(cluster2_cells))/2)),
                                  floor(min(cluster1_batch_freq[batches], cluster2_batch_freq[batches])/2))))
    }
  } else {
    use_batch <- NULL
    # Set balanced sample length
    n_sampled <- min(sample_max,
                     max(1, floor(downsampling_rate*(min(length(cluster1_cells), length(cluster2_cells))/2))))
  }

  # Filter input matrix
  if (proceed == TRUE) {
    # Prepare input matrix
    comparison_input <- input_matrix[c(cluster1_cells, cluster2_cells),]
    comparison_input <- methods::as(comparison_input, "dgCMatrix")
    # Filter to minimum reads per feature
    if (is.null(min_reads)) {
      comparison_input <- comparison_input[, Matrix::colSums(comparison_input) > 0]
    } else {
      feature_sums_cluster1 <- Matrix::colSums(comparison_input[cluster1_cells,])
      feature_sums_cluster2 <- Matrix::colSums(comparison_input[cluster2_cells,])
      min_reads_cluster1 <- length(cluster1_cells)/min_reads
      min_reads_cluster2 <- length(cluster2_cells)/min_reads
      keep_features <- (feature_sums_cluster1 > min_reads_cluster1) | (feature_sums_cluster2 > min_reads_cluster2)
      if (sum(keep_features) == 0) {
        comparison_result <- "merge"
        proceed <- FALSE
        # Add to comparison records
        if (!is.null(comparison_records)) {
          current_comparison <- matrix(rep(NA, ncol(comparison_records)), nrow = 1, ncol = ncol(comparison_records))
          colnames(current_comparison) <- colnames(comparison_records)
          current_comparison <- data.frame(current_comparison)
          current_comparison$comparison <- paste0(cluster1_name, " vs. ", cluster2_name)
          current_comparison$cluster1_size <- length(cluster1_cells)
          current_comparison$cluster2_size <- length(cluster2_cells)
          current_comparison$connectivity <- adjacent
          current_comparison$decision <- "merge: minimum reads"
          if ("root_distance" %in% colnames(comparison_records)) {
            current_comparison$root_distance <- P0_distance
            current_comparison$subtree_distance <- P_i_distance
          }
        }
      } else {
        comparison_input <- comparison_input[, keep_features]
      }
    }
  }

  # Go ahead with comparison
  if (proceed == TRUE) {
    # Run random forest classifiers
    rf_comparison_list <- parallel::mclapply(seq_len(round(n_iterations/2)), FUN = function(i) {
      .runRF(i = i,
             cluster1_cells = cluster1_cells,
             cluster2_cells = cluster2_cells,
             cluster1_cell_batches = cluster1_cell_batches,
             cluster2_cell_batches = cluster2_cell_batches,
             n_trees = n_trees,
             max_repeat_errors = max_repeat_errors,
             collect_all_metrics = collect_all_metrics,
             n_sampled = n_sampled,
             input_matrix = comparison_input,
             use_batch = use_batch) },
      mc.cores = n_cores,
      mc.set.seed = TRUE)

    # Get means of comparison statistics
    # Balanced accuracy
    accuracies <- unlist(do.call(rbind, rf_comparison_list)[, "balanced_accuracy"])
    mean_acc <- mean(accuracies)
    var_acc <- stats::var(accuracies)
    mean_err <- (1 - mean_acc)*n_sampled

    # Permutation accuracies
    permutation_accuracies <- unlist(do.call(rbind, rf_comparison_list)[, "permutation_balanced_accuracy"])
    mean_permutation_acc <- mean(permutation_accuracies)
    var_permutation_acc = stats::var(permutation_accuracies)

    # Percentile of permuted mean via kernel density estimator
    percentile_acc <- 1 - spatstat.univar::CDF(stats::density(permutation_accuracies))(mean_acc)
    # Percentile of permuted variance via bootstrap
    boot_permutation_var <- apply(matrix(seq(1:n_iterations), 1, n_iterations), 2,
                                  function(x) stats::var(sample(permutation_accuracies, n_iterations, replace = TRUE)))
    percentile_var <- spatstat.univar::CDF(stats::density(boot_permutation_var))(var_acc)

    # Repeated errors
    if (max_repeat_errors > 0 | collect_all_metrics == TRUE) {
      # Cluster 1
      error_indices1 <- unlist(do.call(rbind, rf_comparison_list)[, "error_indices1"])
      if (!length(error_indices1)) {
        n_repeat_err1 <- 0
        mean_repeat_err1 <- 0
      } else {
        error_freqs1 <- table(error_indices1)
        n_repeat_err1 <- sum(error_freqs1 > max(n_iterations*0.3, 1))
        mean_repeat_err1 <- (sum(error_freqs1[error_freqs1 > max(n_iterations*0.3, 1)])/n_iterations)*(n_sampled/length(cluster1_cells))
      }
      # Cluster 2
      error_indices2 <- unlist(do.call(rbind, rf_comparison_list)[, "error_indices2"])
      if (!length(error_indices2)) {
        n_repeat_err2 <- 0
        mean_repeat_err2 <- 0
      } else {
        error_freqs2 <- table(error_indices2)
        n_repeat_err2 <- sum(error_freqs2 > max(n_iterations*0.3, 1))
        mean_repeat_err2 <- (sum(error_freqs2[error_freqs2 > max(n_iterations*0.3, 1)])/n_iterations)*(n_sampled/length(cluster2_cells))
      }

      # Modified accuracy excluding repeated errors
      if (max_repeat_errors > 0) {
        mean_modified_acc <- (n_sampled - (mean_err - ((mean_repeat_err1 + mean_repeat_err2)/2)))/n_sampled
        if (mean_modified_acc > 1) {
          mean_modified_acc <- 1
        }
        modified_accuracies <- accuracies + (mean_modified_acc - mean_acc)
        var_modified_acc <- stats::var(modified_accuracies)
        percentile_modified_acc <- 1 - spatstat.univar::CDF(stats::density(permutation_accuracies))(mean_modified_acc)
        percentile_modified_var <- spatstat.univar::CDF(stats::density(boot_permutation_var))(var_modified_acc)
      }
    }
    # Feature importance
    if (collect_all_metrics == TRUE) {
      feature_importance <- do.call(rbind, rf_comparison_list)[, "feature_importance"]
      feature_importance <- do.call(rbind, feature_importance)
      mean_feature_importance <- apply(feature_importance, 2, function(x) mean(x, na.rm = TRUE))
    }

    # Get batch values
    batch_acc <- c()
    batch_var <- c()
    if (!is.null(use_batch)) {
      for (b in 1:length(batches)) {
        batch_inds <- which(use_batch == batches[b])
        batch_acc <- c(batch_acc, mean(accuracies[batch_inds]))
        batch_var <- c(batch_var, stats::var(accuracies[batch_inds]))
      }
    }

    # Create record
    if (!is.null(comparison_records)) {
      current_comparison <- matrix(rep(NA, ncol(comparison_records)), nrow = 1, ncol = ncol(comparison_records))
      colnames(current_comparison) <- colnames(comparison_records)
      current_comparison <- data.frame(current_comparison)
      current_comparison <- dplyr::mutate(current_comparison,
                                          comparison = paste0(cluster1_name, " vs. ", cluster2_name),
                                          cluster1_size = length(cluster1_cells),
                                          cluster2_size = length(cluster2_cells),
                                          sample_size = n_sampled,
                                          mean_accuracy = mean_acc,
                                          var_accuracy = var_acc,
                                          mean_errors = mean_err,
                                          mean_permuted_accuracy = mean_permutation_acc,
                                          var_permuted_accuracy = var_permutation_acc,
                                          percentile_accuracy = percentile_acc,
                                          percentile_variance = percentile_var)
    } else {
      current_comparison <- data.frame(comparison = paste0(cluster1_name, " vs. ", cluster2_name),
                                       cluster1_size = length(cluster1_cells),
                                       cluster2_size = length(cluster2_cells),
                                       sample_size = n_sampled,
                                       mean_accuracy = mean_acc,
                                       var_accuracy = var_acc,
                                       mean_errors = mean_err,
                                       mean_permuted_accuracy = mean_permutation_acc,
                                       var_permuted_accuracy = var_permutation_acc,
                                       percentile_accuracy = percentile_acc,
                                       percentile_variance = percentile_var)
    }
    if (max_repeat_errors > 0 | collect_all_metrics == TRUE) {
      current_comparison <- dplyr::mutate(current_comparison,
                                          n_repeat_errors1 = n_repeat_err1,
                                          n_repeat_errors2 = n_repeat_err2,
                                          mean_repeat_errors1 = mean_repeat_err1,
                                          mean_repeat_errors2 = mean_repeat_err2)
      if (max_repeat_errors > 0) {
        current_comparison <- dplyr::mutate(current_comparison,
                                            mean_modified_accuracy = mean_modified_acc,
                                            var_modified_accuracy = var_modified_acc,
                                            percentile_modified_accuracy = percentile_modified_acc,
                                            percentile_modified_variance = percentile_modified_var)
      }
    }
    if (!is.null(use_batch)) {
      current_comparison <- dplyr::mutate(current_comparison,
                                          batches_used = paste(batches, collapse = "; "),
                                          batch_mean_accuracies = paste(round(batch_acc, 5), collapse = "; "),
                                          batch_var_accuracies = paste(round(batch_var, 5), collapse = "; "))
    }
    if (min_connections > 0 | collect_all_metrics == TRUE) {
      current_comparison <- dplyr::mutate(current_comparison,
                                          connectivity = adjacent)
    }
    if ("root_distance" %in% colnames(comparison_records)) {
      current_comparison <- dplyr::mutate(current_comparison,
                                          root_distance = P0_distance,
                                          subtree_distance = P_i_distance)
    }

    if (collect_all_metrics == TRUE) {
      # Add to feature importance dataframe
      mean_feature_importance <- data.frame(matrix(data = mean_feature_importance,
                                                   nrow = 1,
                                                   ncol = length(mean_feature_importance)))
      mean_feature_importance <- cbind(data.frame("cluster1" = cluster1_name,
                                                  "cluster2" = cluster2_name),
                                       mean_feature_importance)
      colnames(mean_feature_importance) <- c("cluster1", "cluster2",
                                             colnames(comparison_input))
      feature_importance_records <- plyr::rbind.fill(list(feature_importance_records,
                                                          mean_feature_importance))
    }

    # Get comparison result based on conditions
    if (all(current_comparison$mean_accuracy[1] >= min_accuracy,
            current_comparison$percentile_accuracy[1] < alpha,
            (current_comparison$percentile_variance[1] < alpha |
             use_variance == FALSE))) {
      comparison_result <- "split"
      reason <- "split"
      if (use_variance == FALSE) {
        max_p <- current_comparison$percentile_accuracy[1]
      } else {
        max_p <- max(current_comparison$percentile_accuracy[1],
                     current_comparison$percentile_variance[1])
      }
    } else if (max_repeat_errors > 0) {
      if (all(current_comparison$mean_accuracy[1] >= min_accuracy,
              max_repeat_errors > 0,
              current_comparison$percentile_accuracy[1] < alpha*2,
              (current_comparison$percentile_variance[1] < alpha*2 |
               use_variance == FALSE),
              current_comparison$percentile_modified_accuracy[1] < alpha,
              (current_comparison$percentile_modified_variance[1] < alpha |
               use_variance == FALSE),
              current_comparison$n_repeat_errors1[1] <= min(max_repeat_errors, 0.3*length(cluster1_cells)),
              current_comparison$n_repeat_errors2[1] <= min(max_repeat_errors, 0.3*length(cluster2_cells)))) {
        comparison_result <- "split"
        reason <- "split: repeat error"
        if (use_variance == FALSE) {
          max_p <- min(current_comparison$percentile_accuracy[1],
                       current_comparison$percentile_modified_accuracy[1])
        } else {
          max_p <- max(min(current_comparison$percentile_accuracy[1],
                           current_comparison$percentile_modified_accuracy[1]),
                       min(current_comparison$percentile_variance[1],
                           current_comparison$percentile_modified_variance[1]))
        }
      } else if (all(current_comparison$mean_accuracy[1] < min_accuracy,
                     current_comparison$percentile_accuracy[1] < alpha,
                     (current_comparison$percentile_variance[1] < alpha |
                      use_variance == FALSE))) {
        comparison_result <- "merge"
        reason <- "merge: min accuracy"
      } else {
        comparison_result <- "merge"
        reason <- "merge"
      }
    } else {
      comparison_result <- "merge"
      reason <- "merge"
    }

    # Add decision
    if (!is.null(comparison_records)) {
      current_comparison$decision <- reason
    }
  }

  # Time
  if (!is.null(comparison_start_time)) {
    current_comparison$time <- round(difftime(Sys.time(), comparison_start_time, units = "secs"), 2)
  }
  # Append new results to previous records
  if (!is.null(comparison_records)) {
    if (!identical(colnames(comparison_records), colnames(current_comparison))) {
      current_comparison <- current_comparison[, colnames(comparison_records)]
    }
    comparison_records <- rbind(comparison_records, current_comparison)
  }
  output_list <- list("result" = comparison_result,
                      "comparison_records" = comparison_records,
                      "feature_importance_records" = feature_importance_records,
                      "max_p" = max_p)
  return(output_list)
}

# Run a random forest classifier iteration ---------------------------
#
# This function runs random forest classifiers on train/test splits of the data
# with and without permuted cluster labels.
#
# i -- current iteration
# cluster1_cells -- A vector indicating the cells belonging to the first cluster
# cluster2_cells -- A vector indicating the cells belonging to the second cluster
# cluster1_cell_batches -- A vector indicating the batch IDs for cluster 1
# cluster2_cell_batches -- A vector indicating the batch IDs for cluster 2
# n_trees -- A numeric value indicating the number of trees in each random forest
# max_repeat_errors -- A numeric value indicating the maximum number of cells that will be considered as repeated errors
# collect_all_metrics -- A boolean value indicating whether to collect and save additional metrics
# n_sampled -- A numeric indicating number of cells to sample for random forest
# input_matrix -- A matrix of data, on which the random forest classifier will be trained/tested
#
# Returns a list containing the following elements:
# balanced_accuracy -- The prediction accuracy score for the true cluster labels
# permutation_balanced_accuracy -- The prediction accuracy score for the permuted cluster labels
# error_indices1 -- The names of cluster 1 cells that were misassigned
# error_indices2 -- The names of cluster 1 cells that were misassigned
# feature_importance -- A named vector of feature importance scores
.runRF <- function(i = NULL,
                   cluster1_cells,
                   cluster2_cells,
                   cluster1_cell_batches = NULL,
                   cluster2_cell_batches = NULL,
                   n_trees,
                   max_repeat_errors,
                   collect_all_metrics,
                   n_sampled,
                   input_matrix,
                   use_batch = NULL) {
  # Batch-dependent
  if (!is.null(use_batch)) {
    use_batch_i <- use_batch[i]
    cluster1_cells <- cluster1_cells[cluster1_cell_batches == use_batch_i]
    cluster2_cells <- cluster2_cells[cluster2_cell_batches == use_batch_i]
  }

  # Downsample clusters
  cluster2_cells <- sample(cluster2_cells, n_sampled*2)
  cluster1_cells <- sample(cluster1_cells, n_sampled*2)
  # Set true cluster labels
  true_clusters <- c(rep(1, n_sampled*2), rep(2, n_sampled*2))
  # Randomize train/test balanced by class
  case_weights <- c(sample(c(rep(1, n_sampled),
                             rep(0, n_sampled))),
                    sample(c(rep(1, n_sampled),
                             rep(0, n_sampled))))
  # Run random forest classifiers with train & test, then reverse
  if (collect_all_metrics == TRUE) {
    importance <- "impurity"
  } else {
    importance <- "none"
  }
  RF1 <- ranger::ranger(x = input_matrix[c(cluster1_cells, cluster2_cells),],
                        y = as.factor(true_clusters),
                        num.trees = n_trees,
                        case.weights = case_weights,
                        num.threads = 1,
                        holdout = TRUE,
                        importance = importance,
                        verbose = FALSE)
  RF2 <- ranger::ranger(x = input_matrix[c(cluster1_cells, cluster2_cells),],
                        y = as.factor(true_clusters),
                        num.trees = n_trees,
                        case.weights = 1 - case_weights,
                        num.threads = 1,
                        holdout = TRUE,
                        importance = importance,
                        verbose = FALSE)

  # Extract performance
  balanced_accuracy <- c(1- RF1$prediction.error, 1 - RF2$prediction.error)

  # Indices of errors
  if (max_repeat_errors > 0 | collect_all_metrics == TRUE) {
    error_indices1 = c(which(RF1$predictions[1:length(cluster1_cells)] != 1),
                       which(RF2$predictions[1:length(cluster1_cells)] != 1))
    error_indices2 = c(which(RF1$predictions[(length(cluster1_cells) + 1):length(true_clusters)] != 2),
                       which(RF2$predictions[(length(cluster1_cells) + 1):length(true_clusters)] != 2))
  }

  # Feature importance
  if (collect_all_metrics == TRUE) {
    feature_importance <- RF1$variable.importance + RF2$variable.importance
  }

  # Permuted cluster labels
  permuted_clusters <- sample(true_clusters)
  # Run random forest classifiers with train & test, then reverse
  RF1_permuted <- ranger::ranger(x = input_matrix[c(cluster1_cells, cluster2_cells),],
                                 y = as.factor(permuted_clusters),
                                 num.trees = n_trees,
                                 case.weights = case_weights,
                                 num.threads = 1,
                                 holdout = TRUE,
                                 verbose = FALSE)
  RF2_permuted <- ranger::ranger(x = input_matrix[c(cluster1_cells, cluster2_cells),],
                                 y = as.factor(permuted_clusters),
                                 num.trees = n_trees,
                                 case.weights = 1 - case_weights,
                                 num.threads = 1,
                                 holdout = TRUE,
                                 verbose = FALSE)
  # Extract performance
  permutation_balanced_accuracy <- c(1- RF1_permuted$prediction.error,
                                     1 - RF2_permuted$prediction.error)

  # Compile output
  if (collect_all_metrics == TRUE) {
    output <- list("balanced_accuracy" = balanced_accuracy,
                   "permutation_balanced_accuracy" = permutation_balanced_accuracy,
                   "error_indices1" = error_indices1,
                   "error_indices2" = error_indices2,
                   "feature_importance" = feature_importance)
  } else if (max_repeat_errors > 0) {
    output <- list("balanced_accuracy" = balanced_accuracy,
                   "permutation_balanced_accuracy" = permutation_balanced_accuracy,
                   "error_indices1" = error_indices1,
                   "error_indices2" = error_indices2)
  } else {
    output <- list("balanced_accuracy" = balanced_accuracy,
                   "permutation_balanced_accuracy" = permutation_balanced_accuracy)
  }
  return(output)
}

# Check comparison records ---------------------------
#
# Check whether two clusters have already been previously compared and extract
# result.
#
# cluster1_name -- A string indicating the name of the first cluster
# cluster1_cells -- A vector indicating the cells belonging to the first cluster
# cluster2_name -- A string indicating the name of the second cluster
# cluster2_cells -- A vector indicating the cells belonging to the second cluster
# comparison_records -- A dataframe containing the comparison records output of the previous runs of this function
# type -- String ("default" or "bridge") indicating to whether to disregard distance-based records
.checkComparisonRecords <- function(cluster1_name,
                                    cluster1_cells,
                                    cluster2_name,
                                    cluster2_cells,
                                    comparison_records,
                                    type = "default") {
  # Balanced sample length
  n_sampled <- min(length(cluster1_cells), length(cluster2_cells))/2
  # Whether to flip cluster1/cluster2 (only relevant if comparison has already been calculated)
  swap <- FALSE
  comparison_name1 <- paste0(cluster1_name, " vs. ", cluster2_name)
  comparison_name2 <- paste0(cluster2_name, " vs. ", cluster1_name)
  # Check if cluster1/cluster2 comparison has already been calculated:
  if (comparison_name1 %in% comparison_records$comparison) {
    current_comparison <- dplyr::filter(comparison_records,comparison == comparison_name1)
    previously_compared <- TRUE
  } else if (comparison_name2 %in% comparison_records$comparison) {
    current_comparison <- dplyr::filter(comparison_records, comparison == comparison_name2)
    previously_compared <- TRUE
    swap <- TRUE
  } else {
    previously_compared <- FALSE
  }
  # Conditions
  if (previously_compared == TRUE) {
    if (grepl("merge", current_comparison$decision)) {
      comparison_result <- "merge"
    } else if (type == "bridge" & current_comparison$decision == "split: distance") {
      previously_compared <- FALSE
      comparison_result <- NULL
    } else {
      comparison_result <- "split"
    }
  } else {
    comparison_result <- NULL
  }
  # Return
  output_list <- list("result" = comparison_result,
                      "previously_compared" = previously_compared)
  return(output_list)
}

# Check cluster distance ---------------------------
#
# Calculate the average linkage distance between two clusters using a
# distance matrix calculated from the dimensional reduction.
#
# object -- An object of type Seurat, SingleCellExperiment, or ArchRProject
# key -- A character string indicating the name under which data is stored for this run
# cluster1_name -- A string indicating the name of the first cluster
# cluster1_cells -- A vector indicating the cells belonging to the first cluster
# cluster2_name -- A string indicating the name of the second cluster
# cluster2_cells -- A vector indicating the cells belonging to the second cluster
# distance_awareness -- A multiplier used to set the distance threshold
# distance_approx -- Whether to use approximate distance
# use_input_matrix -- Indicates the matrix to use for the current tree/subtree
# dist_matrix -- A cell to cell distance matrix
# reduction -- A dimensionality reduction
# distance_records -- A dataframe containing the distance records output from this run
.checkDistance <- function(object,
                           key,
                           cluster1_name,
                           cluster1_cells,
                           cluster2_name,
                           cluster2_cells,
                           distance_awareness,
                           distance_approx,
                           use_input_matrix,
                           dist_matrix,
                           reduction,
                           distance_records) {
  if (distance_approx == FALSE) {
    .requirePackage("clv", source = "cran")
  }
  # Default output
  distance_conflict <- FALSE
  P0_distance <- NA
  P_i_distance <- NA
  # Don't proceed if distance_awareness is FALSE
  if (methods::is(distance_awareness, "numeric")) {
    # Distance by centroid or average linkage
    if (distance_approx == TRUE) {
      # Calculate P0 distance
      if (is.null(reduction)) {
        current_reduction <- as.matrix(.retrieveData(object, key, "reduction", "P0_reduction"))[c(cluster1_cells, cluster2_cells), ]
      } else {
        current_reduction <- reduction[c(cluster1_cells, cluster2_cells), ]
      }
      P0_distance <- .getCentroidDistance(reduction = current_reduction,
                                          clusters = c(rep(1, length(cluster1_cells)),
                                                       rep(2, length(cluster2_cells))))[1,2]
      # Calculate subtree distance if applicable
      if (use_input_matrix != "P0" & is.null(reduction)) {
        current_reduction <- as.matrix(.retrieveData(object, key, "reduction", paste0(use_input_matrix, "_reduction")))[c(cluster1_cells, cluster2_cells), ]
        P_i_distance <- .getCentroidDistance(reduction = current_reduction,
                                             clusters = c(rep(1, length(cluster1_cells)),
                                                          rep(2, length(cluster2_cells))))[1,2]
      }
      # Clean up
      rm(current_reduction)
    } else {
      # Calculate P0 distance
      if (is.null(dist_matrix)) {
        current_dist_matrix <- as.matrix(.retrieveData(object, key, "reduction", "P0_reduction_dist"))[c(cluster1_cells, cluster2_cells), c(cluster1_cells, cluster2_cells)]
      } else {
        current_dist_matrix <- dist_matrix[c(cluster1_cells, cluster2_cells), c(cluster1_cells, cluster2_cells)]
      }
      P0_distance <- clv::cls.scatt.diss.mx(diss.mx = current_dist_matrix,
                                            clust = as.integer(c(rep(1, length(cluster1_cells)),
                                                                 rep(2, length(cluster2_cells)))))$intercls.average["c1", "c2"]
      # Calculate subtree distance if applicable
      if (use_input_matrix != "P0" & is.null(dist_matrix)) {
        current_dist_matrix <- as.matrix(.retrieveData(object, key, "reduction", paste0(use_input_matrix, "_reduction_dist")))[c(cluster1_cells, cluster2_cells), c(cluster1_cells, cluster2_cells)]
        P_i_distance <- clv::cls.scatt.diss.mx(diss.mx = current_dist_matrix,
                                               clust = as.integer(c(rep(1, length(cluster1_cells)),
                                                                    rep(2, length(cluster2_cells)))))$intercls.average["c1", "c2"]
      }
      # Clean up
      rm(current_dist_matrix)
    }

    # Compare to previous records
    if (!is.null(distance_records)) {
      if (cluster1_name %in% distance_records$cluster_name & cluster2_name %in% distance_records$cluster_name) {
        if (use_input_matrix == "P0") {
          previous_P0_distance <- max(dplyr::filter(distance_records,
                                                    cluster_name == cluster1_name |
                                                      cluster_name == cluster2_name)$min_root_distance, na.rm = TRUE)
          if (P0_distance > (previous_P0_distance*distance_awareness)) {
            distance_conflict <- TRUE
          }
        } else if (use_input_matrix != "P0" & !is.na(P_i_distance)) {
          previous_P_i_distance <- max(dplyr::filter(distance_records,
                                                     cluster_name == cluster1_name |
                                                       cluster_name == cluster2_name)$min_subtree_distance, na.rm = TRUE)
          if (P_i_distance > (previous_P_i_distance*distance_awareness)) {
            distance_conflict <- TRUE
          }
        }
      }
    }
  }
  return(list("distance_conflict" = distance_conflict,
              "P0_distance" = P0_distance,
              "P_i_distance" = P_i_distance))
}

# Add cluster distance to records ---------------------------
#
# If cluster distance is lower than previous minimum cluster difference, add
# to distance records.
#
# cluster1_name -- A string indicating the name of the first cluster
# cluster2_name -- A string indicating the name of the second cluster
# P0_distance -- Distance between two clusters using root tree dimensionality reduction
# P_i_distance -- Distance between two clusters using subtree dimensionality reduction, if applicable
# max_p -- Max p-value for comparison (for filtering if adjusted alpha threshold decreases)
# distance_records -- A dataframe containing the distance records output from this run
.addDistance <- function(cluster1_name,
                         cluster2_name,
                         P0_distance,
                         P_i_distance,
                         max_p,
                         distance_records) {
  # Cluster 1
  if (cluster1_name %in% rownames(distance_records)) {
    distance_records[cluster1_name, "min_root_distance"] <- min(distance_records[cluster1_name, "min_root_distance"], P0_distance, na.rm = TRUE)
    distance_records[cluster1_name, "min_subtree_distance"] <- suppressWarnings(min(distance_records[cluster1_name, "min_subtree_distance"], P_i_distance, na.rm = TRUE))
  } else {
    distance_records <- rbind(distance_records, data.frame(cluster_name = cluster1_name,
                                                           min_root_distance = P0_distance,
                                                           min_subtree_distance = P_i_distance,
                                                           max_pval = max_p))
  }
  # Cluster 2
  if (cluster2_name %in% rownames(distance_records)) {
    distance_records[cluster2_name, "min_root_distance"] <- min(distance_records[cluster2_name, "min_root_distance"], P0_distance, na.rm = TRUE)
    distance_records[cluster2_name, "min_subtree_distance"] <- suppressWarnings(min(distance_records[cluster2_name, "min_subtree_distance"], P_i_distance, na.rm = TRUE))
  } else {
    distance_records <- rbind(distance_records, data.frame(cluster_name = cluster2_name,
                                                           min_root_distance = P0_distance,
                                                           min_subtree_distance = P_i_distance,
                                                           max_pval = max_p))
  }
  rownames(distance_records) <- distance_records$cluster_name
  return(distance_records)
}
