# ---------------------------------------------------------------------------
# General helper functions
# ---------------------------------------------------------------------------

#' Retrieve all CHOIR metadata
#'
#' Retrieve all stored CHOIR data from object
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{ArchRProject} that has undergone CHOIR clustering.
#' @param key The name under which CHOIR-related data for this run is stored in
#' the object. Defaults to “CHOIR”.
#'
#' @return Returns the data stored under the provided key.
#'
#' @export
#'
getRecords <- function(object,
                       key = "CHOIR") {
  # By object type
  if (methods::is(object, "Seurat")) {
    output_data <- object@misc[[key]]
  } else if (methods::is(object, "SingleCellExperiment")) {
    output_data <- object@metadata[[key]]
  } else if (methods::is(object, "ArchRProject")) {
    output_data <- object@projectMetadata[[key]]
  }
  return(output_data)
}

# Retrieve cell IDs ---------------------------
#
# Extract cell IDs/names from provided object
#
# object -- An object of class Seurat, SingleCellExperiment, or ArchRProject
# use_assay -- For Seurat objects, character string/vector indicating assay to use
.getCellIDs <- function(object,
                        use_assay = NULL) {
  # By object type
  if (methods::is(object, "Seurat")) {
    if (is.null(use_assay)) {
      use_assay <- Seurat::DefaultAssay(object)
    }
    if (length(use_assay) > 1) {
      # If multiple assays, check that cell IDs are identical
      cell_IDs <- colnames(object[[use_assay[1]]])
      for (i in 2:length(use_assay)) {
        cell_IDs_i <- colnames(object[[use_assay[i]]])
        if (!identical(cell_IDs, cell_IDs_i)) {
          stop("Cell IDs do not match across provided assays indicated by parameter 'use_assay'. Please supply valid input!")
        }
      }
    } else {
      cell_IDs <- colnames(object[[use_assay]])
    }
  } else if (methods::is(object, "SingleCellExperiment")) {
    cell_IDs <- rownames(object@colData)
  } else if (methods::is(object, "ArchRProject")) {
    cell_IDs <- rownames(object@cellColData)
  }
  return(cell_IDs)
}

# Retrieve matrix ---------------------------
#
# Extract a matrix from provided object
#
# object -- An object of class Seurat, SingleCellExperiment, or ArchRProject
# use_matrix -- If there is a user-supplied matrix, do not retrieve a matrix from the object
# use_assay -- For Seurat or SingleCellExperiment objects, a character string indicating the assay to use
# use_slot -- For Seurat objects, a character string indicating the slot/layer to use
# ArchR_matrix -- For ArchR objects, a character string indicating which matrix to use
# use_features -- A vector of feature names to use to subset the matrix
# exclude_features -- A vector of feature names to exclude from the matrix
# use_cells -- A vector of cell IDs to use to subset the matrix
# verbose -- A Boolean value indicating whether to use verbose output during the execution of this function
.getMatrix <- function(object = NULL,
                       use_matrix = NULL,
                       use_assay = NULL,
                       use_slot = NULL,
                       ArchR_matrix = NULL,
                       use_features = NULL,
                       exclude_features = NULL,
                       use_cells = NULL,
                       verbose) {
  # By object type
  if (is.null(object) | methods::is(object, "Seurat")) {
    use_matrix <- .getMatrix.Seurat(object,
                                    use_matrix,
                                    use_assay,
                                    use_slot,
                                    use_features,
                                    exclude_features,
                                    use_cells,
                                    verbose)
  } else if (methods::is(object, "SingleCellExperiment")) {
    use_matrix <- .getMatrix.SingleCellExperiment(object,
                                                  use_matrix,
                                                  use_assay,
                                                  use_slot,
                                                  use_features,
                                                  exclude_features,
                                                  use_cells,
                                                  verbose)
  } else if (methods::is(object, "ArchRProject")) {
    use_matrix <- .getMatrix.ArchR(object,
                                   use_matrix,
                                   ArchR_matrix,
                                   use_features,
                                   exclude_features,
                                   use_cells,
                                   verbose)
  }
  # Return matrix
  return(use_matrix)
}

.getMatrix.Seurat <- function(object,
                              use_matrix,
                              use_assay,
                              use_slot,
                              use_features,
                              exclude_features,
                              use_cells,
                              verbose) {
  # If matrix is not provided as input
  if (is.null(use_matrix)) {
    # Get assay
    if (is.null(use_assay)) {
      use_assay <- Seurat::DefaultAssay(object)
    } else {
      # Check that input assay is present in object
      .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
    }
    # Determine which slot to use
    if (is.null(use_slot)) {
      if (use_assay %in% c("RNA", "sketch")) {
        use_slot <- "data"
      } else if (use_assay == "SCT" | use_assay == "integrated") {
        use_slot <- "scale.data"
      } else {
        stop("When using a non-standard assay in a Seurat object, please supply a valid input for the slot parameter.")
      }
      # Check that selected slot is present within selected assay in object
      .validInput(use_slot, "use_slot", list(object, use_assay, FALSE, NULL))
    }
    # Extract matrix
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Preparing matrix using '", use_assay, "' assay and '", use_slot, "' slot..")
    if ("Assay5" %in% methods::is(object[[use_assay]])) {
      use_matrix <- object[[use_assay]]@layers[[use_slot]]
      colnames(use_matrix) <- colnames(object[[use_assay]])
      rownames(use_matrix) <- rownames(object[[use_assay]])
    } else {
      use_matrix <- methods::slot(object[[use_assay]], name = use_slot)
    }
  } else {
    # If use_assay is not NULL
    if (!is.null(use_assay)) {
      if (verbose) warning("Input for parameter 'use_assay' is not used when a matrix is provided for parameter 'use_matrix'.")
    }
    # If use_slot is not NULL
    if (!is.null(use_slot)) {
      if (verbose) warning("Input for parameter 'use_slot' is not used when a matrix is provided for parameter 'use_matrix'.")
    }
  }

  # If matrix has no row names
  if (is.null(rownames(use_matrix))) {
    # Stop if trying to subset features
    if (!is.null(use_features) | !is.null(exclude_features)) {
      stop("Provided 'use_matrix' has no row names, therefore, input for parameters 'use_features' and 'exclude_features' cannot be used.")
    }
    rownames(use_matrix) <- seq(1, nrow(use_matrix))
  }
  # Subset matrix by selected features
  if (!is.null(use_features)) {
    unused_features <- use_features[!(use_features %in% rownames(use_matrix))]
    if (verbose & (length(unused_features) > 0)) {
      warning("Could not find the following ",
              length(unused_features),
              " features provided by 'use_features' in 'use_matrix': \n",
              unused_features)
    }
  } else {
    use_features <- rownames(use_matrix)
  }
  use_features <- use_features[!(use_features %in% exclude_features)]
  if (length(use_features) < 1) {
    stop("No remaining features in matrix. Please check input to 'use_features' and/or 'exclude_features'!")
  }

  # If matrix has no column names
  if (is.null(colnames(use_matrix))) {
    # Stop if trying to subset cells
    if (!is.null(use_cells)) {
      stop("Provided 'use_matrix' has no column names, therefore, input for parameter 'use_cells' cannot be used.")
    }
    colnames(use_matrix) <- seq(1, ncol(use_matrix))
  }
  # Subset matrix by selected cells
  if (!is.null(use_cells)) {
    unused_cells <- use_cells[!(use_cells %in% colnames(use_matrix))]
    if (verbose & (length(unused_cells) > 0)) {
      warning("Could not find the following ",
              length(unused_cells),
              " cells provided by 'use_cells' in 'use_matrix': \n",
              unused_cells)
    }
  } else {
    use_cells <- colnames(use_matrix)
  }

  use_matrix <- use_matrix[use_features, use_cells]
  return(use_matrix)
}

.getMatrix.SingleCellExperiment <- function(object,
                                            use_matrix,
                                            use_assay,
                                            use_slot,
                                            use_features,
                                            exclude_features,
                                            use_cells,
                                            verbose) {
  # If matrix is not provided as input
  if (is.null(use_matrix)) {
    # Get assay
    if (is.null(use_assay)) {
      use_assay <- "logcounts"
      .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
    }
    use_matrix <- object@assays@data[[use_assay]]
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Preparing input matrix using '", use_assay, "' assay..")
  } else {
    # If assay is not NULL
    if (!is.null(use_assay)) {
      if (verbose) warning("Input for parameter 'use_assay' is not used when a matrix is provided for parameter 'use_matrix'.")
    }
  }
  # If use_slot is not NULL
  if (!is.null(use_slot)) {
    if (verbose) warning("Input for parameter 'use_slot' is not used when input for parameter 'object' is of type 'SingleCellExperiment'.")
  }

  # If matrix has no row names
  if (is.null(rownames(use_matrix))) {
    # Stop if trying to subset features
    if (!is.null(use_features) | !is.null(exclude_features)) {
      stop("Provided 'use_matrix' has no row names, therefore, input for parameters 'use_features' and 'exclude_features' cannot be used.")
    }
    rownames(use_matrix) <- seq(1, nrow(use_matrix))
  }
  # Subset matrix by selected features
  if (!is.null(use_features)) {
    unused_features <- use_features[!(use_features %in% rownames(use_matrix))]
    if (verbose & (length(unused_features) > 0)) {
      warning("Could not find the following ",
              length(unused_features),
              " features provided by 'use_features' in 'use_matrix': \n",
              unused_features)
    }
  } else {
    use_features <- rownames(use_matrix)
  }
  use_features <- use_features[!(use_features %in% exclude_features)]
  if (length(use_features) < 1) {
    stop("No remaining features in matrix. Please check input to 'use_features' and/or 'exclude_features'!")
  }

  # If matrix has no column names
  if (is.null(colnames(use_matrix))) {
    # Stop if trying to subset cells
    if (!is.null(use_cells)) {
      stop("Provided 'use_matrix' has no column names, therefore, input for parameter 'use_cells' cannot be used.")
    }
    colnames(use_matrix) <- seq(1, ncol(use_matrix))
  }
  # Subset matrix by selected cells
  if (!is.null(use_cells)) {
    unused_cells <- use_cells[!(use_cells %in% colnames(use_matrix))]
    if (verbose & (length(unused_cells) > 0)) {
      warning("Could not find the following ",
              length(unused_cells),
              " cells provided by 'use_cells' in 'use_matrix': \n",
              unused_cells)
    }
  } else {
    use_cells <- colnames(use_matrix)
  }

  use_matrix <- use_matrix[use_features, use_cells]
  return(use_matrix)
}

.getMatrix.ArchR <- function(object,
                             use_matrix,
                             ArchR_matrix,
                             use_features,
                             exclude_features,
                             use_cells,
                             verbose) {
  # If matrix is not provided as input
  if (is.null(use_matrix)) {
    # Matrix to pull from
    if (is.null(ArchR_matrix)) {
      ArchR_matrix <- "GeneScoreMatrix"
    }
    if (ArchR_matrix == "GeneScoreMatrix") {
      gene_score_matrix <- suppressMessages(ArchR::getMatrixFromProject(object, useMatrix = "GeneScoreMatrix"))
      feature_names <- gene_score_matrix@elementMetadata$name
      use_matrix <- gene_score_matrix@assays@data$GeneScoreMatrix
      rownames(use_matrix) <- feature_names
      # Subset to current cells
      if (!is.null(use_cells)) {
        use_matrix <- use_matrix[,use_cells]
      }
      # Subset to current features
      if (!is.null(use_features)) {
        use_matrix <- use_matrix[use_features,]
      }
    } else {
      # Cells to use
      if (is.null(use_cells)) {
        use_cells <- rownames(object@cellColData)
      }
      # Features to use
      if (is.null(use_features)) {
        stop("Please supply input to 'use_features' when extracting ArchR matrix.")
      }
      if (!is.null(exclude_features)) {
        use_features <- dplyr::anti_join(data.frame(use_features), data.frame(exclude_features))
      }
      # Extract matrix & subset
      use_matrix <- suppressMessages(ArchR:::.getPartialMatrix(ArrowFiles = object@sampleColData$ArrowFiles,
                                                               featureDF = use_features,
                                                               cellNames = use_cells,
                                                               useMatrix = ArchR_matrix))
      # Set feature names to chromosome + start
      rownames(use_matrix) <- paste(use_features$seqnames, use_features$start, "_")
    }
  }
  return(use_matrix)
}

# Store matrix ---------------------------
#
# Store a matrix in provided object
#
# object -- An object of class Seurat, SingleCellExperiment, or ArchRProject
# use_matrix -- Matrix to be stored
# use_assay -- For Seurat or SingleCellExperiment objects, a character string indicating the assay to use
# use_slot -- For Seurat objects, a character string indicating the slot/layer to use
# ArchR_matrix -- For ArchR objects, a character string indicating which matrix to use
# verbose -- A Boolean value indicating whether to use verbose output during the execution of this function
.storeMatrix <- function(object,
                         use_matrix,
                         use_assay = NULL,
                         use_slot = NULL,
                         ArchR_matrix = NULL,
                         verbose = TRUE) {
  # By object type
  if (methods::is(object, "Seurat")) {
    object <- .storeMatrix.Seurat(object,
                                  use_matrix,
                                  use_assay,
                                  use_slot,
                                  verbose)
  } else if (methods::is(object, "SingleCellExperiment")) {
    object <- .storeMatrix.SingleCellExperiment(object,
                                                use_matrix,
                                                use_assay,
                                                verbose)
  } else if (methods::is(object, "ArchRProject")) {
    stop("Function '.storeMatrix' does not yet support ArchR objects.")
  }
  # Return object
  return(object)
}

.storeMatrix.Seurat <- function(object,
                                use_matrix,
                                use_assay,
                                use_slot,
                                verbose) {
  # Get assay
  if (is.null(use_assay)) {
    use_assay <- Seurat::DefaultAssay(object)
  } else {
    # Check that input assay is present in object
    .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
  }
  # Determine which slot to use
  if (is.null(use_slot)) {
    stop("For Seurat objects, .storeMatrix requires input for parameter 'use_slot', please supply a valid input for the slot parameter.")
  }
  if (!methods::is(use_slot, "character") | length(use_slot) != 1) {
    stop("Input value for 'use_slot' is not a single value of of class 'character', please supply valid input!")
  }
  # Check that selected slot is NOT already present within selected assay in object
  if ("Assay5" %in% methods::is(object[[use_assay]])) {
    if (use_slot %in% names(object[[use_assay]]@layers)) {
      stop("Layer '", use_slot, "' is already present in assay '", use_assay, "' of provided Seurat v5 object, please supply different input to 'countsplit_suffix' to avoid overwriting data.")
    }
  } else {
    try(slot_exists_1 <- methods::validObject(methods::slot(object[[use_assay]], use_slot)), silent = TRUE)
    if (exists("slot_exists_1")) {
      stop("Slot '", use_slot, "' is already present in assay '", use_assay, "' of provided Seurat object, please supply different input to 'countsplit_suffix' to avoid overwriting data.")
    }
  }
  # Proceed to store matrix
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Storing ", ifelse("Assay5" %in% methods::is(object[[use_assay]]), "layer", "slot"), " '", use_slot, "' under assay '", use_assay, "' in Seurat object.")
  try(object[[use_assay]][use_slot] <- use_matrix)
  # Check that selected slot does now exist within selected assay in object
  if ("Assay5" %in% methods::is(object[[use_assay]])) {
    if (!(use_slot %in% names(object[[use_assay]]@layers))) {
      stop("Layer '", use_slot, "' could not be stored in assay '", use_assay, "' of provided Seurat v5 object.")
    }
  } else {
    try(slot_exists_1 <- methods::validObject(methods::slot(object[[use_assay]], use_slot)), silent = TRUE)
    if (!exists("slot_exists_1")) {
      stop("Slot '", use_slot, "' could not be stored in assay '", use_assay, "' of provided Seurat object, try converting assay to class 'Assay5' before running.")
    } else if (slot_exists_1 == FALSE) {
      stop("Slot '", use_slot, "' could not be stored in assay '", use_assay, "' of provided Seurat object, try converting assay to class 'Assay5' before running.")
    }
  }

  # Return object
  return(object)
}

.storeMatrix.SingleCellExperiment <- function(object,
                                              use_matrix,
                                              use_assay,
                                              verbose) {
  # Check assay
  if (is.null(use_assay)) {
    stop("For SingleCellExperiment objects, .storeMatrix requires input for parameter 'use_assay', please supply a valid input for the slot parameter.")
  } else {
    # Check that input assay is NOT already present in object
    if (use_assay %in% names(object@assays)) {
      stop("Assay '", use_assay, "' provided for parameter '", name, "' is already present in provided SingleCellExperiment object, please supply different input to 'countsplit_suffix' to avoid overwriting data.")
    }
  }
  # Proceed to store matrix
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Storing assay '", use_assay, "' in SingleCellExperiment object.")
  object@assays@data[[use_assay]] <- use_matrix

  # Return object
  return(object)
}


# Store data ---------------------------
#
# Store data in object under specified key
#
# object -- An object of class Seurat, SingleCellExperiment, or ArchRProject
# key -- A string indicating the name under which data is stored for this run
# type -- A string indicating the type of data, so that it can be more easily retrieved
# input_data -- Data to be stored
# name -- Name under which data should be stored
# reduction_method -- For dimensionality reductions, method name is used to store full reductions
# use_assay -- For dimensionality reductions, assay name is used to store full reductions
.storeData <- function(object,
                       key,
                       type,
                       input_data,
                       name,
                       reduction_method = NULL,
                       use_assay = NULL,
                       atac = NULL) {
  # By object type
  if (methods::is(object, "Seurat")) {
    # By type
    if (type == "full_reduction") {
      if (length(use_assay) > 1) {
        for (i in 1:length(use_assay)) {
          reduction_method_i <- .matchArg(reduction_method, i)
          atac_i <- .matchArg(atac, i)
          if (is.null(reduction_method_i)) {
            if (atac_i == FALSE) {
              reduction_method <- "PCA"
            } else {
              reduction_method <- "LSI"
            }
          }
          use_assay_i <- .matchArg(use_assay, i)
          try(object[[paste0(name, "_", use_assay_i)]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = input_data[[i]],
                                                                                                    key = reduction_method_i,
                                                                                                    assay = use_assay_i)), silent = TRUE)
          if (!(paste0(name, "_", use_assay_i) %in% names(object@reductions))) {
            warning("Reduction ", paste0(name, "_", use_assay_i), " is stored only in 'misc' slot under provided CHOIR 'key'.")
          }
        }
      } else {
        try(object[[name]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = input_data,
                                                                        key = reduction_method,
                                                                        assay = use_assay)), silent = TRUE)
        if (!(name %in% names(object@reductions))) {
          warning("Reduction ", name, " is stored only in 'misc' slot under provided CHOIR 'key'.")
        }
      }
    } else if (type == "final_clusters") {
      object@meta.data[, name] <- input_data
    } else {
      object@misc[[key]][[type]][[name]] <- input_data
    }
  } else if (methods::is(object, "SingleCellExperiment")) {
    # By type
    if (type == "full_reduction") {
      if (length(use_assay) > 1) {
        for (i in 1:length(use_assay)) {
          use_assay_i <- .matchArg(use_assay, i)
          input_data_i <- list(input_data[[i]])
          names(input_data_i) <- paste0(name, "_", use_assay_i)
          SingleCellExperiment::reducedDims(object, withDimnames = TRUE) <- input_data_i
        }
      } else {
        input_data <- list(input_data)
        names(input_data) <- name
        SingleCellExperiment::reducedDims(object, withDimnames = TRUE) <- input_data
      }
    } else if (type == "final_clusters") {
      object@colData[, name] <- input_data
    } else {
      object@metadata[[key]][[type]][[name]] <- input_data
    }
  } else if (methods::is(object, "ArchRProject")) {
    # By type
    if (type == "full_reduction") {
      object@reducedDims[[name]] <- input_data
    } else if (type == "final_clusters") {
      object@cellColData[, name] <- input_data
    } else {
      object@projectMetadata[[key]][[type]][[name]] <- input_data
    }
  }
  return(object)
}


# Retrieve data ---------------------------
#
# Retrieve stored data from object
#
# object -- An object of class Seurat, SingleCellExperiment, or ArchRProject
# key -- A string indicating the name under which data is stored for this run
# type -- A string indicating the type of data, so that it can be more easily retrieved
# name -- Name under which data is stored
.retrieveData <- function(object, key, type, name) {
  # By object type
  if (methods::is(object, "Seurat")) {
    # By type
    if (type == "final_clusters" | type == "cell_metadata") {
      output_data <- object@meta.data[, name]
    } else {
      output_data <- object@misc[[key]][[type]][[name]]
    }
  } else if (methods::is(object, "SingleCellExperiment")) {
    # By type
    if (type == "final_clusters" | type == "cell_metadata") {
      output_data <- object@colData[, name]
    } else {
      output_data <- object@metadata[[key]][[type]][[name]]
    }
  } else if (methods::is(object, "ArchRProject")) {
    # By type
    if (type == "final_clusters" | type == "cell_metadata") {
      output_data <- object@cellColData[, name]
    } else {
      output_data <- object@projectMetadata[[key]][[type]][[name]]
    }
  }
  return(output_data)
}

# Match arguments ---------------------------
#
# Match arguments for multiple modalities
#
# var -- Variable
# i -- Index of current modality
.matchArg <- function(var, i) {
  if (length(var) > 1) {
    var_i <- var[i]
  }  else {
    var_i <- var[1]
  }
  return(var_i)
}

# Retrieve parameters ---------------------------
#
# Retrieve parameter values stored in object
#
# input -- Supplied value for parameter
# name -- Name of parameter
# parameter_list -- List of stored parameter values
# default_list -- List of default parameter values
# function_name -- Name of function for which parameter values were stored
.retrieveParam <- function(input,
                           name,
                           parameter_list,
                           default_list,
                           function_name = "buildTree") {
  # If supplied parameter list contains name
  if (name %in% names(parameter_list)) {
    # For any parameters set to NULL, use values from parameter_list
    if (is.null(input)) {
      input <- parameter_list[[name]]
    } else if (input != parameter_list[[name]]) {
      # For any parameters not set to NULL, warn if value does not match parameter_list
      warning("Supplied value for parameter '", name, "' does not match the value used for function '", function_name, "'.")
    }
  }
  # Set features that are still NULL to defaults
  if (is.null(input)) {
    input <- default_list[[name]]
  }
  return(input)
}

# Require package ---------------------------
#
# Load required package
# Source: ArchR code, Jeffrey Granja & Ryan Corces
#
# x -- Name of package
# load -- Whether to load package
# installInfo -- Installation info
# source -- cran/bioc, etc.
.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(utils::installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if (!is.null(source) & is.null(installInfo)) {
      if (tolower(source) == "cran") {
        installInfo <- paste0('install.packages("',x,'")')
      } else if (tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      } else {
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if (!is.null(installInfo)) {
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    } else {
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

# Calculate centroid distances ---------------------------
#
# Calculate centroid distances for clusters using provided dimensionality
# reduction
#
# reduction -- A dimensionality reduction matrix
# clusters -- Cluster labels corresponding to each row of the dimensionality reduction
.getCentroidDistance <- function(reduction,
                                 clusters) {
  unique_clusters <- unique(clusters)
  n_clusters <- length(unique_clusters)
  centroid_distances <- matrix(0, nrow = n_clusters, ncol = n_clusters)

  # Get centroids
  centroids <- lapply(1:n_clusters,
                      FUN = function(c) {
                        colMeans(reduction[clusters == unique_clusters[c], , drop=FALSE])
                      })

  # Calculate pairwise distances
  for (i in 1:n_clusters) {
    for (j in i:n_clusters) {
      centroid_distances[i,j] <- sqrt(sum((centroids[[i]] - centroids[[j]])^2))
      centroid_distances[j,i] <- centroid_distances[i,j]
    }
  }
  rownames(centroid_distances) <- unique_clusters
  colnames(centroid_distances) <- unique_clusters

  return(centroid_distances)
}

# Generate new cluster labels ---------------------------
#
# Generate new cluster labels that are non-redundant.
#
# merge_groups -- A list indicating which clusters will merge
# level -- Level at which new label will be generated
# compiled_labels -- The set of compiled cluster labels that have been previously used

.getNewLabels <- function(merge_groups,
                          level,
                          compiled_labels) {

  # Create new list
  merge_group_labels <- vector(mode = "list", length(merge_groups))

  # For each element in merge_groups
  for (i in 1:length(merge_groups)) {
    if (is.null(merge_groups[[i]])) {
      merge_group_labels[[i]] <- NULL
    } else if (length(merge_groups[[i]]) == 1) {
      merge_group_labels[[i]] <- merge_groups[[i]]
    } else {
      # If two or more clusters will merge
      # Check whether they share the same tree origin
      roots <- unlist(stringr::str_extract_all(merge_groups[[i]], "P\\d*"))
      if (dplyr::n_distinct(roots) == 1) {
        use_root <- roots[1]
      } else {
        use_root <- "P0"
      }
      # Find the max number under that tree among the parent cluster labels with the matching level
      root_compiled_labels <- compiled_labels[grepl(paste0(use_root, "_L", level, "_"), compiled_labels)]
      if (length(root_compiled_labels) == 0) {
        num <- 1
      } else {
        num <- max(as.numeric(unlist(stringr::str_extract(root_compiled_labels, "\\d*$")))) + 1
      }
      # The subsequent number is used and recorded
      new_label <- paste0(use_root, "_L", level, "_", num)
      merge_group_labels[[i]] <- new_label
      compiled_labels <- c(compiled_labels, new_label)
    }
  }
  # Return list
  return(list("merge_group_labels" = merge_group_labels,
              "compiled_cluster_labels" = compiled_labels))
}

# Startup ---------------------------
#
# Adapted from ArchR code, Jeffrey Granja & Ryan Corces

.onAttach <- function(libname, pkgname){
  # ASCII CHOIR logo
  packageStartupMessage("                           ____
                         /  .-. \\
                         \\  \\ / /
    _____  __    __       /\\ \\/        __   _____
  /      ||  |  |  |    /  /\\ \\       |  | |   _  \\
 |  ,----'|  |__|  |  /  /  .----.    |  | |  |_|  |
 |  |     |   __   | |  |  / .---. \\  |  | |      /
 |  `----.|  |  |  | |  |  `-' \\ \\' | |  | |  |\\  \\__
  \\______||__|  |__|  `. ` .____| |/  |__| | _| `.___|
                         ` -----| |
                         /`.___/ /
                         `------'
")
  # package startup
  v <- utils::packageVersion("CHOIR")
  packageStartupMessage("CHOIR : Version ", v,
                        "\nFor more information see our website : www.CHOIRclustering.com\nIf you encounter a bug please report : https://github.com/CorcesLab/CHOIR/issues")
}
