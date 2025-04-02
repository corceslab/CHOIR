#' Generate root dimensionality reduction
#'
#' This function generates the preliminary dimensionality reduction that can
#' subsequently be used by CHOIR to construct a hierarchical clustering tree.
#'
#' For multi-modal data, optionally supply parameter inputs as vectors/lists
#' that sequentially specify the value for each modality.
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{ArchRProject}. For multi-omic data, we recommend using
#' \code{ArchRProject} objects.
#' @param key The name under which CHOIR-related data for this run is stored in
#' the object. Defaults to “CHOIR”.
#' @param normalization_method A character string or vector indicating which
#' normalization method to use. In general, input data should be supplied to
#' CHOIR after normalization, except when the user wishes to use
#' \code{Seurat SCTransform} normalization. Permitted values are “none” or
#' “SCTransform”. Defaults to “none”. Because CHOIR has not been tested
#' thoroughly with \code{SCTransform} normalization, we do not recommend this
#' approach at this time. For multi-omic datasets, provide a vector with a value
#' corresponding to each provided value of \code{use_assay} or
#' \code{ArchR_matrix} in the same order.
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
#' minus 2. This function uses parallelization only for \code{ArchR} objects.
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
#'   \item{parameters}{Record of parameter values used}
#'   \item{reduction}{Cell embeddings for the calculated dimensionality
#'   reduction}
#'   \item{var_features}{Variable features the calculated dimensionality
#'   reduction}
#'   }
#'
#' @export
#'
runRootReduction <- function(object,
                      key = "CHOIR",
                      normalization_method = "none",
                      reduction_method = NULL,
                      reduction_params = list(),
                      n_var_features = NULL,
                      batch_correction_method = "none",
                      batch_correction_params = list(),
                      batch_labels = NULL,
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

  # Progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                       " : (Step 1/2) Checking inputs and preparing object..")

  .validInput(object, "object", "buildTree")
  .validInput(key, "key", list("buildTree", object))
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

  .validInput(ArchR_depthcol, "ArchR_depthcol", list(object, n_modalities))
  .validInput(reduction, "reduction", list("buildTree", object))
  .validInput(var_features, "var_features", reduction)
  .validInput(atac, "atac", n_modalities)
  .validInput(normalization_method, "normalization_method", list(object, n_modalities, use_assay))
  .validInput(reduction_method, "reduction_method", list(object, n_modalities))
  .validInput(reduction_params, "reduction_params", list(object, ArchR_matrix, use_assay, reduction_method))
  .validInput(n_var_features, "n_var_features", n_modalities)
  .validInput(batch_correction_method, "batch_correction_method", n_modalities)
  .validInput(batch_correction_params, "batch_correction_params", list(object, ArchR_matrix, use_assay, batch_correction_method))
  .validInput(batch_labels, "batch_labels", object)
  .validInput(n_cores, "n_cores")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
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
    }
  } else {
    countsplit_suffix <- c("", "")
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
                       "\n - Normalization method: ", normalization_method,
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
                       "\n - # of cores: ", n_cores,
                       "\n - Random seed: ", random_seed,
                       "\n")

  # ---------------------------------------------------------------------------
  # Step 2: Initial dimensionality reduction
  # ---------------------------------------------------------------------------

  # Run dimensionality reduction if not supplied by user
  if (is.null(reduction)) {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/2) Running initial dimensionality reduction..")
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
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : (Step 2/2) Setting initial dimensionality reduction..")
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
  # Store data in object
  # ---------------------------------------------------------------------------

  # Record parameters used and add to original object
  parameter_list <- list("normalization_method" = normalization_method,
                         "reduction_method" = reduction_method,
                         "reduction_params" = reduction_params,
                         "n_var_features" = n_var_features,
                         "batch_correction_method" = batch_correction_method,
                         "batch_correction_params" = batch_correction_params,
                         "batch_labels" = batch_labels,
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

  object <- .storeData(object, key, "parameters", parameter_list, "runRootReduction_parameters")

  # Return object with new additions
  return(object)
}
