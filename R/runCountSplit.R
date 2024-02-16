#' Run countsplitting
#'
#' Fetches specified count matri(ces) from the provided input object, generates
#' countsplit matrices using function \code{countsplit::countsplit()} from A.
#' Neufeld, and stores these matrices back in the object with suffixes
#' provided by 'countsplit_suffix'.
#'
#' @param object An object of class 'Seurat', 'SingleCellExperiment', or
#' 'ArchRProject'.
#' @param key The name under which CHOIR-related data for this run is stored in
#' the object. Defaults to 'CHOIR'.
#' @param use_assay For Seurat or SingleCellExperiment objects, a character
#' string or vector indicating the assay(s) to use in the provided object.
#' Default = \code{NULL} will choose the current active assay for Seurat objects
#' and the \code{counts} assay for SingleCellExperiment objects.
#' @param use_slot For Seurat objects, a character string or vector indicating
#' the layers(s) — previously known as slot(s) — to use in the provided object.
#' Default = \code{NULL} will use the 'counts' slot.
#' @param ArchR_matrix For ArchR objects, a character string or vector
#' indicating which matri(ces) to use in the provided object. Default =
#' \code{NULL} will use the 'TileMatrix' for ATAC-seq data or the
#' 'GeneExpressionMatrix' for RNA-seq data.
#' @param countsplit_suffix A character vector indicating the suffixes
#' that distinguish the two countsplit matrices to be used. Suffixes are
#' appended onto input string/vector for \code{use_slot} for Seurat objects,
#' \code{use_assay} for SingleCellExperiment objects, or \code{ArchR_matrix} for
#' ArchR objects. Defaults to suffixes "_1" and "_2".
#' @param countsplit_params A list of additional parameters to be passed to
#' \code{countsplit::countsplit()}.
#' @param normalization_method A character string or vector indicating which
#' normalization method to apply after countsplitting. Permitted values are
#' 'none', 'log', or 'tfidf'. Defaults to 'log'.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Can be set to \code{FALSE} for a
#' cleaner output.
#'
#' @return Returns the object including the newly generated countsplit matrices
#' and the following added data stored under the provided key: \describe{
#'   \item{parameters}{Record of parameter values used}
#'   }
#'
#' @export
#'
runCountSplit <- function(object,
                          key = "CHOIR",
                          use_assay = NULL,
                          use_slot = NULL,
                          ArchR_matrix = NULL,
                          countsplit_suffix = c("_1", "_2"),
                          countsplit_params = list(),
                          normalization_method = "log",
                          verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  # Progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                       " : Checking inputs and preparing object..")

  # Require package countsplit
  .requirePackage("countsplit", source = "cran")

  .validInput(object, "object", "runCountSplit")
  .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
  .validInput(use_slot, "use_slot", list(object, use_assay, FALSE, NULL))
  .validInput(ArchR_matrix, "ArchR_matrix", list(object, FALSE, NULL))
  .validInput(countsplit_suffix, "countsplit_suffix", TRUE)
  .validInput(countsplit_params, "countsplit_params")
  .validInput(verbose, "verbose")

  # Number of modalities & object type
  if (methods::is(object, "ArchRProject")) {
    # n_modalities <- max(length(ArchR_matrix), 1)
    # object_type <- "ArchRProject"
    # .requirePackage("ArchR", installInfo = "Instructions at archrproject.com")
    stop("Countsplitting using function 'runCountSplit' is not yet supported for ArchR objects.")
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
      if (is.null(use_slot)) {
        use_slot <- "counts"
        .validInput(use_slot, "use_slot", list(object, use_assay, FALSE, NULL))
      }
    } else if (methods::is(object, "SingleCellExperiment")) {
      object_type <- "SingleCellExperiment"
      .requirePackage("SingleCellExperiment", source = "bioc")
    }
  }

  # ---------------------------------------------------------------------------
  # Countsplitting
  # ---------------------------------------------------------------------------

  for (m in 1:n_modalities) {
    # Fetch matrix
    # Progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                         " : Fetching matrix ", m, " of ", n_modalities, "..")
    # Match input arguments
    use_assay_m <- .matchArg(use_assay, m)
    use_slot_m <- .matchArg(use_slot, m)
    ArchR_matrix_m <- .matchArg(ArchR_matrix, m)
    normalization_method_m <- .matchArg(normalization_method, m)
    matrix_m <- .getMatrix(object = object,
                           use_assay = use_assay_m,
                           use_slot = use_slot_m,
                           ArchR_matrix = ArchR_matrix_m,
                           verbose = TRUE)
    # Warn if not integer values
    # Check up to 100 random cells
    sampled_cells <- sample(colnames(matrix_m), min(100, ncol(matrix_m)))
    sample_matrix <- as(matrix_m[,sampled_cells], "dgCMatrix")
    if (!all(sample_matrix == floor(sample_matrix))) {
      warning("Countsplitting is not intended for non-integer matrices, counts will be rounded.")
    }
    # Generate countsplit matrices
    # Progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                         " : Countsplitting matrix..")
    matrix_m_split <- do.call(countsplit::countsplit,
                              c(list("X" = matrix_m),
                                countsplit_params))
    # Order by original row & column names
    matrix_m_1 <- matrix_m_split[[1]][rownames(matrix_m), colnames(matrix_m)]
    matrix_m_2 <- matrix_m_split[[2]][rownames(matrix_m), colnames(matrix_m)]

    # Store matrices
    # Progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                         " : Storing countsplit matrices..")
    object <- .storeMatrix(object = object,
                           use_matrix = matrix_m_1,
                           use_assay = `if`(object_type == "SingleCellExperiment", paste0(use_assay_m, countsplit_suffix[1]), use_assay_m),
                           use_slot = `if`(object_type == "Seurat", paste0(use_slot_m, countsplit_suffix[1]), use_slot_m),
                           ArchR_matrix = `if`(object_type == "ArchRProject", paste0(ArchR_matrix_m, countsplit_suffix[1]), ArchR_matrix_m),
                           verbose = TRUE)
    object <- .storeMatrix(object = object,
                           use_matrix = matrix_m_2,
                           use_assay = `if`(object_type == "SingleCellExperiment", paste0(use_assay_m, countsplit_suffix[2]), use_assay_m),
                           use_slot = `if`(object_type == "Seurat", paste0(use_slot_m, countsplit_suffix[2]), use_slot_m),
                           ArchR_matrix = `if`(object_type == "ArchRProject", paste0(ArchR_matrix_m, countsplit_suffix[2]), ArchR_matrix_m),
                           verbose = TRUE)
    # Normalize matrices
    if (normalization_method_m != "none") {
      # Progress
      if (verbose == TRUE) message(format(Sys.time(), "%Y-%m-%d %X"),
                                   " : Running ", normalization_method_m,
                                   " normalization on countsplit matrices..")
      if (normalization_method_m == "log") {
        matrix_m_1_norm <- Seurat::NormalizeData(matrix_m_1, verbose = FALSE)
        matrix_m_2_norm <- Seurat::NormalizeData(matrix_m_2, verbose = FALSE)
      } else if (normalization_method_m == "tfidf") {
        matrix_m_1_norm <- Signac::RunTFIDF(matrix_m_1, verbose = FALSE)
        matrix_m_2_norm <- Signac::RunTFIDF(matrix_m_2, verbose = FALSE)
      }
      # Store normalized matrices
      # Progress
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                           " : Storing normalized countsplit matrices..")
      object <- .storeMatrix(object = object,
                             use_matrix = matrix_m_1_norm,
                             use_assay = `if`(object_type == "SingleCellExperiment",
                                                paste0(use_assay_m, "_", normalization_method_m, countsplit_suffix[1]),
                                                use_assay_m),
                             use_slot = `if`(object_type == "Seurat",
                                               paste0(use_slot_m, "_", normalization_method_m, countsplit_suffix[1]),
                                               use_slot_m),
                             ArchR_matrix = `if`(object_type == "ArchRProject",
                                                   paste0(ArchR_matrix_m, "_", normalization_method_m, countsplit_suffix[1]),
                                                   ArchR_matrix_m),
                             verbose = TRUE)
      object <- .storeMatrix(object = object,
                             use_matrix = matrix_m_2_norm,
                             use_assay = `if`(object_type == "SingleCellExperiment",
                                                paste0(use_assay_m, "_", normalization_method_m, countsplit_suffix[2]),
                                                use_assay_m),
                             use_slot = `if`(object_type == "Seurat",
                                               paste0(use_slot_m, "_", normalization_method_m, countsplit_suffix[2]),
                                               use_slot_m),
                             ArchR_matrix = `if`(object_type == "ArchRProject",
                                                   paste0(ArchR_matrix_m, "_", normalization_method_m, countsplit_suffix[2]),
                                                   ArchR_matrix_m),
                             verbose = TRUE)
    }
  }
  # Record parameters used and add to original object
  parameter_list <- list("use_assay" = use_assay,
                         "use_slot" = use_slot,
                         "ArchR_matrix" = ArchR_matrix,
                         "countsplit_suffix" = countsplit_suffix,
                         "countsplit_params" = countsplit_params,
                         "normalization_method" = normalization_method)
  object <- .storeData(object, key, "parameters", parameter_list, "runCountSplit_parameters")

  # Standard object size warning
  warning("Adding countsplit matrices may dramatically increase object size. Consider removing matrices that are not needed, as CHOIR will run faster on smaller objects.")

  # Return object with added matrices
  return(object)
}
