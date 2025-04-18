% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runCountSplit.R
\name{runCountSplit}
\alias{runCountSplit}
\title{Run countsplitting}
\usage{
runCountSplit(
  object,
  key = "CHOIR",
  use_assay = NULL,
  use_slot = NULL,
  ArchR_matrix = NULL,
  countsplit_suffix = c("_1", "_2"),
  countsplit_params = list(),
  normalization_method = "log",
  verbose = TRUE
)
}
\arguments{
\item{object}{An object of class \code{Seurat}, \code{SingleCellExperiment},
or \code{ArchRProject}. For multi-omic data, we recommend using
\code{ArchRProject} objects.}

\item{key}{The name under which CHOIR-related data for this run is stored in
the object. Defaults to “CHOIR”.}

\item{use_assay}{For \code{Seurat} or \code{SingleCellExperiment} objects, a
character string or vector indicating the assay(s) to use in the provided
object. The default value, \code{NULL}, will choose the current active assay
for \code{Seurat} objects and the \code{logcounts} assay for
\code{SingleCellExperiment} objects.}

\item{use_slot}{For Seurat objects, a character string or vector indicating
the layers(s) — previously known as slot(s) — to use in the provided object.
Default = \code{NULL} will use the 'counts' slot.}

\item{ArchR_matrix}{For \code{ArchR} objects, a character string or vector
indicating which matrix or matrices to use in the provided object. The
default value, \code{NULL}, will use the “GeneScoreMatrix” for ATAC-seq data
or the “GeneExpressionMatrix” for RNA-seq data. For multi-omic datasets,
provide a vector with a value corresponding to each modality.}

\item{countsplit_suffix}{A character vector indicating the suffixes that
distinguish the two count split matrices to be used. Suffixes are appended
onto the input string/vector for parameter \code{use_slot} for \code{Seurat}
objects, \code{use_assay} for \code{SingleCellExperiment} objects, or
\code{ArchR_matrix} for \code{ArchR} objects. Default value \code{NULL} uses
suffixes "_1" and "_2".}

\item{countsplit_params}{A list of additional parameters to be passed to
\code{countsplit::countsplit()}.}

\item{normalization_method}{A character string or vector indicating which
normalization method to apply after countsplitting. Permitted values are
'none', 'log', or 'tfidf'. Defaults to 'log'. For multi-omic datasets,
provide a vector with a value corresponding to each provided value of
\code{use_assay} or \code{ArchR_matrix} in the same order.}

\item{verbose}{A Boolean value indicating whether to use verbose output
during the execution of CHOIR. Defaults to \code{TRUE}, but can be set to
\code{FALSE} for a cleaner output.}
}
\value{
Returns the object including the newly generated countsplit matrices
and the following added data stored under the provided key: \describe{
\item{parameters}{Record of parameter values used}
}
}
\description{
Fetches specified count matri(ces) from the provided input object, generates
countsplit matrices using function \code{countsplit::countsplit()} from A.
Neufeld, and stores these matrices back in the object with suffixes
provided by 'countsplit_suffix'.
}
