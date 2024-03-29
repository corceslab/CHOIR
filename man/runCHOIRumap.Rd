% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingUtils.R
\name{runCHOIRumap}
\alias{runCHOIRumap}
\title{Run UMAP on CHOIR dimensionality reductions}
\usage{
runCHOIRumap(object, key = "CHOIR", reduction = NULL, verbose = TRUE)
}
\arguments{
\item{object}{An object of class 'Seurat', 'SingleCellExperiment', or
'ArchRProject' that has undergone CHOIR clustering.}

\item{key}{The name under which CHOIR-related data for this run is retrieved
from the object. Defaults to 'CHOIR'.}

\item{reduction}{A character vector indicating which CHOIR subtree
dimensionality reduction to run UMAP on (e.g., 'P0_reduction',
'P1_reduction'). Default = \code{NULL} will run UMAP on all of the
dimensionality reductions generated by CHOIR stored under the provided 'key'.}

\item{verbose}{A boolean value indicating whether to use verbose output
during the execution of this function. Can be set to \code{FALSE} for a
cleaner output.}
}
\value{
Returns the object with the following added data stored under the
provided key: \describe{
\item{reduction}{Cell embeddings for all calculated dimensionality
reductions}
}
}
\description{
Simplifies running \code{Seurat::RunUMAP()} after CHOIR clustering by
automatically fetching the pre-generated dimensionality reductions.
}
