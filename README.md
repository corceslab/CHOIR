<br>
<a href ="https://www.CHOIRclustering.com"><img src="man/figures/logo.png" width="200px" align="right" /></a>

<!-- badges: start -->
<!-- badges: end -->

**CHOIR** (**c**lustering **h**ierachy **o**ptimization by **i**terative **r**andom forests) is a clustering algorithm for single-cell data. CHOIR applies a framework of permutation tests and random forest classifiers across a hierarchical clustering tree to statistically identify clusters that represent distinct populations.

<br>

## Installation

CHOIR is designed to be run on Unix-based operating systems such as macOS and linux.

CHOIR installation currently requires `remotes` and `BiocManager` for installation of GitHub and Bioconductor packages. Run the following commands to install the various dependencies used by CHOIR:

First, install remotes (for installing GitHub packages) if it isn’t already installed:
```{r, eval = FALSE}
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
```

Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

Then, install CHOIR:
```{r, eval = FALSE}
remotes::install_github("corceslab/CHOIR", ref="main", repos = BiocManager::repositories(), upgrade = "never")
```

Notes:

* Installation should complete in under 2 minutes.
* This package is supported for macOS and Linux. 
* CHOIR depends heavily on the Seurat package, which has been undergoing many changes in recent months. It has been tested successfully with Seurat version 4.3.0 and with SHA c666647.
* Other package dependencies can be found in the "DESCRIPTION" file.

## Usage

Please follow the [vignette](https://www.choirclustering.com/articles/CHOIR.html). The vignette takes less than 10 minutes to run on a standard laptop.

## How CHOIR works

CHOIR is a hierarchical clustering algorithm that uses permutation testing for cluster
identification by statistical inference. 

<p align="left"><img src="man/figures/Fig1a.png" alt="" width="600"></a></p>

CHOIR identifies clusters that should be merged by applying a permutation test approach to assess the accuracy of random forest classifiers in predicting cluster assignments from a normalized feature matrix.

<p align="left"><img src="man/figures/Fig1b.png" alt="" width="600"></a></p>

CHOIR constructs and iteratively prunes a hierarchical clustering tree using statistical inference to prevent underclustering and overclustering.

<hr>

<p align="left"><a href ="https://www.corceslab.com/"><img src="man/figures/CorcesLab_logo.png" alt="" width="300"></a></p>

CHOIR is developed and maintained by the [Corces Lab](https://www.corceslab.com/) at the Gladstone Institutes.
