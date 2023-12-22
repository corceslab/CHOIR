<br>
<a href ="https://www.CHOIRclustering.com"><img src="man/figures/logo.png" width="200px" align="right" /></a>

<!-- badges: start -->
<!-- badges: end -->

**CHOIR** (**c**lustering **h**ierachy **o**ptimization by **i**terative **r**andom forests) is a clustering algorithm for single-cell data. CHOIR applies a framework of permutation tests and random forest classifiers across a hierarchical clustering tree to statistically identify clusters that represent distinct populations.

<br>

## Installation

To install the package, please use the following which requires the [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) package to be installed:
``` r
remotes::install_github("corceslab/CHOIR", ref="main", repos = BiocManager::repositories(), upgrade = "never")
```

Notes:

* Installation should complete in under 2 minutes.
* This package is supported for macOS and Linux. 
* CHOIR depends heavily on the Seurat package, which has been undergoing many changes in recent months. It has been tested successfully with Seurat version 4.3.0 and with SHA c666647.
* Other package dependencies can be found in the "DESCRIPTION" file.

## Usage

Please follow the [vignette](https://www.choirclustering.com/articles/CHOIR.html). The vignette takes less than 10 minutes to run on a standard laptop.

<hr>

<p align="left"><a href ="https://www.corceslab.com/"><img src="man/figures/CorcesLab_logo.png" alt="" width="300"></a></p>

CHOIR is developed and maintained by the [Corces Lab](https://www.corceslab.com/) at the Gladstone Institutes.
