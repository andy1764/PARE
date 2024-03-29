PARE: PARtial Embeddings for deconfounded dimension reduction
================

**Maintainer**: Andrew Chen, <andrewac@pennmedicine.upenn.edu>

**License**: Artistic License 2.0

Partial embeddings (PAREs) provide a framework for removal of
confounding effects from any distance-based dimension reduction method
including $t$-SNE, UMAP, Laplacian Eigenmaps, diffusion map embeddings,
and others. For Euclidean distances, PAREs are equivalent to applying
dimension reduction to the principal coordinates after regressing out
nuisance covariates.

## 1. Installation

The R package can be installed via devtools by running the following
code

```
# install.packages("devtools")
devtools::install_github("andy1764/PARE", build_vignettes = FALSE)
```

Then, you can load this package via

```
library(PARE)
```

## 2. Usage

`PARE` leverages the existing R implementations for dimension reduction
methods. Below is an example call for partial UMAP in the `iris`
dataset, which removes the species effects:

``` r
library(PARE)
library(umap)

pe <- pare(iris[,-5], ~ Species, data = iris, umap)
```

A vignette is provided for the `pare` function, which also covers other
dimension reduction methods. Once the suggested dependencies are installed, please run:

```
devtools::install_github("andy1764/PARE", build_vignettes = TRUE, force = TRUE)
vignette("pare")
```

## 3. Citations

If using the PARE methodology, please cite the following preprint:

> Chen, A. A., Clark, K., Dewey, B., DuVal, A., Pellegrini, N., Nair,
> G., Jalkh, Y., Khalil, S., Zurawski, J., Calabresi, P., Reich, D.,
> Bakshi, R., Shou, H., Shinohara, R. T., Initiative, the A. D. N., &
> Cooperative, the N. A. I. in M. S. (2023). Deconfounded Dimension
> Reduction via Partial Embeddings (p. 2023.01.10.523448). bioRxiv.
> <https://doi.org/10.1101/2023.01.10.523448>
