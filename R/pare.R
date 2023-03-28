#' PAREs: PARtial Embeddings for deconfounded dimension reduction
#'
#' `pare` implements a framework for removing nuisance covariates from any
#' distance-based dimension reduction method including \link[Rtsne]{Rtsne} and
#' \link[umap]{umap}. The `pare` function takes either data or a distance matrix
#' as an input then passes a deconfounded distance matrix to a dimension
#' reduction method. The output is a PARE with respect to the specified
#' confounder.
#'
#' @param Y Input data or distance matrix. If `Y` does not have class "dist",
#'   then Euclidean distance is used to construct the distance matrix
#' @param formula Formula for regressing out confounding variables
#' @param data Data frame containing variables in the model specified by `formula`
#' @param embed Dimension reduction method. Must be a function such as
#'   \link[Rtsne]{Rtsne} or \link[umap]{umap}
#' @param type Type of input into `embed`. Refer to documentation of embedding
#'   function
#' @param ... Additional arguments passed to `embed`
#'
#' @return `pare` returns the output of `embed` and the format depends on the
#'   method specified
#' @export
#'
#' @examples
#' pare(iris[,-5], ~ Species, data = iris, cmdscale, type = "dist")
pare <- function(Y, formula, data = NULL, embed, type = c("MDS", "dist"), ...) {
  type <- match.arg(type)

  if (inherits(Y, "dist")) {
    n <- attr(D, "Size")
    D <- Y
  } else {
    n <- dim(Y)[1]
    D <- dist(Y)
  }

  switch(type,
    "MDS" = {
      Z <- cmdscale(D, add = TRUE, k = n-2)$points
      z_formula <- update(formula, Z ~ .)
      data$Z <- Z

      pZ <- lm(z_formula, data)$residuals

      return(do.call(embed, list(pZ, ...)))
    },
    "dist" = {
      D <- as.matrix(D)
      X <- model.matrix(formula, data = data)
      H <- X %*% solve(t(X) %*% X) %*% t(X)
      pD <- as.dist((diag(n) - H) %*% D)

      return(do.call(embed, list(pD, ...)))
    }
  )
}
