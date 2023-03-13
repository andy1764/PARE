pare <- function(Y, formula, data = NULL, embed, type = c("MDS", "dist"),
                 out = FALSE, ...) {
  type <- match.arg(type)
  n <- dim(Y)[1]
  D <- dist(Y)

  switch(type,
    "MDS" = {
      Z <- cmdscale(D, add = TRUE, k = n-2)$points
      z_formula <- update(formula, Z ~ .)

      pZ <- lm(z_formula, data)$residuals

      if (out) { return(pZ) }

      return(
        do.call(embed, list(pZ, ...))
      )
    },
    "dist" = {
      X <- model.matrix(formula, data = data)
      H <- X %*% solve(t(X) %*% X) %*% t(X)
      pD <- as.dist((diag(n) - H) %*% D)

      if (out) { return(pD) }

      return(
        do.call(embed, list(pD, ...))
      )
    }
  )
}
