#' Isometric log-ratio transformations
#'
#' `ilr_transform()` maps positive compositions to ILR coordinates using the
#' Helmert sub-matrix. `ilr_inverse()` maps ILR coordinates back to positive
#' compositions with the requested closure constant.
#'
#' @param x A vector or matrix of strictly positive compositions. Observations
#'   are rows and parts are columns.
#' @param ilr A vector or matrix of ILR coordinates. Observations are rows.
#' @param kappa A positive closure constant, or a vector of closure constants
#'   with one value per row.
#'
#' @return
#' `ilr_transform()` returns a matrix with one fewer column than the number of
#' compositional parts. `ilr_inverse()` returns a matrix of positive
#' compositions.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 3, 2, 1), ncol = 3, byrow = TRUE)
#' y <- ilr_transform(y)
#' ilr_inverse(y)
#'
#' @export
ilr_transform <- function(x)
{
  x <- if(is.vector(x)) matrix(x, nrow = 1) else as.matrix(x)
  D <- ncol(x)
  # CLR transformation
  clr <- log(x) - rowMeans(log(x))
  # Helmert sub-matrix (D x D-1)
  H <- helmert_basis(D)
  # isometric log-ratios (ILR)
  ilr <- clr %*% H
  dimnames(ilr) <- list(rownames(x), paste0("ilr", seq_len(D - 1L)))
  ilr
}

#' @rdname ilr_transform
#' @export
ilr_inverse <- function(ilr, kappa = 1)
{
  ilr <- if(is.vector(ilr)) matrix(ilr, nrow = 1) else as.matrix(ilr)
  D <- ncol(ilr)+1
  kappa <- as.numeric(kappa)
  # Helmert sub-matrix (D x D-1)
  H <- helmert_basis(D)
  # recover CLR
  clr <- ilr %*% t(H)
  # inverse CLR
  x <- exp(clr)
  x <- x / rowSums(x) * kappa
  dimnames(x) <- NULL
  return(x)
}

helmert_basis <- function(D)
{
  H <- stats::contr.helmert(D)
  scale <- 1 / sqrt(colSums(H^2))
  H %*% diag(scale, nrow = length(scale))
}
