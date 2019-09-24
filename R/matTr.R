matTr <- function(a,b) {
  # trace of a %*% b; dimensions of a and b should match.
  if (!identical(dim(a), dim(b))){
    stop("Dimensions do not match!")
  } else {
    crossprodCpp(as.vector(a), as.vector(t(b)))
  }
}

matTr2 <- function(z) sum(diag(z))
