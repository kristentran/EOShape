#' Subset a borealis object
#'
#' `[i, j]` subset borealis object by individuals (i) or landmarks (j)
#'
#' @param x A `borealis` object.
#' @param i,j Subsetting indices.
#' @param ...,drop Additional arguments passed to the method.
#'
#' @return A subset of `x`.
#'
#' @export
#' @method [ borealis

`[.borealis` = function(x, i, j, ..., drop = FALSE) {

  # set i as individuals and j as landmarks
  if (!missing(i)) {

    # if character vector provided, then convert to logical vector
    if(is.character(i)){
      if("identifiers" %in% names(x)){
        i = x$identifiers %in% i
      } else {stop("Please use derive_meta to detect identifiers or instead, provide a numbered or logical index for subsetting.")}
    }

    x$lmks = x$lmks[,,i, drop = FALSE]
    x$meta = x$meta[i,]

    # pca
    if("pca" %in% names(x))
      x$pca$x = x$pca$x[i,]

    x$identifiers = x$identifiers[i]
  }

  if (!missing(j)) {
    x$lmks = x$lmks[j,,, drop = FALSE]
  }
  return(x)
}
