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
  if ("pca" %in% names(x))
    warning("PC space was constructed using original landmarks/sample. Rerun process_shapes if PCA should be redone with new subset.")

  # set i as individuals and j as landmarks
  if (!missing(i)) {

    # if character vector provided, then convert to logical vector
    if(is.character(i)){
      if("identifiers" %in% names(x)){
        i = x$identifiers %in% i
      } else {stop("Please use derive_meta to detect identifiers or provide a numbered or logical index for subsetting.")}
    }

    # 3D
    x$lmks = x$lmks[,,i, drop = FALSE]

    # 2D
    item2d = c("pca", "shapes", "meta")
    for (ite in item2d){
      if(ite %in% names(x)){
        if(ite != "pca")
          x[[ite]] = x[[ite]][i,, drop = FALSE]
        else
        x$pca$x = x$pca$x[i,]}}

    # 1D
    item1d = c("size", "identifiers")
    for (ite in item1d)
      if(ite %in% names(x))
          x[[ite]] = x[[ite]][i, drop = FALSE]

  }

  if (!missing(j)) {
    p = dim(x$lmks)[1]
    k = dim(x$lmks)[2]
    if (is.character(j)) {
      lmk_names = dimnames(x$lmks)[[1]]
      if (is.null(lmk_names))
        stop("Landmark names are not set on lmks; use a numeric or logical index for j.")
      j_lmk = match(j, lmk_names)
      if (anyNA(j_lmk))
        stop("Unknown landmark names in j.")
    } else {
      j_lmk = seq_len(p)[j]
    }
    # two.d.array columns: lmk1x, lmk1y, lmk1z, lmk2x, ...
    j_cols = unlist(lapply(j_lmk, function(jj) ((jj - 1L) * k + 1L):(jj * k)), use.names = FALSE)

    x$lmks = x$lmks[j_lmk, , , drop = FALSE]

    if ("shapes" %in% names(x) && !is.null(x$shapes) &&
        length(dim(x$shapes)) == 2 && !x$analysis_params$use_PC_scores)
      x$shapes = x$shapes[, j_cols, drop = FALSE]
  }
  return(x)
}
