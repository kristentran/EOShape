#' Measure linear distances and angles
#'
#' Wrapper around \code{MusPO::MeasureAll()}. Landmarks are \code{eos$lmks} at
#' the current stage, and \code{group} is an independent variable set by
#' \code{derive_meta()}. \code{process_shapes()} is not required. The function
#' can be run at any stage once metadata have been derived.
#'
#' @param eos eos object. \code{derive_meta()} must already have been run.
#' @param independent independent variable to pass as \code{group}. \code{NULL}
#'   uses the only independent variable. Supply a name when
#'   \code{derive_meta()} set more than one. The \code{independent_} prefix is
#'   optional.
#' @param silence if TRUE, suppress messages from \code{MusPO::MeasureAll()}.
#' @param ... further arguments passed to \code{MusPO::MeasureAll()}, such as
#'   \code{control}, \code{traits}, \code{region}, \code{covariates},
#'   \code{alpha}, \code{n_perm}, and \code{seed}. \code{shapes} and
#'   \code{group} are set from the eos object.
#'
#' @return The list returned by \code{MusPO::MeasureAll()}.
#'
#' @export

caliper = function(eos, independent = NULL, silence = FALSE, ...){
  if(!inherits(eos, "eos"))
    stop("Object provided is not an eos object. Please create using initialise_eos().")
  if(!is.logical(silence) || length(silence) != 1 || is.na(silence))
    stop("silence must be a single logical value.")
  if(any(names(list(...)) %in% c("shapes", "group")))
    stop("shapes and group are taken from the eos object. Pass other MeasureAll arguments through ...")
  if(!requireNamespace("MusPO", quietly = TRUE))
    stop("Package 'MusPO' is required for caliper(). Install it from ~/Documents/MusPO.")

  if(!"meta" %in% names(eos) || !is.data.frame(eos$meta))
    stop("Please run derive_meta() before caliper(). process_shapes() is not required.")
  if(!"identifiers" %in% names(eos))
    stop("Please run derive_meta() before caliper(). process_shapes() is not required.")
  if(is.null(eos$lmks))
    stop("eos$lmks is missing.")

  ind_cols = colnames(eos$meta)[grepl("^independent_", colnames(eos$meta))]
  ind_names = sub("^independent_", "", ind_cols)
  if(!length(ind_cols))
    stop("No independent variable in eos$meta. Set one in derive_meta() before caliper().")

  if(is.null(independent)){
    if(length(ind_cols) != 1)
      stop("More than one independent variable is set (", paste(ind_names, collapse = ", "), "). Pass independent = one of these names.")
    independent_col = ind_cols
  } else {
    if(!is.character(independent) || length(independent) != 1 || is.na(independent) || !nzchar(independent))
      stop("independent must be a single variable name.")
    independent = sub("^independent_", "", independent)
    if(!independent %in% ind_names)
      stop("independent must be one of: ", paste(ind_names, collapse = ", "), ".")
    independent_col = paste0("independent_", independent)
  }

  shapes = eos$lmks
  if(length(dim(shapes)) == 3){
    n = dim(shapes)[3]
    specimens = dimnames(shapes)[[3]]
  } else if(length(dim(shapes)) == 2 && ncol(shapes) == 3){
    n = 1L
    specimens = NULL
  } else {
    stop("eos$lmks must be a p x 3 matrix or a p x 3 x n array.")
  }

  if(!is.null(specimens)){
    if(is.null(rownames(eos$meta)) || anyNA(match(specimens, rownames(eos$meta))))
      stop("Specimen names on eos$lmks are not all present in eos$meta.")
    group = eos$meta[specimens, independent_col]
  } else if(nrow(eos$meta) != n){
    stop("eos$lmks and eos$meta have different numbers of specimens.")
  } else {
    group = eos$meta[[independent_col]]
  }

  measure = function() MusPO::MeasureAll(shapes, group = group, ...)
  if(silence) suppressMessages(measure()) else measure()
}
