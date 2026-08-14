#' Initialize data object for automated analysis
#'
#' Initialize borealis object with a list object containing "lmks" as 3D array and "meta" as a data.frame.
#' Constant variables in meta will be dropped (i.e., variables with the same value across all individuals).
#'
#' @param data list object containing "lmks" as 3D array and "meta" as a data frame.
#' @param atlas_mesh mesh of atlas for visualization. Provide if using package from R.
#' @param atlas_lmks landmarks associated with atlas_mesh. Provide if using package from R and if landmarks =/= vertices from mesh.
#'
#' @return borealis object for subsequent analyses using borealis package
#'
#' @import dplyr
#'
#' @export
#'

initialise_borealis = function(data, atlas_mesh = NULL, atlas_lmks = NULL, silence = FALSE, analysis_params = NULL){
  if(!is.list(data) || !all(c("lmks", "meta") %in% names(data)))
    stop("data must be a list object that contains meta (data frame) and lmks (array)")
  if(!is.null(atlas_lmks) && !is.matrix(atlas_lmks))
    stop("Atlas landmarks must be in a matrix")
  if(!is.null(atlas_mesh) && "mesh3d" %in% class(atlas_mesh))
    stop("Atlas mesh must be in a matrix")

  require(dplyr)
  require(crayon)

  class(data) = "borealis"

  # store analysis parameters
  if(is.null(analysis_params)){
    analysis_params = list(
      aligned = FALSE,
      use_PC_scores = FALSE,
      imputed_meta = FALSE,
      remove_size = TRUE,
      assess_allometry = FALSE,
      log_size = FALSE
    )
  } else {
    analysis_params$aligned = FALSE
    analysis_params$imputed_meta = FALSE
  }

  if(analysis_params$remove_size && analysis_params$assess_allometry)
    stop("Parameters remove_sizet' and assess_allometry cannot both be TRUE. Please reset analysis parameters.")

  # this will store analysis stuff to be passed from API
  analysis_log = character(0)


  # drop constant and NA-only meta variables
  constant_cols = names(select(data$meta, where(~ n_distinct(., na.rm = TRUE) == 1 | all(is.na(.)))))
  if (length(constant_cols) && !silence)
    cat(crayon::blue("Removing constant and NA variables:"), paste(constant_cols, collapse = ", "), "\n")

  data$meta = select(data$meta, -any_of(constant_cols))

  data$lmks = convert(data$lmks, silence = silence)

  # store things inside data
  if(!is.null(atlas_mesh)){
    data[["atlas_mesh"]] = atlas_mesh
    if(is.null(atlas_lmks)){
      if(ncol(atlas_mesh$vb) == dim(data$lmks)[1]){
        analysis_log = c(analysis_log, "Atlas landmarks taken from mesh's vertices")
        data[["atlas_lmks"]] = t(atlas_mesh$vb[1:3,])
      } else {stop("Atlas_mesh provided without corresponding landmarks and has different number of vertices than landmark data. Please provide corresponding landmarks for atlas.")}
    }
  }
  if(!is.null(atlas_lmks)){
    data[["atlas_lmks"]] = atlas_lmks}

  data[["analysis_params"]] = analysis_params
  data[["analysis_log"]] = analysis_log

  return(data)
}
