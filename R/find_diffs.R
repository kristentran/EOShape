#' Find differences
#'
#' Are there differences in csize and shape by independent variable(s)? Loops
#' \code{scrutinise()} through shape and size; see \code{scrutinise()} for
#' details. If only one independent variable is provided and it is categorical,
#' then post-hoc testing will be completed (if that variable has more than 2
#' levels). * note that this is an analysis fx, so the object returned are just
#' results, not the eos object
#'
#' TODO: get landmark coordinates for HOW they differ if significant. If
#' continuous, morph along axis from procd.lm obj. If categorical, CVA?
#' TODO: If >1 independent var, use evaluate model to get assess sig. variables.
#' Then only model the significant terms
#'
#' If differences exist, then how are they different? If continuous, then
#'
#' @param eos eos object after metadata detection
#' @param independent default to 'All' specified independent variables. Specify variable if only one is desired.
#' @param interactions if TRUE and there are multiple independent variables, ...
#' @param silence if TRUE, suppress printed messages from process_shapes and internal calls.
#'
#' @return eos object for subsequent analyses
#'
#' @export

find_diffs = function(eos, independent = "All", interactions = TRUE, silence = FALSE, ...){
  if(!inherits(eos, "eos"))
    stop("Object provided is not an eos object. Please create using initialise_eos().")
  if(!is.logical(silence) || length(silence) > 1)
    stop("silence parameters must be single value logicals")
  if(independent != "All"){
    independents = gsub("independent_", "", colnames(eos$meta)[grepl("independent", colnames(eos$meta))])
    if(!any(gsub("independent_", "", independent) %in% independents))
      stop("independent argument should be 'All' or the name of variable to analyse (e.g. Experimental_Group, or Treatment_Dosage).")
    rm(independents)
  }

  # get column names from meta data for independent variables
  if(independent == "All"){
    independents = colnames(eos$meta)[grepl("independent", colnames(eos$meta))]
  } else {
    if(!grepl("independent_", independent))
      independent = paste0("independent_", independent)
    independents = independent
  }

  # get data frame of only independent variables
  independent_df = eos$meta[,independents, drop = FALSE]

  # analyse shape and size differences
  dependents = c("shapes", "size")

  # PC scores need the mean shape and rotation to visualise differences as landmarks
  if(eos$analysis_params$use_PC_scores){
    mean_shape = eos$mshape
    rotation = eos$pca$rotation
  } else {mean_shape = rotation = NULL}

  # object to store results
  analysis_results = list()
  vis_diffs = NULL

  # loop through size and shape as Y
  for(dependent in dependents){
    covs = NULL
    if(dependent != "size" && eos$analysis_params$assess_allometry)
      covs = data.frame(size = eos$size)

    res = scrutinise(eos[[dependent]], independents = independent_df, covs = covs, interactions = interactions,
                     mean_shape = mean_shape, pca_rotation = rotation, silence = silence)
    analysis_results[[dependent]] = res[c("full_aov", "post_hoc")]
    if(dependent == "shapes")
      vis_diffs = res$landmark_diffs
  }
  return(list(anovas = analysis_results, landmark_diffs = vis_diffs))
}
