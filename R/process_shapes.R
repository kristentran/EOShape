#' Process shape data
#'
#' Process shape data before testing. First, symmetrize data using
#' \code{detect_symmetry()} (see below for details); Next, superimpose shapes
#' using procSym. If symmetric pairs are given, pass to procSym and extract
#' symmetric landmarks. If use_PC_scores = TRUE, PC scores will be extracted as
#' shape variables. Lastly, remove shape effects of confounders (special
#' variables included: sex, age, and centroid size).
#'
#' @param borealis borealis object after metadata detection
#' @param symmetrize NULL, "Detect" to use detect.symmetry from MVG's package, or a data frame|matrix with two columns that specify landmark pairings by row. This is provided to the pairedLM argument of Morpho::procSym.
#' @param correction "None", "Pooled", or variable (as column name) to split data.
#' @param silence if TRUE, suppress printed messages from process_shapes and internal calls.
#'
#' @return borealis object for subsequent analyses
#'
#' @export
#'

process_shapes = function(borealis, symmetrize = NULL, correction = "None", silence = FALSE, ...){
  if(!inherits(borealis, "borealis"))
    stop("Object provided is not a borealis object. Please create using initialise_borealis().")
  if(!is.null(symmetrize) && !any(symmetrize == "detect", is.matrix(symmetrize), is.data.frame(symmetrize)))
    stop("symmetrize must be 'detect' or a data frame/matrix with two columns that specify landmark pairings by row.")
  independants = gsub("independant_", "", colnames(borealis$meta)[grepl("independant", colnames(borealis$meta))])
  correction = gsub("independant_", "", correction)
  if(!any(correction %in% c("None", "Pooled", independants)))
    stop("correction argument should be 'None', 'Pooled', or the column name of the group to split data.")
  if(!is.logical(silence) || length(silence) > 1)
    stop("silence parameters must be single value logicals")

  # Symmetry
  # if NULL, don't provide any argument to procSym
  # if "detect" use detect_symmetry
  # if data frame or matrix with two columns, directly send to procSym
  # if 748, get_symmetric_pairs()
  if(!is.null(symmetrize) && symmetrize == "detect"){
    if(dim(borealis$lmks)[1] == 748){
      sympair = sympair_748
      if(!silence)
        cat("748 landmarks detected... using BH lab's symmetric pairings")
    } else {
      # detect
      if(!silence)
        cat("Detecting landmarks... please be sure to evaluate results.")
      sympair = detect.symmetry(borealis$lmks, plot = FALSE)$pairedLM
    }
  } else if(!is.null(symmetrize) && (is.data.frame(symmetrize) || is.matrix(symmetrize()))){
    sympair = symmetrize
  } else {
    sympair = NULL
  }

  # align and extract aligned configurations
  gpa = Morpho::procSym(borealis$lmks, pairedLM = as.matrix(sympair))
  if(is.null(sympair)) aligned = gpa$rotated else aligned = gpa$Sym
  borealis$mshape = geomorph::mshape(aligned)

  borealis$shapes = borealis$lmks = aligned
  borealis$analysis_params$aligned = TRUE

  # if we want to remove size, then add to metadata
  if(borealis$analysis_params$remove_size)
    borealis$meta$size = log(gpa$size)

  # if PCA true, replace lmks with pc_scores
  if(borealis$analysis_params$use_PC_scores){
    borealis$pca = reduce_shapes(borealis$lmks, silence = silence, method = borealis$analysis_params$reduce_shape_dimensions, fromWeb = TRUE)
    borealis$shapes = borealis$pca$x
  } else {
    borealis$pca = reduce_shapes(borealis$lmks, silence = silence, method = 98, fromWeb = TRUE)
  }

  if(correction != "None"){
    if("age" %in% colnames(borealis$meta)){
      borealis$meta$age2 = borealis$meta$age^2
      borealis$meta$age3 = borealis$meta$age^3
    }

    # done this way (as opposed to !independants) to ensure a consistent order of covariates
    covs = cbind(borealis$meta[,colnames(borealis$meta) %in% c("sex", "size"), drop = FALSE],
                 borealis$meta[,grepl("age", colnames(borealis$meta)), drop = FALSE],
                 borealis$meta[,grepl("confounder", colnames(borealis$meta)), drop = FALSE])

    group_split = if(correction == "Pooled") NULL else borealis$meta[,paste0("independant_", correction)]
    mod = residualise(borealis$shapes, covs = covs, interactions = TRUE, group = group_split, silence = silence, ...)

    borealis$shapes = mod$corrected_shapes
    if(correction == "Pooled" || is.null(mod$groups)){
      borealis$analysis_params$correction_model = mod$model
    } else {
      assignments = mod$groups$Assignment
      group_names = as.character(mod$groups$Group)
      models = lapply(assignments, function(a) mod[[a]]$model)
      names(models) = group_names
      borealis$analysis_params$correction_model = models
    }
  }

  borealis$size = if(borealis$analysis_params$log_size) log(gpa$size) else gpa$size

  return(borealis)
}
