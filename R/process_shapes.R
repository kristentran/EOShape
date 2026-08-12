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
#' @param symmetrize borealis object after metadata detection
#' @param automated if TRUE, ...
#' @param silence if TRUE, suppress printed messages from process_shapes and internal calls.
#'
#' @return borealis object for subsequent analyses
#'
#' @export
#'

process_shapes = function(borealis, symmetrize = NULL, correction = "None", automated = FALSE, silence = FALSE, ...){
  if(!inherits(borealis, "borealis"))
    stop("Object provided is not a borealis object. Please create using initialise_borealis().")
  if(!is.null(symmetrize) && !any(symmetrize == "detect", is.matrix(symmetrize), is.data.frame(symmetrize)))
    stop("symmetrize must be 'detect' or a data frame/matrix with two columns that specify landmark pairings by row.")
  dependants = gsub("dependant_", "", colnames(borealis$meta)[grepl("dependant", colnames(borealis$meta))])
  correction = gsub("dependant_", "", correction)
  if(!any(correction %in% c("None", "Pooled", dependants)))
    stop("correction argument should be 'None', 'Pooled', or the column name of the group to split data.")
  if(!is.logical(automated) || !is.logical(silence) || length(automated) > 1 || length(silence) > 1)
    stop("automated and silence parameters must be single value logicals")

  # symmetry
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
    borealis$meta$logsize = log(gpa$size)

  # if PCA true, replace lmks with pc_scores
  if(borealis$analysis_params$use_PC_scores){
    borealis$pca = reduce_shapes(borealis$lmks, silence = silence, ...)
    borealis$shapes = borealis$pca$x
  }

  if(correction != "None"){
    if("age" %in% colnames(borealis$meta)){
      borealis$meta$age2 = borealis$meta$age^2
      borealis$meta$age3 = borealis$meta$age^3
    }

    # done this way (as opposed to !dependants) to ensure a consistent order of covariates
    covs = cbind(borealis$meta[,colnames(borealis$meta) %in% c("sex", "logsize"), drop = FALSE],
                 borealis$meta[,grepl("age", colnames(borealis$meta)), drop = FALSE],
                 borealis$meta[,grepl("confounder", colnames(borealis$meta)), drop = FALSE])

    group_split = if(correction == "Pooled") NULL else borealis$meta[,paste0("dependant_", correction)]
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

  borealis$logsize = log(gpa$size)

  return(borealis)
}
