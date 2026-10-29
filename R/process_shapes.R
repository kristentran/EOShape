#' Process shape data
#'
#' Process shape data before testing. First, symmetrize data using
#' \code{detect_symmetry()} (see below for details); Next, superimpose shapes
#' using procSym. If symmetric pairs are given, pass to procSym and extract
#' symmetric landmarks. If use_PC_scores = TRUE, PC scores will be extracted as
#' shape variables. Lastly, remove shape effects of confounders (special
#' variables included: sex, age, and centroid size).
#'
#' @param eos eos object after metadata detection
#' @param symmetrize NULL, "Detect" to use detect.symmetry from MVG's package, or a data frame|matrix with two columns that specify landmark pairings by row. This is provided to the pairedLM argument of Morpho::procSym.
#' @param correction "None", "Pooled", or variable (as column name) to split data.
#' @param silence if TRUE, suppress printed messages from process_shapes and internal calls.
#'
#' @return eos object for subsequent analyses
#'
#' @export
#'

process_shapes = function(eos, symmetrize = NULL, correction = "None", silence = FALSE, ...){
  if(!inherits(eos, "eos"))
    stop("Object provided is not an eos object. Please create using initialise_eos().")
  if(!is.null(symmetrize) && !any(symmetrize == "detect", is.matrix(symmetrize), is.data.frame(symmetrize)))
    stop("symmetrize must be 'detect' or a data frame/matrix with two columns that specify landmark pairings by row.")
  independants = gsub("independant_", "", colnames(eos$meta)[grepl("independant", colnames(eos$meta))])
  correction = gsub("independant_", "", correction)
  if(!any(correction %in% c("None", "Pooled", independants)))
    stop("correction argument should be 'None', 'Pooled', or the column name of the group to split data.")
  if(!is.logical(silence) || length(silence) > 1)
    stop("silence parameters must be single value logicals")

  log = c()

  if(!eos$analysis_params$aligned){
  ## SYMMETRIZE ##
  # if NULL, don't provide any argument to procSym
  # if "detect" use detect_symmetry
  # if data frame or matrix with two columns, directly send to procSym
  # if 748, get_symmetric_pairs()
    if(!is.null(symmetrize) && symmetrize == "detect"){
      if(dim(eos$lmks)[1] == 748){
        sympair = as.matrix(sympair_748)
        if(!silence)
          cat("748 landmarks detected... using BH lab's symmetric pairings")
        log = c(log, "748 landmarks detected... used BH lab's symmetric pairings to symmetrize data.")
      } else {
        # detect
        if(!silence)
          cat("Detecting landmarks... please be sure to evaluate results.")
        sympair = as.matrix(detect.symmetry(eos$lmks, plot = FALSE)$pairedLM)
        log = c(log, "Symmetric landmark pairs automatically detected to symmetrize data.")
      }
    } else if(!is.null(symmetrize) && (is.data.frame(symmetrize) || is.matrix(symmetrize()))){
      sympair = as.matrix(symmetrize)
      log = c(log, "Provided symmetric pairs used to symmetrize data.")
    } else {
      sympair = NULL
      log = c(log, "Data was NOT symmetrized.")
    }

    ## PROCRUSTES SUPERIMPOSITION ##
    prealigned = eos$lmks
    # if atlas_lmks exist and are the same landmark set as data, then align to atlas during superimposition
    if(!is.null(eos$atlas_lmks) && all(dim(eos$atlas_lmks) == dim(eos$lmks)[1:2])){
      prealigned = abind::abind(eos$atlas_lmks, prealigned, along = 3)
      log = c(log, "Shapes aligned to atlas.")
      pcAlign = FALSE # Remember that the Morpho documentation seems to be backwards for this setting
    } else {pcAlign = TRUE}

    # align and extract aligned configurations
    gpa = Morpho::procSym(prealigned, pairedLM = sympair, pcAlign = pcAlign)
    if(is.null(sympair)) aligned = gpa$rotated else aligned = gpa$Sym

    if(!is.null(eos$atlas_lmks) && all(dim(eos$atlas_lmks) == dim(eos$lmks)[1:2]))
      aligned = aligned[,,-1]

    eos$analysis_params$aligned = TRUE
  } else {
    log = c(log, "Data marked as aligned... skipping symmetrization and Procrustes superimposition.")
    aligned = eos$lmks
  }

  # store Grand mean shape and aligned landmarks
  eos$mshape = geomorph::mshape(aligned)
  eos$lmks = aligned
  eos$shapes = geomorph::two.d.array(aligned) # $shapes is expected in 2D format; pc scores or coordinates

  # save size according to parameters, untransformed vs log-transformed
  if(exists("gpa")){
    eos$size = if(eos$analysis_params$log_size) log(gpa$size) else gpa$size

    # if size is a confounder, then add to metadata
    if(eos$analysis_params$remove_size)
      eos$meta$size = eos$size
  }

  ## PCA ##
  # if use PC scores if true, then replace shapes with pc_scores
  if(eos$analysis_params$use_PC_scores){
    # dim reduction is changable from analysis params
    eos$pca = reduce_shapes(eos$lmks, silence = silence, method = eos$analysis_params$reduce_shapes_method, fromWeb = TRUE)
    eos$shapes = eos$pca$x
    log = c(log, paste0("Shape space reduced to ", ncol(eos$pca$x), "PCs."))
  } else {
    # automatic PC selection is to truncate at 98% cumulative variance - silence this
    eos$pca = reduce_shapes(eos$lmks, silence = TRUE, method = 98, fromWeb = TRUE)
  }

  # get landmarks along first 4 PCs (capped at as many PCs exist)
  vpca = visualise_pcs(eos$pca, pcs = 1:min(4, ncol(eos$pca$x)), mean_shape = eos$mshape, fromWeb = TRUE)
  eos$pca$morphs = vpca$morphs

  # Set up correction
  if(correction != "None"){
    # Age needs polynomial terms
    if("age" %in% colnames(eos$meta)){
      eos$meta$age2 = eos$meta$age^2
      eos$meta$age3 = eos$meta$age^3
    }

    # Create covariate data frame
    # done this way (as opposed to !independants) to ensure a consistent order of covariates
    covs = cbind(eos$meta[,colnames(eos$meta) %in% c("sex", "size"), drop = FALSE],
                 eos$meta[,grepl("age", colnames(eos$meta)), drop = FALSE],
                 eos$meta[,grepl("confounder", colnames(eos$meta)), drop = FALSE])

    # set up groups if correction is not set to pooled
    group_split = if(correction == "Pooled") NULL else eos$meta[,paste0("independant_", correction)]
    mod = residualise(eos$shapes, covs = covs, interactions = TRUE, group = group_split, silence = silence, ...)
    # save corrected shapes and log from residualise
    eos$shapes = mod$corrected_shapes
    log = c(log, mod$log)

    # Save model to eos object
    if(correction == "Pooled" || is.null(mod$groups)){
      eos$correction_model = mod$model
    } else {
      assignments = mod$groups$Assignment
      group_names = as.character(mod$groups$Group)
      models = lapply(assignments, function(a) mod[[a]]$model)
      names(models) = group_names
      eos$correction_model = models
    }
  }

  eos$analysis_log = c(eos$analysis_log, log)
  return(eos)
}
