#' Evaluate linear models
#'
#' Correct shapes by chosen covariates using a linear model. Uses workhorse
#' function evaluate_model() to determine which covariates are significant and
#' will be included in the model. If there are ANY significant terms in the
#' final evaluated model, then the mean shape will be added back to the
#' residuals of the model.
#' TODO: fix the force to run through evaluate_model to allow for more control.
#'
#' @param shapes 2d matrix or 3d array
#' @param covs data frame with identifiers as rownames and all covariates to be included in the model as columns.
#' @param interactions if TRUE, include all possible interactions.
#' @param p_threshold numerical between 0 and 1. Determine required significance of covariates before removing from the model
#' @param model procD.lm object. if provided, correct based provided model (bypass evaluate_model())
#' @param force if TRUE, then the correction will always occur after evaluate_model().
#' @param group vector of same length as shapes and cov. If provided, performs shape correction independently (within groups).
#' @param silence if TRUE, suppress printed messages from evaluate_model().
#'
#' @return model model used for shape correction.
#' @return corrected_shapes residuals of model added back to the mean shape
#'
#' @import geomorph
#'
#' @export

# Different covariates - either loop function, or combine covs

# Thanks to Eva for the fun name
residualise = function(shapes, covs, interactions = FALSE, p_threshold = 0.05, model = NULL, force = FALSE,
                       group = NULL, silence = FALSE){
  stopifnot(is.data.frame(covs),is.logical(interactions),is.array(shapes)|is.matrix(shapes), length(dim(shapes)) %in% c(2,3))
  if(!is.null(group)){stopifnot(dim(shapes)[3] == length(group), is.vector(group))}

  if(length(dim(shapes)) == 3)
    shapes = geomorph::two.d.array(shapes)

  # ensure specimen ids exist so within-group results can be restored to input order
  if(is.null(rownames(shapes))){
    if(!is.null(rownames(covs)) && nrow(shapes) == nrow(covs)){
      rownames(shapes) = rownames(covs)
    } else {
      rownames(shapes) = as.character(seq_len(nrow(shapes)))
    }
  }
  input_ids = rownames(shapes)

  # get group index
  if(is.null(group)){
    groups = 1; group = rep(1, nrow(covs))
  } else {
    groups = unique(group)
  }

  results = list()
  for (i in 1:length(groups)){
    ind = as.logical(group == groups[i])
    shape = shapes[ind, , drop = FALSE]
    cov = covs[ind, , drop = FALSE]

    cov = na.omit(cov) # note: if groups have diff covs, this may remove a whole group - do it after filtering
    shape = shape[rownames(cov), , drop = FALSE]

    # evaluate per group when model is not supplied; do not reuse a prior group's fit
    if(is.null(model)){
      this_model = evaluate_model(shapes = shape, covs = cov, interactions = interactions, p_threshold = p_threshold, silence = silence)$final_model
    } else {
      stopifnot(any(class(model) %in% c("procD.lm", "lm.rrpp")))
      this_model = model
    }
    residuals = this_model$residuals

    if(!force && all(this_model$aov.table[setdiff(rownames(this_model$aov.table), c("Residuals", "Total")),"Pr(>F)"] > 0.05)){
      corrected_shapes = shape
      warning("No shape effects in the model were significant. Shapes were not corrected.")
      this_model = "Not corrected"
    } else {
      mean_shape = colMeans(shape)

      # if(length(dim(shape)) == 3){
        mean_shape = matrix(mean_shape, nrow(residuals), ncol(residuals), byrow = TRUE)
      # } else {
      #   residuals = convert(residuals)
      #   mean_shape = array(convert(mean_shape), dim = dim(residuals))}
      corrected_shapes = residuals + mean_shape
      rownames(corrected_shapes) = rownames(shape)
    }
    result = list(corrected_shapes = corrected_shapes, model = this_model)
    results[[paste0("group", i)]] = result
  }

  if(length(groups) > 1){
    results[["groups"]] = data.frame(Assignment = paste0("group", 1:length(groups)), Group = groups)

    # paste corrected shapes back together in the original input order
    corrected_by_group = lapply(results[paste0("group", seq_along(groups))], `[[`, "corrected_shapes")
    corrected_all = do.call(rbind, corrected_by_group)
    retained = input_ids[input_ids %in% rownames(corrected_all)]
    results$corrected_shapes = corrected_all[retained, , drop = FALSE]
  } else {
    results = unlist(results, recursive = FALSE, use.names = TRUE)
    names(results) = gsub(".*\\.", "", names(results))
  }

  class(results) = "residualised"
  return(results)
}

# create custom residualise
print.residualised = function(x, ...){
  require(crayon)
  if(any(grepl("groups", names(x)))){
    groups = names(x)
    groups = groups[grepl("group\\d+", groups)]
    cat("Shapes corrected within groups:","\n")
    # print(x$groups)
    for(i in 1:length(groups)){
      cat(bold(blue("\n   ", x$groups[i,"Group"], "\n")))
      print(summary(x[[groups[i]]][["model"]]))
    }
  } else {
    cat("Shapes corrected using the following model:","\n")
    summary(x$model)
  }
}

