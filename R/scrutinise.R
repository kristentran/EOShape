#' Scrutinise differences in shape or size
#'
#' Are there differences in shape (or size) by independent variable(s)? Uses
#' geomorph::procD.lm as the workhorse. If only one independent variable is
#' provided, then:
#' \itemize{
#'   \item if it is categorical with more than 2 levels and the full model is
#'   significant, pairwise post-hoc Procrustes ANOVAs are performed.
#'   \item if it is categorical, \code{y} is multivariate, and the full model is
#'   significant, a CVA (Morpho::CVA) is performed and shapes along each CV
#'   (or group means if there are 2 groups) are returned.
#'   \item if it is numeric, \code{y} is multivariate, and the full model is
#'   significant, predicted shapes at -3 and +3 standard deviations from the
#'   mean of the independent variable are returned.
#' }
#'
#' @param y dependent variable. Shapes as a 2d matrix (coordinates or PC
#'   scores) or 3d array, a PCA object from \code{trim()} or \code{prcomp()}
#'   (PC scores are used; \code{pca_rotation} and \code{mean_shape} are taken
#'   from the object if not provided), or a numeric vector (e.g. centroid size).
#' @param independents data frame of independent variables (columns), or a
#'   vector for a single independent variable.
#' @param covs optional data frame of additional covariates to include in the
#'   model (e.g. size to assess allometry). Covariates are not used for
#'   post-hoc tests or to visualise differences.
#' @param interactions if TRUE and there are multiple terms, include all
#'   pairwise interactions.
#' @param p_threshold numerical between 0 and 1. Significance required before
#'   post-hoc testing and visualising differences.
#' @param mean_shape mean shape (p x k matrix) corresponding to
#'   \code{pca_rotation}. Provide with \code{pca_rotation} if \code{y} are PC
#'   scores so that landmark_diffs are returned as landmarks. If not provided,
#'   multivariate \code{y} is assumed to be landmark coordinates.
#' @param pca_rotation rotation matrix from the PCA used to produce \code{y}.
#'   Provide with \code{mean_shape}.
#' @param iter number of permutations passed to geomorph::procD.lm.
#' @param silence if TRUE, suppress printed progress from procD.lm.
#'
#' @return full_aov ANOVA table of the full model.
#' @return post_hoc data frame of pairwise ANOVAs, or NULL if not performed.
#' @return cva Morpho::CVA object, or NULL if not performed.
#' @return landmark_diffs 3d array of CV morphs (categorical) or a list of
#'   \code{neg} and \code{pos} shapes (numeric); NULL if not performed.
#'
#' @import geomorph
#'
#' @export

scrutinise = function(y, independents, covs = NULL, interactions = TRUE, p_threshold = 0.05,
                      mean_shape = NULL, pca_rotation = NULL, iter = 999, silence = FALSE){
  if(!is.data.frame(independents))
    independents = data.frame(group = independents)
  if(!is.null(covs)) stopifnot(is.data.frame(covs))

  # PCA object (e.g. from trim): use PC scores; prcomp's center is the mean shape
  if(inherits(y, "prcomp")){
    if(is.null(pca_rotation)) pca_rotation = y$rotation
    if(is.null(mean_shape) && !isFALSE(y$center)) mean_shape = y$center
    y = y$x
  }
  stopifnot(is.numeric(y), is.logical(interactions), is.numeric(p_threshold))
  if(!is.logical(silence) || length(silence) > 1)
    stop("silence parameters must be single value logicals")
  if(xor(is.null(mean_shape), is.null(pca_rotation)))
    stop("mean_shape and pca_rotation must be provided together.")

  if(length(dim(y)) == 3)
    y = geomorph::two.d.array(y)
  multivariate = !is.null(dim(y)) && ncol(y) > 1
  n = if(multivariate) nrow(y) else length(y)
  if(nrow(independents) != n || (!is.null(covs) && nrow(covs) != n))
    stop("y, independents, and covs must have the same number of specimens.")
  pc_scores = multivariate && !is.null(pca_rotation)
  if(pc_scores){
    mean_shape = convert(mean_shape)
    if(ncol(y) > ncol(pca_rotation))
      stop("y has more PCs than pca_rotation.")
    pca_rotation = as.matrix(pca_rotation)[, seq_len(ncol(y)), drop = FALSE]
  }
  # landmarks can only be visualised from coordinates (p x 3) or from PC scores with their rotation
  can_visualise = multivariate && (pc_scores || ncol(y) %% 3 == 0)
  if(multivariate && !can_visualise)
    warning("y does not look like 3D landmark coordinates; if y are PC scores, provide mean_shape and pca_rotation to visualise differences.")

  # right side formula
  ind_terms = c(colnames(independents), colnames(covs))
  if(interactions & length(ind_terms) > 1){int_terms = apply(combn(ind_terms, 2), 2, paste, collapse = ":")
  } else {int_terms = c()}
  rightside_formula = paste0(c(ind_terms, int_terms), collapse = " + ")

  if(is.null(covs)) model_df = independents else model_df = cbind(independents, covs)
  lm_df = geomorph::geomorph.data.frame(y = y, model_df)
  aov = geomorph::procD.lm(as.formula(paste0("y ~ ", rightside_formula)), data = lm_df,
                           iter = iter, print.progress = !silence)$aov.table
  significant = aov[1, "Pr(>F)"] < p_threshold

  post_hoc = cva = landmark_diffs = NULL

  # post hoc testing and visualising differences only if there is one independent variable
  if(ncol(independents) == 1 && significant){
    independent = independents[,1]

    if(is.factor(independent) || is.character(independent)){ # if independent var is categorical
      n_tests = combn(as.character(unique(independent)), 2)

      if(ncol(n_tests) > 1){
        post_hoc = data.frame()
        for(i in 1:ncol(n_tests)){
          index = independent %in% n_tests[,i]
          posthoc_df = geomorph::geomorph.data.frame(
            y = if(multivariate) y[index, , drop = FALSE] else y[index],
            group = independent[index]
          )
          aov_posthoc = geomorph::procD.lm(y ~ group, data = posthoc_df, iter = iter, print.progress = !silence)$aov.table
          post_hoc = rbind(post_hoc, data.frame(
            Group = paste0(n_tests[,i], collapse = " vs "),
            Df = aov_posthoc[1, "Df"], SumSq = aov_posthoc[1, "SS"], MeanSq = aov_posthoc[1, "MS"],
            Fvalue = aov_posthoc[1, "F"], Zvalue = aov_posthoc[1, "Z"], Pvalue = aov_posthoc[1, "Pr(>F)"]
          ))
        }
      }

      # CVA to see how groups differ
      if(multivariate){
        cva = Morpho::CVA(y, groups = independent)
        if(can_visualise){
          morphs = visualise_cvs(cva, mean_shape = mean_shape, pca_rotation = pca_rotation, fromWeb = TRUE)$morphs
          landmark_diffs = abind::abind(morphs, along = 3)
        }
      }
    } else if(is.numeric(independent) && can_visualise){ # if independent var is numeric, visualise along the axis
      coef = geomorph::procD.lm(y ~ independent, data = geomorph::geomorph.data.frame(y = y, independent = independent),
                                iter = 1, print.progress = FALSE)$coefficients
      selected_scores = mean(independent) + c(neg = -3, pos = 3)*sd(independent)

      landmark_diffs = list()
      for(i in names(selected_scores)){
        pred = c(1, selected_scores[[i]]) %*% coef
        if(pc_scores){
          pred = mean_shape + convert(t(pca_rotation %*% pred[1,]))
        } else {
          pred = convert(pred)
        }
        landmark_diffs[[i]] = pred
      }
    }
  }

  return(list(full_aov = aov, post_hoc = post_hoc, cva = cva, landmark_diffs = landmark_diffs))
}
