#' Find differences
#'
#' Are there differences in csize and shape by independant variable(s)? Uses
#' geomorph::procD.lm as the workhorse. If only one independant variable is
#' provided and it is categorical, then post-hoc testing will be completed (if
#' that variable has more than 2 levels). * note that this is an analysis fx, so
#' the object returned are just results, not the borealis object
#'
#' TODO: get landmark coordinates for HOW they differ if significant. If
#' continuous, morph along axis from procd.lm obj. If categorical, CVA?
#' TODO: If >1 independant var, use evaluate model to get assess sig. variables.
#' Then only model the significant terms
#'
#' If differences exist, then how are they different? If continuous, then
#'
#' @param borealis borealis object after metadata detection
#' @param independant default to 'All' specified independant variables. Specify variable if only one is desired.
#' @param interactions if TRUE and there are multiple independant variables, ...
#' @param silence if TRUE, suppress printed messages from process_shapes and internal calls.
#'
#' @return borealis object for subsequent analyses
#'
#' @export

find_diffs = function(borealis, independant = "All", interactions = TRUE, silence = FALSE, ...){
  if(!inherits(borealis, "borealis"))
    stop("Object provided is not a borealis object. Please create using initialise_borealis().")
  if(!is.logical(silence) || length(silence) > 1)
    stop("silence parameters must be single value logicals")
  if(independant != "All"){
    independants = gsub("independant_", "", colnames(borealis$meta)[grepl("independant", colnames(borealis$meta))])
    if(!any(gsub("independant_", "", independant) %in% independants))
      stop("independant argument should be 'All' or the name of variable to analyse (e.g. Experimental_Group, or Treatment_Dosage).")
    rm(independants)
  }

  # get column names from meta data for independant variables
  if(independant == "All"){
    independants = colnames(borealis$meta)[grepl("independant", colnames(borealis$meta))]
  } else {
    if(!grepl("independant_", independant))
      independant = paste0("indepenant_", independant)
    independants = independant
  }

  # get data frame of only independant variables
  independant_df = borealis$meta[,independants, drop = FALSE]

  # analyze shape and size differences
  dependants = c("shapes", "size")

  # object to store results
  analysis_results = list()

  # loop through size and shape as Y
  for(dependant in dependants){
    # get right side formula
    ind_terms = colnames(independant_df)
    if(dependant != "size" && borealis$analysis_params$assess_allometry)
      ind_terms = c(ind_terms, "size")
    if(interactions & length(ind_terms) > 1){int_terms = apply(combn(ind_terms, 2), 2, paste, collapse = ":")
    } else {int_terms = c()}
    rightside_formula = paste0(c(ind_terms, int_terms), collapse = " + ")

    lm_df = geomorph::geomorph.data.frame(
      y = borealis[[dependant]],
      independant_df,
      size = borealis$size
    )

    aov = geomorph::procD.lm(as.formula(paste0("y ~ ", rightside_formula)), data = lm_df)$aov
    analysis_results[[dependant]][["full_aov"]] = aov

    # post hoc testing if multiple levels in group; only available if independant != all and the independant variable is categorical
    if(ncol(independant_df) == 1){ # proceed if only one independant var
      if(is.factor(independant_df[,1]) || is.character(independant_df[,1])){ # if independant var is categorical
        n_tests = combn(unique(independant_df[,1]), 2)

        if(ncol(n_tests) > 1 && aov[1, "Pr(>F)"] < 0.05){
          vis_diffs = list()
          posthoc = data.frame(Group = rep(NA, ncol(n_tests)), Df = NA, SumSq = NA, MeanSq = NA, Fvalue = NA, Zvalue = NA, Pvalue = NA)

          for(i in 1:ncol(n_tests)){
            index = independant_df[,1] %in% n_tests[,i]
            lm_df = data.frame(
              y = borealis[[dependant]],
              group = independant_df[,1][index]
            )
            aov_posthoc = geomorph::procD.lm(y ~ group, data = lm_df)$aov
            posthoc[i,] = c(paste0(n_tests[,i], collapse = " vs "),
                            aov_posthoc[1, c("Df", "SS", "MS", "F", "Z", "Pr(>F)")])
          }
          analysis_results[[dependant]][["post_hoc"]] = posthoc
        }

        # perform CVA if on shape analysis and if the full model was significant (and also if only one independant categorical var)
        if(dependant == "shapes" && aov[1, "Pr(>F)"] < 0.05){
          cva = Morpho::CVA(borealis[[dependant]], groups = independant_df[,1])
          if(borealis$analysis_params$use_PC_scores){
            mean_shape = borealis$mshape
            rotation = borealis$pca$rotation
          } else {mean_shape = rotation = NULL}
          vis_diffs = abind::abind(visualise_cvs(cva, mean_shape = mean_shape, pca_rotation = rotation, fromWeb = TRUE)$morphs, along =3)
        }
      } else if(is.numeric(independant_df[,1])){ # if independant var is numeric
        # visualise along the axis
        if(dependant == "shapes" && aov[1, "Pr(>F)"] < 0.05){
          vis_diffs = list()
          coef = geomorph::procD.lm(as.formula(paste0("y ~ ", colnames(independant_df))), data = lm_df, iter = 1)$coefficients
          selected.scores <- c(mean(independant_df[,1]) - (sd(independant_df[,1])*3), mean(independant_df[,1]) + (sd(independant_df[,1])*3))
          names(selected.scores) = c("neg", "pos")

          for(i in 1:length(selected.scores)){
            datamod <- ~ selected.scores[i]
            mat <- model.matrix(datamod)
            pred <- mat %*% coef

            if(borealis$analysis_params$use_PC_scores){
              pred = borealis$mshape + convert(t(as.matrix(borealis$pca$rotation) %*% pred[1,]))
            }
            vis_diffs[[names(selected.scores)[i]]] = pred
          }
        }
      }
    }
  }
  return(list(anovas = analysis_results, landmark_diffs = vis_diffs))
}
