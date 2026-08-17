#' Evaluate linear models
#'
#' Model shapes ~ covariates using procD.lm. Identify which covariates are significant and provide final model.
#' If \code{force = TRUE}, skip term selection and use the full covariate set.
#'
#' @param shapes 2d matrix or 3d array of shape variables.
#' @param covs data frame of covariates.
#' @param interactions if TRUE, include all possible pairwise interactions.
#' @param p_threshold numerical between 0 and 1. Required significance to retain covariates.
#' @param silence if TRUE, suppress printed model evaluation messages.
#' @param force if TRUE, skip stepwise term selection and return the starting
#'   model (all covariates, and pairwise interactions if \code{interactions = TRUE})
#'   as the final model.
#'
#' @return starting_model procD.lm object with all covariates provided.
#' @return final_model procD.lm object with only significant covariates, or the
#'   starting model when \code{force = TRUE}.
#'
#' @import geomorph
#'
#' @export

# shapes = data$lmks
# covs = data$meta[,"sex", drop= FALSE]
# interactions = TRUE
# p_threshold = 0.05

evaluate_model = function(shapes, covs, interactions = TRUE, p_threshold = 0.05, silence = FALSE, force = FALSE){
  gdf = geomorph.data.frame(shapes = shapes, covs)
  ind_terms = colnames(covs)
  if(interactions & length(ind_terms) > 1){int_terms = apply(combn(ind_terms, 2), 2, paste, collapse = ":")
  } else {int_terms = c()}
  starting_formula = paste0(c(ind_terms, int_terms), collapse = " + ")
  starting_model = procD.lm(as.formula(paste0("shapes ~ ", starting_formula)), data = gdf, print.progress = !silence)
  if(!silence)
    cat(blue(bold("Starting model:\n")), "shapes ~", starting_formula, "\n")

  if(force){
    if(!silence)
      cat(blue(bold("Force = TRUE: skipping term selection; using starting model as final model.\n")))
    return(list(final_model = starting_model, starting_model = starting_model))
  }

  model = starting_model$aov.table; mod_ind_terms = ind_terms; mod_int_terms = int_terms
  last_formula = starting_formula; new_formula = c("empty")
  max_iter = length(c(ind_terms, int_terms));iter = 0

  if(all(model[setdiff(rownames(model), c("Residuals", "Total")),"Pr(>F)"] > 0.05))
    return(list(final_model = starting_model, starting_model = starting_model))

  while(last_formula != new_formula && max_iter >= iter){
    last_formula = new_formula
    iter = iter + 1
    mod_int_terms = na.omit(mod_int_terms[ifelse(model[mod_int_terms,"Pr(>F)"] <= p_threshold, TRUE, FALSE)])

    if(interactions){
      index = c()
      for(i in mod_ind_terms){
        if(model[i,"Pr(>F)"] > p_threshold & !any(grepl(i, mod_int_terms))){index = c(index, FALSE)} else {index = c(index, TRUE)}
      }; mod_ind_terms = mod_ind_terms[index]
    } else {mod_ind_terms = na.omit(mod_ind_terms[ifelse(model[mod_ind_terms,"Pr(>F)"] <= p_threshold, TRUE, FALSE)])}

    new_formula = paste0(c(mod_ind_terms, mod_int_terms), collapse = " + ")
    if(last_formula != new_formula){
      if(!silence)
        cat("\rEvaluating model", rep(".", iter))
      model = procD.lm(as.formula(paste0("shapes ~ ", paste0(c(mod_ind_terms, mod_int_terms), collapse = " + "))), data = gdf, print.progress = !silence)$aov.table
    }
  }
  if(max_iter == iter){warning("\nEntered infinite loop; please review model manually.")}
  final_model = procD.lm(as.formula(paste0("shapes ~ ", paste0(c(mod_ind_terms, mod_int_terms), collapse = " + "))), data = gdf, print.progress = !silence)
  if(!silence)
    cat(blue(bold("\rEvaluation complete! Final model:\n")), "shapes ~", new_formula, "\n")
  return(list(final_model = final_model, starting_model = starting_model))
}
