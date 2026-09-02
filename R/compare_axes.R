#' Measure shape vector alignment - INCOMPLETE
#'
#' Each Model's shape vector is defined as the vector between the mutant and control. Alignment is measured by either cosine similarity (default) or correlation.
#'
#' @param shapes currently only compatible with n x (p*k) or PC data (not p x k x n).
#' @param meta data frame. Must contain columns Group (experimental group e.g., control or mutant) and Model (e.g. Hras, Bmp2, Longshanks, etc.). Rownames should be biosample names. Other variables in meta will be used to match for controls if there are insufficient controls within the respective dataset.
#' @param models2compare model(s) to which other models are compared. Must exist in meta$Model.
#' @param method cos_similarity or cor (correlation).
#' @param inclusion_threshold alpha threshold for each Model to determine inclusion in the analysis (Procrustes ANOVA(shapes ~ Group) p-value < inclusion_threshold). If controls and mutants are wildly imbalanced (e.g. more than double the number of controls compared to mutants), available controls are randomly sampled to the same number of mutants.
#' @param group_by_colour groups to split by colour.
#' @param col_manual colours to provide for group_by_colour. If null, use default colours.
#' @param show_legend if TRUE, show legend on right side.
#'
#' @return results data frame with measurement (specified in method) reflecting shape vector alignment of each models2compare to all other models.
#' @return plots bar plots displaying results for each models2compare.
#'
#' @import ggplot2
#' @import Morpho

#TODO: make compatible with 3D format
compare_axes = function(shapes,
                        meta,
                        models2compare,
                        method = c("cos_similarity"),
                        inclusion_threshold = 0.2,
                        group_by_colour = NULL,
                        col_manual = NULL,
                        show_legend = NULL
){
  stopifnot(is.data.frame(meta), all(c("Group", "Model") %in% colnames(meta)), is.matrix(shapes)|is.array(shapes), !is.null(rownames(meta)), is.vector(models2compare))

  model_list = unique(na.omit(meta$Model))
  results = data.frame(matrix(NA, nrow = length(model_list), ncol = length(models2compare)), row.names = model_list)
  colnames(results) = models2compare
  results = cbind(results, p_value = NA)

  control_covs = setdiff(colnames(meta), c("Group", "Model"))
  controls_all = meta[meta$Group == "control",]

  base_axes = matrix(NA,ncol = ncol(shapes),nrow = length(models2compare))
  rownames(base_axes) = models2compare

  for(i in rownames(base_axes)){
    imodel = meta[meta$Model == i,]
    test_gdf = geomorph::geomorph.data.frame(
      shape = shapes[rownames(imodel),],
      group = ifelse(imodel$Group == "control", 0, 1))
    base_axes[i,] = geomorph::procD.lm(shape ~ group, data = test_gdf, iter = 1000, RRPP = TRUE)$coefficients["group",]}

  for (j in 1:length(model_list)){
    jmutant = meta[!is.na(meta$Model) & meta$Model == model_list[j],]

    control_shape = na.omit(rownames(jmutant)[jmutant$Group == "control"])
    mutant_shape = na.omit(rownames(jmutant)[jmutant$Group == "mutant"])

    # TODO: detect a column that might reflect age; if so, match by range
    if(length(control_shape) < 3){
      possible_controls = list()
      for(k in control_covs){possible_controls[[k]] = rownames(controls_all[controls_all[,k] %in% jmutant[,k],])}
      possible_controls <- Reduce(intersect, possible_controls)
      control_shape = unique(c(control_shape, possible_controls))
    }

    control_shape = na.omit(control_shape)
    mutant_shape = na.omit(mutant_shape)

    if(length(control_shape) < 3){cat("\nNot enough suitable controls for",model_list[j])}
    if(length(mutant_shape) < 3){cat("\nNot enough mutants for ", model_list[j])}

    if(length(control_shape) >= 3 && length(mutant_shape) >= 3){
      if (length(mutant_shape)/length(control_shape) < 0.5){
        control_shape = sample(control_shape, size = length(mutant_shape))}

      test_gdf = geomorph::geomorph.data.frame(
        shape = shapes[c(control_shape, mutant_shape),],
        group = c(rep(0, length(control_shape)),rep(1, length(mutant_shape)))  # set 0 to control, 1 to mutant
      )

      test_means = geomorph::procD.lm(shape ~ group, data = test_gdf, iter = 1000, RRPP = TRUE)
      results$p_value[j] = test_means$aov.table["group", "Pr(>F)"]

      if(test_means$aov.table["group", "Pr(>F)"] <= inclusion_threshold){
        # Calculate axes
        jaxis = test_means$coefficients["group",]

        for(i in rownames(base_axes)){
          if(model_list[j] != i){
            if(method == "cor") results[j,i] = cor(jaxis, base_axes[i,])
            if(method == "cos_similarity") results[j,i] = cos(angle.calc(jaxis, base_axes[i,]))
          }
        }
      }
    }
  }

  my_theme <- get_theme()
  col_palette = get_palette()
  plots = list()

  for (i in 1:length(models2compare)){
    iresults = na.omit(results[setdiff(rownames(results),models2compare[i]),models2compare[i], drop = FALSE])
    iresults = cbind(rownames(iresults), iresults); colnames(iresults) = c("Model", "Score")
    if(length(models2compare) > 1) iresults = na.omit(iresults[c(setdiff(rownames(iresults), models2compare), models2compare),])
    if(is.null(group_by_colour)){colour_by_group = ifelse(rownames(iresults) %in% models2compare, "Highlighted", "Other")} else {colour_by_group = group_by_colour}

    ylims = max(abs(iresults$Score))

    plots[[models2compare[i]]] = plotbars(x = reorder(iresults$Model, -abs(iresults$Score)), iresults$Score, group_by_colour = colour_by_group, title = models2compare[i], lim = c(-ylims, ylims), show_legend = show_legend) +
      geom_hline(yintercept = 0, col = my_theme$panel.border@colour, linewidth = 0.3)
  }

  return(list(results = results, plots = plots))
}

