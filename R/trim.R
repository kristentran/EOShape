#' Reduce dimensionality of coordinate data
#'
#' Performs PCA with stats package (prcomp). Reduce the numbers of PCs via Horn's parallel analysis, analysis of the scree plot, or by a variance explained threshold.
#' If method is set to "parallel", then Horn's parallel analysis is performed using PCAtools::parallelPCA. If "scree", then the function will produce a scree plot and prompt the user to select up to which PC to retain. If a numeric is provided, then this will be used as the threshold for cumulative variance explained as a percentage.
#'
#' @param coords aligned coordinate data.
#' @param method chosen method to reduce dimensionality of PCs. See description of function for more details.
#' @param reduce if TRUE, truncate all PCA outputs to the threshold PC.
#' @param threshold PC threshold. If provided, the function will skip the dimensionality reduction method and truncate the PCs directly. Useful if parallel or scree analysis was already performed for your data.
#' @param fromWeb If TRUE, do not return any visualization objects
#' @param atlas_mesh If provided, mesh3d objects will be created to view morphs along the threshold PC.
#' @param atlas_lmks If provided with atlas_mesh, mesh3d objects will be created to view morphs along the threshold PC. Only necessary for mouse data.
#' @param silence if TRUE, suppress printed messages and scree plot display.
#'
#' @return var variance.
#' @return var_explained variance explained.
#' @return var_cumulative cumulative variance explained.
#' @return threshold integer.
#' @return scree plot showing threshold for retaining PCs.
#' @return threshold_min mesh3d object of model for shape along the threshold PC; score of 3 standard deviations below the mean.
#' @return threshold_max mesh3d object of model for shape along the threshold PC; score of 3 standard deviations above the mean.
#'
#' @import ggplot2
#' @import Morpho
#' @import geomorph
#' @import PCAtools
#'
#' @export

trim <- function(
    coords,
    method = "parallel",
    reduce = TRUE,
    threshold = NULL,
    fromWeb = FALSE,
    atlas_mesh = NULL,
    atlas_lmks = NULL,
    silence = FALSE
){

  # if coords is 3d, change to 2d
  if(length(dim(coords)) == 3){coords = two.d.array(coords)}

  # method for dimensionality reduction
  if((grepl("[A-Za-z]", method) & !method %in%  c("parallel", "scree")) | (is.integer(method) & (method<1 | method>100))){
    stop("method must be \"parallel\", \"scree\", or if reducing by cumulative variance, then method must be between 0-100 (%).")
  }

  # run pca
  pca = prcomp(coords)

  # variance
  pca$var = pca$sdev^2
  pca$var_explained = (pca$sdev^2)/sum(pca$sdev^2)
  # cumulative variance
  pca$var_cumulative = cumsum(pca$var_explained)

  rainbow = get_palette()
  custom_theme = get_theme()

  # dimensionality reduction
  if(is.null(threshold)){ # only run if threshold is not provided
    if(method == "parallel"){
      n_pc = PCAtools::parallelPCA(coords, transposed = TRUE)$n
    } else if (method == "scree"){
      n_pc = NA
      # plot scree plot for user and prompt for #
      pcstoplot = min(length(pca$sdev), 70) # maximum 70 PCs to display
      screeplotdata = data.frame(PC = factor(1:pcstoplot, 1:pcstoplot), Variance = pca$var[1:pcstoplot])
      if(!silence)
        print(ggplot(screeplotdata, aes(x = PC, y = Variance)) +
                geom_bar(stat = "identity", width = 0.5, fill = rainbow[8], alpha = 0.8, color = rainbow[8], linewidth = 0.4) +
                custom_theme +
                theme(axis.text = element_text(size = 7)))

      while(!is.integer(n_pc)){
        n_pc = readline(prompt = "Please examine scree plot and enter the number of PCs to retain (integer only): ")
        if(grepl("[0-9]+", n_pc))
          n_pc = as.integer(n_pc)
      }
    } else if (is.numeric(method)){ # Cumulative variance threshold
      n_pc = sum(pca$var_cumulative*100 <= method)
    }
  } else {
    pca$threshold = n_pc = threshold
  }
  pca$threshold = n_pc

  # save a scree plot to pca object
  if(!fromWeb){
    screeplotdata <- data.frame(
      PC = factor(1:(n_pc+10), 1:(n_pc+10)),
      Variance = pca$var[1:(n_pc+10)]
    )
    pca$scree = ggplot(screeplotdata, aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", width = 0.5, fill = rainbow[8], alpha = 0.8, color = rainbow[8], linewidth = 0.4) +
      custom_theme +
      xlab("Principal Components") +
      ylab("Variance Explained") +
      geom_vline(xintercept = n_pc+0.5, linewidth = 1, linetype = 3, col = rainbow[10]) +
      theme(axis.text = element_text(size = 8))
    if(!silence)
      print(pca$scree)
    if(!reduce && !silence)
      cat(blue(bold(n_pc, "PCs to be retained, explaining", round(pca$var_cumulative[n_pc]*100,digits = 1),"% of total variance.")))
  }

  # provide min and max landmark coordinates for the last PC included
  if(!is.null(atlas_mesh)){
    # dummy matrix
    pc_dummies <- matrix(0, nrow = 2, ncol = length(pca$sdev))
    pc_dummies[1,n_pc] <- mean(pca$x[,n_pc]) + pca$sdev[n_pc]*3
    pc_dummies[2,n_pc] <- mean(pca$x[,n_pc]) - pca$sdev[n_pc]*3
    mean_shape = convert(as.vector(suppressWarnings(mshape(coords))))
    PCmin = mean_shape + convert(t(as.matrix(pca$rotation) %*% pc_dummies[1,]))
    PCmax = mean_shape + convert(t(as.matrix(pca$rotation) %*% pc_dummies[2,]))

    if(!fromWeb){
      PCmin = shapeshifter(PCmin, atlas_mesh, atlas_lmks)
      PCmax = shapeshifter(PCmax, atlas_mesh, atlas_lmks)
    }
    pca[[paste0("threshold_min")]] = PCmin
    pca[[paste0("threshold_max")]] = PCmax
  }

  if(reduce == TRUE){
    if(!silence)
      cat(blue(bold("Retaining", n_pc, "PCs, explaining", round(pca$var_cumulative[n_pc]*100,digits = 1),"% of total variance.")))
    pca[c("x", "rotation")] <- lapply(pca[c("x", "rotation")], \(m) m[, 1:n_pc, drop = FALSE])
    pca[c("sdev", "var", "var_cumulative")] <- lapply(pca[c("sdev", "var", "var_cumulative")], \(v) v[1:n_pc])}

  return(pca)
}
