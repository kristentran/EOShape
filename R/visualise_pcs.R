#' Visualise CVA
#'
#' @param pca pca object from prcomp or reduce_shapes.
#' @param pcs vector specifying the PCs to visualize. If NULL, all PCs will be visualized.
#' @param group_by_colour groups to split by colour in biplots and density plots (vector of length n).
#' @param col_manual colours for unique values of group_by_colour. If NULL, use default colours.
#' @param group_by_shape groups to split by shape in biplots (vector of length n).
#' @param show_legend if TRUE, show legend on right side.
#' @param mean_shape if provided with atlas_lmks, create mesh3d objects to visualize each PC specified in pcs.
#' @param atlas_mesh if provided with mean_shape, create mesh3d objects to visualize each PC specified in pcs.
#' @param atlas_lmks if provided with mean_shape and atlas_mesh, create mesh3d objects to visualize each PC specified in pcs. Necessary if the atlas_lmks != atlas_mesh vertices.
#'
#' @return plots contains density and biplots for all PCs.
#' @return morphs mesh3d objects of shapes reflecting scores of -3 and +3 along each PC.
#'
#' @import ggplot2
#' @import ggpubr
#'
#' @export

visualise_pcs = function(pca,
                         pcs = 1:2,
                         group_by_colour = NULL,
                         col_manual = NULL,
                         group_by_shape = NULL,
                         show_legend = TRUE,
                         mean_shape = NULL,
                         atlas_mesh = NULL,
                         atlas_lmks = NULL
){
  stopifnot(is.numeric(pcs), length(pcs) >= 2)

  evenpcs = pcs; if(evenpcs[length(evenpcs)]%%2 != 0){evenpcs = evenpcs[evenpcs-1]}
  xaxes = evenpcs[seq(1, length(evenpcs), 2)]; yaxes = evenpcs[seq(2, length(evenpcs), 2)]
  plots = list()

  custom_theme = get_theme(); shapes_manual = c(16,17,18,15,5,8,3,4,2,0)
  if(is.null(group_by_colour)){
    col_manual = get_palette()[1]
  } else if(is.null(col_manual)){
    col_manual = get_palette(group_by_colour = group_by_colour)
  } else if(length(col_manual) >= length(unique(group_by_colour))){
    col_manual = col_manual[1:length(unique(group_by_colour))]
  } else {stop("Not enough colours provided for col_manual")}

  # biplots for every two pcs
  for(i in 1:length(xaxes)){
    xindex = xaxes[i]; yindex = yaxes[i]
    xvar = round(pca$sdev[xindex]^2/sum(pca$sdev^2)*100, 1)
    yvar = round(pca$sdev[yindex]^2/sum(pca$sdev^2)*100, 1)

    plots[[paste0("PC",xindex,"vs",yindex)]] = suppressWarnings(biplot(pca$x[,xindex], pca$x[,yindex], group_by_colour = group_by_colour, col_manual = col_manual, group_by_shape = group_by_shape,
                                                                       x_title = paste0("PC", xindex, " (", xvar,"%)"), y_title = paste0("PC", yindex, " (", yvar,"%)")))}
  for (i in pcs){
    plots[[paste0("PC",i)]] = plotdensity(pca$x[,i], group_by_colour = group_by_colour, col_manual = col_manual, title = paste0("PC",i, " (", round(pca$sdev[i]^2/sum(pca$sdev^2)*100, 1), "%)"), show_legend = TRUE)}

  # PC morphs
  morphs = list()
  if(!is.null(atlas_mesh) & !is.null(mean_shape)){
    mean_shape = convert(mean_shape)
    for(i in pcs){
      pc_dummies <- matrix(0, nrow = 2, ncol = length(pca$sdev))
      pc_dummies[1,i] <- mean(pca$x[,i]) + pca$sdev[i]*3
      pc_dummies[2,i] <- mean(pca$x[,i]) - pca$sdev[i]*3
      PCmax = mean_shape + convert(t(as.matrix(pca$rotation) %*% pc_dummies[1,]))
      PCmin = mean_shape + convert(t(as.matrix(pca$rotation) %*% pc_dummies[2,]))
      morphs[[paste0("PC",i,"_min")]] = shapeshifter(PCmin, atlas_mesh, atlas_lmks)
      morphs[[paste0("PC",i,"_max")]] = shapeshifter(PCmax, atlas_mesh, atlas_lmks)
    }
  }
  return(list(plots = plots, morphs = morphs))
}
