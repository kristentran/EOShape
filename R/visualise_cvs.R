#' Visualise CVA
#'
#' @param cva CVA object from Morpho::CVA
#' @param cvs vector specifying the CVs to visualize. If NULL, all CVs will be visualized.
#' @param group_by_colour groups to split by colour in biplots and density plots (vector of length n).
#' @param col_manual colours for unique values of group_by_colour. If NULL, use default colours.
#' @param group_by_shape groups to split by shape in biplots (vector of length n).
#' @param show_legend if TRUE, show legend on right side.
#' @param mean_shape if provided with atlas_lmks, create mesh3d objects to visualize each CV specified in cvs.
#' @param pca_rotation rotation matrix from PCA to convert back to real shapes for mesh3d objects. Necessary if the data input to CVA were PC scores.
#' @param atlas_mesh if provided with mean_shape, create mesh3d objects to visualize each CV specified in cvs.
#' @param atlas_lmks if provided with mean_shape and atlas_mesh, create mesh3d objects to visualize each CV specified in cvs. Necessary if the atlas_lmks != atlas_mesh vertices.
#'
#' @return plots contains density and biplots for all CVs.
#' @return morphs mesh3d objects of shapes reflecting scores of -3 and +3 along each CV.
#'
#' @import ggplot2
#' @import ggpubr
#'
#' @export

visualise_cvs = function(cva,
                         cvs = NULL,
                         group_by_colour = NULL,
                         col_manual = NULL,
                         group_by_shape = NULL,
                         show_legend = TRUE,
                         mean_shape = NULL,
                         pca_rotation = NULL,
                         atlas_mesh = NULL,
                         atlas_lmks = NULL
){
  if(is.null(cvs)){
    cvs = 1:ncol(cva$CV)
  } else {stopifnot(is.numeric(cvs))}

  evencvs = cvs; if(evencvs[length(evencvs)]%%2 != 0){evencvs = evencvs[evencvs-1]}
  if(length(evencvs) != 0) {xaxes = evencvs[seq(1, length(evencvs), 2)]; yaxes = evencvs[seq(2, length(evencvs), 2)]}

  plots = list()
  morphs = list()

  if(is.null(group_by_colour))
    group_by_colour = cva$groups
  custom_theme = get_theme()

  if(is.null(col_manual)){
    col_manual = get_palette(group_by_colour = group_by_colour)
  } else if(length(col_manual) >= length(unique(group_by_colour))){
    col_manual = col_manual[1:length(unique(group_by_colour))]
  } else {stop("Not enough colours provided for col_manual")}

  if(length(evencvs) != 0) {
    for(i in 1:length(xaxes)){
      xindex = xaxes[i]; yindex = yaxes[i]
      if(!is.na(xindex) & !is.na(yindex)){
        plots[[paste0("CV",xindex,"vs",yindex)]] = biplot(cva$CVscores[,xindex], cva$CVscores[,yindex], group_by_colour = group_by_colour, col_manual = col_manual, group_by_shape = group_by_shape,
                                                          x_title = paste0("CV", xindex), y_title = paste0("CV", yindex))}}}

  for (i in cvs){
    plots[[paste0("CV",i)]] = plotdensity(x = cva$CVscores[,i], group_by_colour = group_by_colour, col_manual = col_manual, title = paste0("CV",i), show_legend = TRUE)

    if(!is.null(atlas_mesh)){
      if(length(levels(cva$groups)) == 2){
        vis_score = sort(c(mean(cva$CVscores[cva$groups == levels(cva$groups)[1]]), mean(cva$CVscores[cva$groups == levels(cva$groups)[2]])))
      } else {
        vis_score = c(-3, 3)*sd(cva$CVscores[,i])
      }
      if(!is.null(ncol(cva$Grandm))){
        CVmin <- vis_score[1]*matrix(cva$CVvis[,i],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm
        CVmax <- vis_score[2]*matrix(cva$CVvis[,i],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm
      } else if(!is.null(mean_shape) & !is.null(pca_rotation)){
        CVmin = vis_score[1]*cva$CVvis[,i] + cva$Grandm
        CVmax = vis_score[2]*cva$CVvis[,i] + cva$Grandm
        CVmin = mean_shape + convert(t(as.matrix(pca_rotation) %*% CVmin))
        CVmax = mean_shape + convert(t(as.matrix(pca_rotation) %*% CVmax))
      }
      if(exists("CVmin") & exists("CVmax")){
        morphs[[paste0("CV",i,"_min")]] = shapeshifter(CVmin, atlas_mesh, atlas_lmks)
        morphs[[paste0("CV",i,"_max")]] = shapeshifter(CVmax, atlas_mesh, atlas_lmks)}
    } else {vis_score = NULL}
  }

  return(list(plots = plots, morphs = morphs, vis_score = vis_score))
}
