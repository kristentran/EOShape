#' Visualise CVA
#'
#' @param cva CVA object from Morpho::CVA
#' @param cvs vector specifying the CVs to visualise. If NULL, all CVs will be visualised.
#' @param group_by_colour groups to split by colour in biplots and density plots (vector of length n).
#' @param col_manual colours for unique values of group_by_colour. If NULL, use default colours.
#' @param group_by_shape groups to split by shape in biplots (vector of length n).
#' @param show_legend if TRUE, show legend on right side.
#' @param mean_shape mean shape corresponding to PCA rotation. Provide
#' @param pca_rotation rotation matrix from PCA to convert back to real shapes for mesh3d objects. Necessary if the data input to CVA were PC scores.
#' @param atlas_mesh if provided with mean_shape, create mesh3d objects to visualise each CV specified in cvs.
#' @param atlas_lmks if provided with mean_shape and atlas_mesh, create mesh3d objects to visualise each CV specified in cvs. Necessary if the atlas_lmks != atlas_mesh vertices.
#'
#' @return plots contains density and biplots for all CVs.
#' @return morphs mesh3d objects of shapes reflecting scores of -3 and +3 along each CV.
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
                         atlas_lmks = NULL,
                         fromWeb = FALSE
){
  if(is.null(cvs)){
    cvs = 1:ncol(cva$CV)
  } else {stopifnot(is.numeric(cvs))}

  evencvs = cvs; if(evencvs[length(evencvs)]%%2 != 0){evencvs = evencvs[evencvs-1]}
  if(length(evencvs) != 0) {xaxes = evencvs[seq(1, length(evencvs), 2)]; yaxes = evencvs[seq(2, length(evencvs), 2)]}

  plots = list()
  morphs = list()

  # Only create ggplots if not from tool
  if(!fromWeb){
    if(is.null(group_by_colour))
      group_by_colour = cva$groups
    custom_theme = get_theme()

    if(is.null(col_manual)){
      col_manual = get_palette(group_by_colour = group_by_colour)
    } else if(length(col_manual) >= length(unique(group_by_colour))){
      col_manual = col_manual[1:length(unique(group_by_colour))]
    } else {stop("Not enough colours provided for col_manual")}

    # create biplots for CVs
    if(length(evencvs) != 0) {
      for(i in 1:length(xaxes)){
        xindex = xaxes[i]; yindex = yaxes[i]
        if(!is.na(xindex) & !is.na(yindex)){
          plots[[paste0("CV",xindex,"vs",yindex)]] = biplot(cva$CVscores[,xindex], cva$CVscores[,yindex], group_by_colour = group_by_colour, col_manual = col_manual, group_by_shape = group_by_shape,
                                                            x_title = paste0("CV", xindex), y_title = paste0("CV", yindex))}}}
  }
  for (i in cvs){
    # Create CV density plots if not from tool
    if(!fromWeb)
      plots[[paste0("CV",i)]] = plotdensity(x = cva$CVscores[,i], group_by_colour = group_by_colour, col_manual = col_manual, title = paste0("CV",i), show_legend = TRUE)

    # CALCULATING SHAPES
    # Decide scores to visualise: if there are 2 groups, visualise group mean scores, else +/- 3SD
    if(length(levels(cva$groups)) == 2){ # 2 levels = 1 CV therefore it doesn't matter that groupmeans is all CVs compared to per CV in the else condition.
      # save mean CV score by group to vis_score
      vis_score = tapply(cva$CVscores[, i], cva$groups, mean)
      # group with the lower mean CV score is "min", the higher is "max"
      ordered_groups = names(sort(vis_score))
      for(extreme in c("min", "max")){
        group = ordered_groups[ifelse(extreme == "min", 1, 2)]
        tmpshape = cva$groupmeans[group,]
        if(!is.null(mean_shape) & !is.null(pca_rotation)){
          # if pc scores provided
          tmpshape = mean_shape + convert(t(as.matrix(pca_rotation) %*% tmpshape))
        } else {
          tmpshape = convert(tmpshape)
        }
        if(!is.null(atlas_mesh)){
          tmpshape = shapeshifter(tmpshape, atlas_mesh, atlas_lmks)
        }
        morphs[[paste0("CV",i,"_",extreme)]] = tmpshape
      }
    } else {
      vis_score = c(-3, 3)*sd(cva$CVscores[,i]) + mean(cva$CVscores[,i])

      # Calculate shapes from vis_scores:
      for(extreme in c("min", "max")){
        vis_ind = ifelse(extreme == "min", 1, 2)
        if(!is.null(mean_shape) & !is.null(pca_rotation)){
          # if pc scores provided
          CVextreme = vis_score[vis_ind]*cva$CVvis[,i] + cva$Grandm
          CVextreme = mean_shape + convert(t(as.matrix(pca_rotation) %*% CVextreme))
        } else {
          # if coordinates provided
          if(!is.null(ncol(cva$Grandm))){ # if shapes in 3d format
            CVextreme = vis_score[vis_ind]*matrix(cva$CVvis[,i],nrow(cva$Grandm),ncol(cva$Grandm))+cva$Grandm
          } else { #if shapes in 2d format
            CVextreme = matrix((vis_score[vis_ind]*cva$CVvis[,i]+cva$Grandm), nrow = length(cva$Grandm)/3, ncol = 3, byrow = TRUE)
          }
        }

        if(exists("CVextreme")){
          if(!is.null(atlas_mesh)){
            CVextreme = shapeshifter(CVextreme, atlas_mesh, atlas_lmks)
            # CVmax = shapeshifter(CVmax, atlas_mesh, atlas_lmks)
          }
          # store morphs as landmarks or meshes
          morphs[[paste0("CV",i,"_",extreme)]] = CVextreme
          rm(CVextreme)
        }
      } # end of max min loop
    }
  } # end of CV loop
  return(list(plots = plots, morphs = morphs, vis_score = vis_score))
}
