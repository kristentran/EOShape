#' Biplot and density plots using ggplot2
#'
#' @param x data to plot on x axis.
#' @param y data to plot on y axis.
#' @param group_by_colour groups to split by colour (vector of length n).
#' @param col_manual colours for unique values of group_by_colour. If NULL, use default colours.
#' @param group_by_shape groups to split by shape (vector of length n).
#' @param x_title x axis title.
#' @param y_title y axis title,
#' @param show_legend if TRUE, show legend on right side.
#'
#' @return ggplot2 object
#'
#' @export

# TO DO show_legend

plotscatter = function(x, y, group_by_colour = NULL, col_manual = NULL, group_by_shape = NULL, x_title = NULL, y_title = NULL, show_legend = TRUE, show_axes = TRUE){
  require(ggplot2)
  # create biplots - maybe have this as a separate function
  df = data.frame(x = x, y = y)
  custom_theme = get_theme(show_axes = show_axes); shapes_manual = c(16,17,15,18,5,8,3,4,2,0)

  if(is.null(group_by_colour)){
    col_manual = get_palette()[8]
  } else if(is.null(col_manual)){
    col_manual = get_palette(group_by_colour = group_by_colour)
  } else if(length(col_manual) >= length(unique(group_by_colour))){
    col_manual = col_manual[1:length(unique(group_by_colour))]
  } else {stop("Not enough colours provided for col_manual")}

  if(!is.null(group_by_shape)){
    if(length(unique(group_by_shape)) > 10){stop("Too many groups for number of available shapes")}
    shapes_manual = shapes_manual[1:length(unique(group_by_shape))]}

  if(is.null(x_title))
    x_title = " "
  if(is.null(y_title))
    y_title = " "

  scatter = ggplot(df, aes(x = x, y = y, col = group_by_colour, fill = group_by_colour, shape = group_by_shape)) +
    scale_color_manual(values = col_manual) + scale_fill_manual(values = col_manual) + scale_shape_manual(values = shapes_manual) +
    xlab(x_title) + ylab(y_title) + custom_theme +
    theme(legend.position = "right",
          legend.box.background = element_rect(colour = custom_theme$panel.border@colour),
          legend.box.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))

  if(is.null(group_by_colour)){
    scatter = scatter + geom_point(alpha = 0.7, col = col_manual)
  } else {
    scatter = scatter + geom_point(alpha = 0.7) #+ stat_ellipse(type = "norm", geom = "polygon", alpha = 0.1, col = NA)
  }

  return(scatter)
}
