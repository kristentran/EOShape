#' Density plots using ggplot2
#'
#' @param x variable to plot
#' @param group_by_colour groups to split by colour
#' @param col_manual colours to provide for group_by_colour. If null, use default colours.
#' @param lim range of data to plot
#' @param title axis title
#' @param show_legend if TRUE, show legend on right side
#'
#' @return ggplot2 object
#'
#' @import ggplot2
#' @import ggpubr
#' @export

plotdensity = function(x, group_by_colour = NULL, col_manual = NULL, lim = NULL, title = NULL, show_legend = FALSE){
  df = data.frame(x = x)
  custom_theme = get_theme()

  if(is.null(group_by_colour)){
    col_manual = get_palette()[8]
  } else if(is.null(col_manual)){
    col_manual = get_palette(group_by_colour = group_by_colour)
  } else if(length(col_manual) >= length(unique(group_by_colour))){
    col_manual = col_manual[1:length(unique(group_by_colour))]
  } else {stop("Not enough colours provided for col_manual")}

  plot = ggplot(df, aes(x = x, fill = group_by_colour, col = group_by_colour)) + ylab(" ") +
    scale_color_manual(values = alpha(col_manual, 0.2)) + scale_fill_manual(values = alpha(col_manual, 0.5)) +
    custom_theme

  if(is.null(group_by_colour)){
    plot = plot + geom_density(fill = alpha(col_manual, 0.3), col = alpha(col_manual, 0.4))
  } else {plot = plot + geom_density()}

  if(!is.null(lim))
    plot = plot + xlim(lim)
  if(is.null(title))
    plot = plot + theme(axis.title.x = element_blank())
  if(!is.null(title))
    plot = plot + xlab(title)
  if(show_legend)
    plot = plot + theme(legend.position = "right",
                        legend.box.background = element_rect(colour = custom_theme$panel.border@colour),
                        legend.box.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))
  if(!show_legend) # setting for the biplots
    plot = plot + theme(panel.border = element_blank())

  return(plot)
}
