#' Barplots using ggplot2
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
#' @export

plotbars = function(x, y = NULL, group_by_colour = NULL, col_manual = NULL, lim = NULL, title = NULL, show_legend = FALSE){
  require(EOS)
  if(!is.null(y)){
    df = data.frame(x = x, y = y)
  } else {
    df = data.frame(x = x)
  }

  custom_theme = get_theme()

  if(is.null(group_by_colour)){
    col_manual = get_palette()[8]
  } else if(is.null(col_manual)){
    col_manual = get_palette(group_by_colour = group_by_colour)
  } else if(length(col_manual) > length(unique(group_by_colour))){
    col_manual = col_manual[1:length(unique(group_by_colour))]
  } else if(length(col_manual) < length(unique(group_by_colour))){stop("Not enough colours provided for col_manual")}

  if(!is.null(y)){
    plot = ggplot(df, aes(x = x, y = y, fill = group_by_colour, col = group_by_colour)) +
    geom_col(width = 0.5, alpha = 0.7, linewidth = 0.3) +
    scale_color_manual(values = col_manual) + scale_fill_manual(values = col_manual) +
    custom_theme + theme(axis.title.x = element_blank(),
                         axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
                         plot.margin = unit(c(1,2,1,2), "mm"))
  } else {
    plot = ggplot(df, aes(x = x, fill = group_by_colour, col = group_by_colour)) +
      geom_bar(width = 0.5, alpha = 0.7, linewidth = 0.3) +
      scale_color_manual(values = col_manual) + scale_fill_manual(values = col_manual) +
      custom_theme + theme(axis.title.x = element_blank(),
                           axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
                           plot.margin = unit(c(1,2,1,2), "mm"))
    }

  if(!is.null(lim))
    plot = plot + ylim(lim)
  if(is.null(title))
    plot = plot + theme(axis.title.y = element_blank())
  if(!is.null(title))
    plot = plot + ylab(title)
  if(show_legend)
    plot = plot + theme(legend.position = "right",
                        legend.box.background = element_rect(colour = custom_theme$panel.border@colour, linewidth = 0.3),
                        legend.box.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))

  return(plot)
}
