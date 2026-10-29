#' Get ggplot2 theme for consistent plots
#' @import ggplot2

get_theme = function(show_axes = FALSE, legend_size = NULL, axis_title_size = NULL){
  require(ggplot2)
  my_theme = theme_light() + theme(
    axis.title.x = element_text(margin = margin(t = 10), face = "bold", size = axis_title_size),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold", size = axis_title_size),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = legend_size),
    plot.title = element_text(hjust = 0, margin = margin(b = 10)),
    # text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

  )
  if(!show_axes){
    my_theme = my_theme + theme(
      strip.text = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  }

  return(my_theme)
}
