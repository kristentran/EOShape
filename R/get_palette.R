#' Get colour palette
#'
#' Get unique colours for ggplot2
#'
#' @param group_by_colour vector of group membership for each individual to determine how many unique groups need colours.
#' @param colorblind if TRUE, choose colours from the TOL palette.
#'
#' @return vector of colours, length of unique group_by_colour values.
#'
#' @export

get_palette = function(group_by_colour = NULL, colorblind = FALSE){
  rainbow = c("#e63232", "#f3722c", "#f8961e", "#ffd043", "#7fc96b", "#43aa8b", "#277da1", "#3b498e", "#66418a", "#a01a7d")
  if(colorblind){rainbow = c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")}

  if(is.null(group_by_colour)){
    col_manual = rainbow
  } else {
    if(length(unique(group_by_colour)) > 5){col_manual = rainbow[1:length(unique(group_by_colour))]} else {col_manual = rainbow[c(1,8,4,6,10)]}
    col_manual = col_manual[1:length(unique(group_by_colour))]}
  return(col_manual)}
