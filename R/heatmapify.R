#' Wrapper for meshDist from Morpho
#'
#' Use either Blue-White-Red or Green-White-Magenta colors to show 3D shape effects.
#'
#' @param x mesh3D object to be morphed.
#' @param y mesh3D object to show changes from x.
#' @param direction "normal" or "reverse"; if normal, heatmap shows changes from x to y (default meshDist). If reversed, heatmap shows changes from y to x.
#' @param palette "BluRed" or "GreMag" for heatmap palette.
#'
#' @return meshDist object
#'
#' @import Morpho
#'
#' @export

heatmapify = function(x,y, direction = "normal", palette = "BluRed", dists = NULL){
  require(Morpho)

  # check if arguments are the correct class
  stopifnot(class(x)[1] == "mesh3d", class(y)[1] == "mesh3d",
            is.character(direction),
            is.character(palette),
            is.null(dists) || is.numeric(dists))

  # make sure arguments are the correct choices
  direction <- match.arg(direction, c("normal", "reverse"))
  palette <- match.arg(palette, c("GreMag", "BluRed"))

  # Create heatmap distances, if not provided - this is needed to ensure that the heatmap is symmetric
  if(is.null(dists)){
    dists <- max(abs(meshDist(x,y)$dists))}

  # Choose colour palette: Green and Magenta | Blue and Red
  if (palette == "GreMag"){
    # Save two different colour schemes
    coloursofthewind <- colorRampPalette(colors = c(
      "#300033", "#450e49", "#591e5d", "#6c2a71", "#7e3584", "#8f4096", "#a04aa8", "#b054b9", "#bf5dca", "#ce66da", "#dd6eeb", "#eb77fa", "#ed8cf8", "#eca1f4", "#ecb3f1", "#edc2f0", "#eecef1", "#f0d9f2", "#f3e3f4", "#f6ebf7", "#f8f2f9", "#fcfbfc", "#ffffff", "#f5f9f9", "#e5f1f1", "#d9edec", "#cde8e7", "#bce3e1", "#abdddc", "#94d7d5", "#77d1cf", "#4bcbc8", "#2cc1bf", "#2db5b3", "#2da8a5", "#2c9a98", "#2a8c8a", "#277d7c", "#246f6e", "#1f5f5e", "#194f4f", "#123e3e", "#082c2c", "#001919"
    ))(100); coloursofthewind <- rev(coloursofthewind) # green should be inward vs pink should be outward
  } else if(palette == "BluRed") {
    coloursofthewind <- c("#053061", "#083669", "#0B3C72", "#0E427A", "#114882", "#154E8B", "#185493", "#1B5A9B", "#1E60A4", "#2166AC", "#2166AC", "#256BAF", "#2970B1", "#2C75B4", "#307AB6", "#347FB9", "#3884BB", "#3B89BE", "#3F8EC0", "#4393C3", "#4393C3", "#4C99C6", "#559EC9", "#5DA4CC", "#66A9CF", "#6FAFD2", "#78B4D5", "#80BAD8", "#89BFDB", "#92C5DE", "#92C5DE", "#99C9E0", "#A0CCE2", "#A7D0E4", "#AED3E6", "#B5D7E8", "#BCDAEA", "#C3DEEC", "#CAE1EE", "#D1E5F0", "#D1E5F0", "#D5E7F1", "#D9E9F2", "#DEEBF2", "#E2EDF3", "#E6EFF4", "#EAF1F5", "#EFF3F5", "#F3F5F6", "#F7F7F7", "#F7F7F7", "#F8F4F2", "#F8F1EC", "#F9EEE7", "#FAEBE2", "#FAE7DC", "#FBE4D7", "#FCE1D2", "#FCDECC", "#FDDBC7", "#FDDBC7", "#FCD5BF", "#FBCFB8", "#FAC9B0", "#F9C3A8", "#F8BDA1", "#F7B799", "#F6B191", "#F5AB8A", "#F4A582", "#F4A582", "#F19D7C", "#ED9676", "#EA8E70", "#E7866A", "#E37F65", "#E0775F", "#DD6F59", "#D96853", "#D6604D", "#D6604D", "#D25849", "#CE5045", "#CA4842", "#C6403E", "#C2383A", "#BE3036", "#BA2833", "#B6202F", "#B2182B", "#B2182B", "#AA152A", "#A11328", "#991027", "#910D26", "#880B24", "#800823", "#780522", "#6F0320", "#67001F")
  }

  # Reverse colours; i.e. so that the heatmap shows the phenotypic changes from morph y (target) to morph x (ref)
  # contrary to the default; heatmap shows the changes from ref to target
  # this is useful if you want to show morph y AND changes
  if (direction == "reverse") {
    coloursofthewind <- rev(coloursofthewind)}

  return(meshDist(x, y, steps = 100, rampcolors = coloursofthewind, from = -dists, to = dists))
}
