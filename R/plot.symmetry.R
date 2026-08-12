#' @name plot.symmetry
#' @title plot.symmetry
#' @alias plot.symmetry
#'
#' @description
#'   This function plots paired landmarks in a dataset with object symmetry, similarly to MorphoJ (Klingenberg, 2011). This is a copy of plot.symmetry from MVG's symmetry package.
#'
#' @usage
#'   plot.symmetry(A, sym.pairs, sym.plane = NULL)
#' @param item{A} An object of the class "array" that contains Procrustes shape variables for 'n' specimens
#'   }
#'   \item{sym.pairs}{
#'     A matrix defining paired landmarks
#'   }
#'   \item{sym.pairs}{
#'     An optional parameter defining whether to plot the symmetry plane
#'   }
#'
#' @value
#'   This function returns a plot with the paired landmarks
#'
#'
#' @examples
#'     left <- c(58, 15, 19, 21, 23, 25, 27, 31, 33, 35, 37, 43, 39, 45, 47, 49, 41, 29, 17, 6, 4, 8, 10, 12, 63, 62, 61, 60, 59, 51)
#'     right <- c(57, 14, 18, 20, 22, 24, 26, 30, 32, 34, 36, 42, 38, 44, 46, 48, 40, 28, 16, 7, 5, 9, 11, 13, 56, 55, 54, 53, 52, 50)
#'     sym.pairs <- cbind(left, right)
#'     plot.sym(A, sym.pairs)
#'     plot.sym(A, sym.pairs, sym.plane = "yz")
#'     plot.sym(A, sym.pairs, sym.plane = "xy")
#'
#' @author Marta Vidal-Garcia, \email{marta.vidalga@@gmail.com}
#' @references Klingenberg, C.P. (2011) MorphoJ: an integrated software package for geometric morphometrics. Molecular ecology resources. 11: 353-357
#' @export


plot.symmetry <- function(A, sym.pairs, sym.plane = NULL){
  require(rgl)
  specimen_m <- as.data.frame(geomorph:::mshape(geomorph::gpagen(A, print.progress = FALSE)$coords))
  side.1 <- sym.pairs[,1]
  side.2 <- sym.pairs[,2]
  pairs <- c(side.1, side.2)
  all <- as.numeric(row.names(specimen_m))
  diff <- setdiff(all, pairs)
  par3d(windowRect = c(20, 30, 800, 800))
  plot3d(x = specimen_m[pairs,1], y = specimen_m[pairs,2], z = specimen_m[pairs,3], col = "blue", type="s", aspect = 'iso', size=0.6,
         xlab = "x", ylab = "y", zlab = "z")
  plot3d(x = specimen_m[diff,1], y = specimen_m[diff,2], z = specimen_m[diff,3], col = "black", type="s", aspect = 'iso', size=0.6, add = TRUE,
         xlab = "x", ylab = "y", zlab = "z")
  text3d(specimen_m[,1], specimen_m[,2], specimen_m[,3], texts = paste(1:(dim(specimen_m)[1])), adj=1.5)
  if (!is.null(sym.plane)) {
    if (sym.plane == "xz"){
      planes3d(a=0,b=1,c=0,d=0, alpha=0.2, col="gray") # xz-plane
    }
    if (sym.plane == "yz"){
      planes3d(a=1,b=0,c=0,d=0, alpha=0.2, col="gray") # yz-plane
    }
    if (sym.plane == "xy"){
      planes3d(a=0,b=0,c=1,d=0, alpha=0.2, col="gray") # xy-plane
    }
  }
  else{
    sym.plane <- NULL
  }
  pairs_plot <- data.frame(side.1, side.2)
  for (i in 1:length(side.1)){
    segments3d(specimen_m[c(side.1[i],side.2[i]), ], col="red", lwd = 2)
  }
}
