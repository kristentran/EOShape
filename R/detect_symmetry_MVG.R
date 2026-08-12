#' @name detect.symmetry
#' @title detect.symmetry
#' @alias detect.symmetry
#'
#' @description
#'   This function detects paired landmarks in a dataset with object symmetry, similarly to MorphoJ (Klingenberg, 2011). This is a copy of detect.symmetry from MVG's symmetry package.
#'
#' @usage
#'   detect.symmetry(A, sym.plane = "xz", non.sym = NULL, plot = TRUE)
#' @param item{A} An object of the class "array" that contains Procrustes shape variables for 'n' specimens
#'   }
#'   \item{non.sym}{
#'     An optional variable defining non-paired landmarks
#'   }
#'
#' @value
#'   This function returns a list with three objects. The first object is the paired landmarks on side 1, the second object is the paired landmarks on side 2, and the third object is a matrix with paired landmarks from both sides. By default, the estimated paired landmarks and non-paired landmarks are plotted.
#'
#' @details
#'   This function returns the paired landmarks.
#'
#'
#' @examples
#'     detect.symmetry(A, sym.plane = "xz", non.sym = NULL, plot = TRUE)
#'
#' @author Marta Vidal-Garcia, \email{marta.vidalga@@gmail.com}
#' @references Klingenberg, C.P. (2011) MorphoJ: an integrated software package for geometric morphometrics. Molecular ecology resources. 11: 353-357
#' @export


detect.symmetry <- function(A, sym.plane = "xz", non.sym = NULL, plot = TRUE){
  if (sym.plane == "xz"){
    specimen_m <- as.data.frame(geomorph:::mshape(geomorph::gpagen(A, print.progress = FALSE)$coords)) # estimate the mean shape
    specimen_refl <- data.frame(specimen_m[,1], specimen_m[,2][sapply(specimen_m[,2], is.numeric)] * -1, specimen_m[,3])
    side.1 <- which(specimen_m[,2] < 0)
    side.2 <- as.integer(setdiff(row.names(specimen_m), side.1))
    non <- abs(length(side.2) - length(side.1)) # estimated number of non-symetrical points
    h <- abs(specimen_m[,2])
    names(h) <- row.names(specimen_m)
    non_sym <- names(sort(h)[1:non])
    paired_sym <- matrix(data = NA, ncol = 2, nrow = dim(specimen_m)[1])
    for (i in side.1) {
      distance <- numeric(length = length(side.2))
      for (j in side.2){
        distance[j] <- sqrt((specimen_m[j,1] - specimen_refl[i,1])^2 + (specimen_m[j,2] - specimen_refl[i,2])^2 + (specimen_m[j,3] - specimen_refl[i,3])^2)
      }
      distance[distance == 0] <- 999
      paired_sym[i,1] <- i
      paired_sym[i,2] <- which.min(distance)
    }
    paired_sym <- as.data.frame(na.omit(paired_sym))
    colnames(paired_sym) <- c("side.1", "side.2")
  }
  if (sym.plane == "xy"){
    specimen_m <- as.data.frame(mshape(gpagen(A)$coords)) # estimate the mean shape
    specimen_refl <- data.frame(specimen_m[,1], specimen_m[,2], specimen_m[,2][sapply(specimen_m[,3], is.numeric)] * -1)
    side.1 <- which(specimen_m[,3] < 0)
    side.2 <- as.integer(setdiff(row.names(specimen_m), side.1))
    non <- abs(length(side.2) - length(side.1)) # estimated number of non-symetrical points
    h <- abs(specimen_m[,2])
    names(h) <- row.names(specimen_m)
    non_sym <- names(sort(h)[1:non])
    paired_sym <- matrix(data = NA, ncol = 2, nrow = dim(specimen_m)[1])
    for (i in side.1) {
      distance <- numeric(length = length(side.2))
      for (j in side.2){
        distance[j] <- sqrt((specimen_m[j,1] - specimen_refl[i,1])^2 + (specimen_m[j,2] - specimen_refl[i,2])^2 + (specimen_m[j,3] - specimen_refl[i,3])^2)
      }
      distance[distance == 0] <- 999
      paired_sym[i,1] <- i
      paired_sym[i,2] <- which.min(distance)
    }
    paired_sym <- as.data.frame(na.omit(paired_sym))
    colnames(paired_sym) <- c("side.1", "side.2")
  }
  if (sym.plane == "yz"){
    specimen_m <- as.data.frame(mshape(gpagen(A)$coords)) # estimate the mean shape
    specimen_refl <- data.frame(specimen_m[,2][sapply(specimen_m[,1], is.numeric)] * -1, specimen_m[,2], specimen_m[,3])
    side.1 <- which(specimen_m[,1] < 0)
    side.2 <- as.integer(setdiff(row.names(specimen_m), side.1))
    non <- abs(length(side.2) - length(side.1)) # estimated number of non-symetrical points
    h <- abs(specimen_m[,2])
    names(h) <- row.names(specimen_m)
    non_sym <- names(sort(h)[1:non])
    paired_sym <- matrix(data = NA, ncol = 2, nrow = dim(specimen_m)[1])
    for (i in side.1) {
      distance <- numeric(length = length(side.2))
      for (j in side.2){
        distance[j] <- sqrt((specimen_m[j,1] - specimen_refl[i,1])^2 + (specimen_m[j,2] - specimen_refl[i,2])^2 + (specimen_m[j,3] - specimen_refl[i,3])^2)
      }
      distance[distance == 0] <- 999
      paired_sym[i,1] <- i
      paired_sym[i,2] <- which.min(distance)
    }
    paired_sym <- as.data.frame(na.omit(paired_sym))
    colnames(paired_sym) <- c("side.1", "side.2")
  }
  dups <- as.data.frame(table(paired_sym$side.2))
  n.dups <- as.numeric(as.character(dups[dups$Freq > 1,]$Var1))
  check <- vector("list", length(n.dups))
  best.pair <- vector("list", length(n.dups))
  for (i in 1:length(n.dups)){
    check[[i]] <- c(paired_sym[,1][which(paired_sym[,2] == n.dups[i])])
    # which landmark pair best with the duplicated ones?
    dist.sec <- vector("numeric", length(check[[i]]))
    for (j in check[[i]]){
      dist.sec[j] <- sqrt((specimen_m[j,1] - specimen_refl[n.dups[i],1])^2 + (specimen_m[j,2] - specimen_refl[n.dups[i],2])^2 + (specimen_m[j,3] - specimen_refl[n.dups[i],3])^2)
    }
    best.pair[[i]] <- check[[i]][which.min(dist.sec[c(check[[i]])])]
    paired_sym$side.2[paired_sym$side.2 == n.dups[i]] <- "NA"
    paired_sym$side.2[which(paired_sym$side.1 == best.pair[[i]])] <- n.dups[i]
  }
  non_sym <- unique(c(non_sym, paired_sym$side.1[which(paired_sym$side.2 == "NA")]))
  pairedLM <- na.omit(suppressWarnings(transform(paired_sym, side.2 = as.numeric(side.2))))
  side.1 <-pairedLM$side.1
  side.2 <-pairedLM$side.2
  non.sym <- setdiff(row(specimen_m)[,1], c(side.1, side.2))
  if(plot)
    plot.symmetry(A, sym.pairs = pairedLM, sym.plane)
  return(list("side.1" = side.1, "side.2" = side.2, "non.sym" = non.sym, "pairedLM" = as.matrix(pairedLM)))
}
