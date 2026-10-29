symmetrise_wrapper <- function(lmks){
  rem <- c(337:430,  94, 325) # These are duplicated or extra landmarks.

  # First, check how many landmarks
  if(dim(lmks)[1] == 844){
    lmks <- lmks[-rem,,]
    rownames(lmks) <- 1:748
  } else if(dim(lmks)[1] != 748){
    stop("Check the number of landmarks")
  }

  # Symmetrise landmarks...
  sym <- readRDS("inst/extdata/MM_landmark_sym.Rds")
  require(Morpho)
  sym_lmks <- Morpho::symmetrize(lmks, as.matrix(sym$sym.pairings))
  rm(rem); return(sym_lmks)
}
