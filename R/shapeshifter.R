#' Deform atlas mesh to target landmarks.
#'
#' If the atlas' vertices are equal to the number of landmarks, the vertices are simply replaced. If not, then this serves as a wrapper function for Morpho::tps3D.
#'
#' @param target_lmks vector, matrix (2D), or 3D array of target landmarks to which the atlas is deformed.
#' @param atlas_mesh mesh3d object of your atlas to deform.
#' @param atlas_lmks vector, matrix (2D), or 3D array of landmarks associated with atlas_mesh. Only necessary if the atlas vertices are not equal to the target landmarks.
#'
#' @return mesh3D object of deformed atlas.
#'
#' @import Morpho
#'
#' @export

shapeshifter = function(target_lmks, atlas_mesh, atlas_lmks = NULL){
  stopifnot(any("mesh3d" %in% class(atlas_mesh)), is.matrix(target_lmks)|is.array(target_lmks)|is.numeric(target_lmks))
  if(!is.null(atlas_lmks)) stopifnot(is.matrix(atlas_lmks)|is.array(atlas_lmks)|is.numeric(atlas_lmks))

  target_lmks = convert(target_lmks)

  if(length(dim(target_lmks)) == 2){
    target_lmks = array(target_lmks, dim = c(nrow(target_lmks), ncol(target_lmks), 1))
    dimnames(target_lmks)[[3]] = "onlyone"}

  shifter_internal = function(target_lmks, atlas_mesh, atlas_lmks = NULL){
    if(is.null(atlas_lmks)){
      stopifnot(dim(target_lmks) == dim(t(atlas_mesh$vb[1:3,])))
      target_mesh = atlas_mesh
      target_mesh$vb[1:3,] = t(target_lmks)
    } else {
      target_mesh = tps3d(atlas_mesh, atlas_lmks, target_lmks)
    }
    return(target_mesh)
  }

  results = list()
  for(i in 1:dim(target_lmks)[3]){
    results[[dimnames(target_lmks)[[3]][i]]] = shifter_internal(target_lmks = target_lmks[,,i], atlas_mesh = atlas_mesh, atlas_lmks = atlas_lmks)
  }
  if(length(results) == 1){results = results[["onlyone"]]}
  return(results)
}

