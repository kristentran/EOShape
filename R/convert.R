#' Convert landmarks between useful formats
#'
#' Convert three-dimensional landmarks between 1D (row), 2D (matrix), and 3D (array). This function automatically converts to p x k (landmarks as rows x dimensions as columns) format. If convert_always is TRUE, then p x k formats will be converted to p*k format (e.g. 1x, 1y, 1z,...).
#'
#' @param lmks vector, matrix (2D), or 3D array of landmarks to convert.
#' @param convert_always if TRUE, convert to p*k format.
#' @param silence if FALSE, show how the landmarks are detected and how they will be converted.
#'
#' @return landmarks converted to new format.
#'
#' @importFrom geomorph two.d.array
#'
#' @export

convert = function(lmks, convert_always = FALSE, silence = TRUE){
  if(is.vector(lmks) || (is.matrix(lmks) & nrow(lmks) == 1)){
    if(!silence){cat(blue("Landmarks are 1D, converting vector to p x k matrix.\n"))}
    p = length(lmks)/3
    n = 1
    lmks = matrix(lmks, byrow = TRUE, ncol = 3)
  } else if(length(dim(lmks)) == 2){
    if(ncol(lmks) == 3){
      if(convert_always){
        if(!silence){cat(blue("Landmarks detected as p x k matrix for 1 specimen... converting to vector.\n"))}
        lmks = as.vector(t(lmks))
      } else {
        if(!silence){cat(blue("Landmarks detected as p x k matrix for 1 specimen. Nothing to convert.\n"))}}
    } else {
      if(!silence){cat(blue("Landmarks detected as n x (p*k) matrix... converting to 3D array\n"))}
      p = ncol(lmks)/3
      n = nrow(lmks)
      lmks3d = array(NA, dim = c(p,3,n))
      for (i in 1:n){
        lmks3d[,,i] = matrix(lmks[i,], byrow = TRUE, ncol = 3)
      };
      if(!is.null(rownames(lmks))){
        n_names = rownames(lmks)
        dimnames(lmks3d)[[3]] = n_names}
      lmks = lmks3d
    }
  } else if(length(dim(lmks)) == 3){
    if(convert_always){
      if(!silence){cat(blue("Landmarks detected as p x k x n array. Converting to n x (p*k) matrix.\n"))}
      lmks = geomorph::two.d.array(lmks)
    } else {
      if(!silence){cat(blue("Landmarks detected as p x k x n array. Nothing to convert.\n"))}
    }
  } else {
    stop("Landmark format not recognized")
  }
  return(lmks)
}

