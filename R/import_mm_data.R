#' Import .tag files to 3D array
#'
#' Internal function to import mouse data from BH Lab servers. This function is minimally adapted from morpho.tools.GM::tag2array
#'
#' @param dir directory to find .tag files
#' @param string_del delete parts of .tag file names
#' @param names vector of samples to read. If provided, only read in the provided names
#'
#' @return landmarks as p x k x n

import_mm_data <- function (dir = NULL, ID = NULL, string_del = NULL, propagated = FALSE,
                            save.txt = FALSE, names = NULL){
  if (is.null(dir) == TRUE) {
    path <- getwd()
  }
  else {
    path <- dir
  }
  tag_list <- dir(path = path, pattern = "*.tag", full.names = TRUE)

  # filter taglist
  if(!is.null(names)){
    tag_list <- tag_list[basename(tag_list) %in% paste0(names, "_Cranium_Landmarks.tag")]
  }

  n_land <- vector("numeric", length = length(tag_list))
  for (i in 1:length(tag_list)) {
    if (isTRUE(propagated) == TRUE) {
      n_land[i] <- length(count.fields(tag_list[[i]])) -
        4
    }
    else {
      n_land[i] <- length(count.fields(tag_list[[i]])) -
        3
    }
  }
  if (is.null(ID) == TRUE) {
    if (is.null(string_del) == TRUE) {
      dimnames_tag <- gsub(".tag", "", tag_list)
    }
    else {
      # debugged path argument to function without set working directory
      dimnames_tag <- gsub(path, "", gsub(string_del, "", gsub(".tag","", tag_list)))
      dimnames_tag <- stringr::str_remove_all(dimnames_tag, pattern = "/")
    }
  }
  else {
    dimnames_tag <- ID
  }
  if (length(unique(n_land)) != 1) {
    stop("Specimens have different number of landmarks.")
  }
  LM_array <- array(data = NA, dim = c(n_land[1], 3, length(tag_list)))
  dimnames(LM_array)[[3]] <- dimnames_tag
  LM_type <- vector(mode = "character", length(tag_list))
  if (isTRUE(propagated) == TRUE) {
    for (i in 1:length(tag_list)) {
      # Add progress bar
      LM_array[, 1, i] <- suppressWarnings(read.table(file = tag_list[[i]],
                                                      skip = 5, sep = " ", header = F, fill = T))[, 2]
      LM_array[, 2, i] <- suppressWarnings(read.table(file = tag_list[[i]],
                                                      skip = 5, sep = " ", header = F, fill = T))[, 3]
      LM_array[, 3, i] <- suppressWarnings(read.table(file = tag_list[[i]],
                                                      skip = 5, sep = " ", header = F, fill = T))[, 4]
      cat("\r", i, "of", length(tag_list))
    }
  }
  else { # if not propagated
    for (i in 1:length(tag_list)) {
      LM_array[, 1, i] <- suppressWarnings(read.table(file = tag_list[[i]],
                                                      skip = 4, sep = " ", header = F, fill = T))[, 2]
      LM_array[, 2, i] <- suppressWarnings(read.table(file = tag_list[[i]],
                                                      skip = 4, sep = " ", header = F, fill = T))[, 3]
      LM_array[, 3, i] <- suppressWarnings(read.table(file = tag_list[[i]],
                                                      skip = 4, sep = " ", header = F, fill = T))[, 4]
      cat("\rRead", i, "of", length(tag_list))
    }
  }
  if (isTRUE(save.txt) == TRUE) {
    for (i in 1:dim(LM_array)[3]) {
      write.table(LM_array[, , i], paste0(dimnames(LM_array)[[3]][i],
                                          ".txt"), col.names = FALSE, row.names = FALSE)
    }
  }

  # Return warning if not numeric
  if(!is.numeric(LM_array)){warning("landmark coordinates are not numeric")}
  return(LM_array)
}
