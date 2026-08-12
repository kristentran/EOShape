#' Set DAMorphometrics MM data path in ~/.Renviron
#'
#' Writes or updates DAMORPHOMETRICS_MM_DATA_PATH in the user's ~/.Renviron
#' so it is available in all future R sessions.
#'
#' @param path Server path
#' @param storage hs or mm data path.
#' @param renviron_path Optional path to Renviron file (default: "~/.Renviron").
#' @param create Logical; create file if missing (default TRUE).
#' @return Invisibly returns the path written.
#' @export
#'

set_servers = function(path, storage, renviron_path = "~/.Renviron", create = TRUE) {
  stopifnot(is.character(path), length(path) == 1, nzchar(path))
  storage <- tolower(storage)

  if(storage == "mm"){
    key = "CRANIOMORPH_MM_DATA_PATH"
    suffix = "//Hallgrimsson//Mutant_Data//Adult"
  } else if(storage == "hs"){
    key = "CRANIOMORPH_HS_DATA_PATH"
    suffix = "//FaceBase_3//Data//Basic//Complete"
  } else {
    stop("storage must be either 'mm' or 'hs'")
  }

  # Allow either a server root or a fully specified data path.
  path <- trimws(path)
  if (!endsWith(path, suffix)) {
    path <- paste0(path, suffix)
  }

  # Normalize path for consistency (keep as entered if normalization fails)
  norm_path <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE),
                        error = function(e) path)

  renv <- path.expand(renviron_path)
  newline <- paste0(key, "=", norm_path)

  if (!file.exists(renv)) {
    if (!isTRUE(create)) {
      stop("Renviron file does not exist: ", renv)
    }
    dir.create(dirname(renv), recursive = TRUE, showWarnings = FALSE)
    file.create(renv)
  }

  lines <- readLines(renv, warn = FALSE)
  key_pattern <- paste0("^\\s*", key, "=")

  if (any(grepl(key_pattern, lines))) {
    lines[grepl(key_pattern, lines)] <- newline
  } else {
    # Keep file tidy: add a blank line before appending if needed
    if (length(lines) > 0 && tail(lines, 1) != "") {
      lines <- c(lines, "")
    }
    lines <- c(lines, newline)
  }

  writeLines(lines, renv, useBytes = TRUE)

  message("Updated ", renv, " with ", key, ". Restart R to apply.")
  invisible(norm_path)
}
