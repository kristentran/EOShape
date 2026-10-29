#' Query data from BHLAB servers
#'
#' Choose data for human datasets or mouse datasets. Servers need to be mounted first.
#'
#' @param species select from mouse (mm) or human (hs) database.
#' @param dataset from selected database, choose specific dataset to query or "all"; if not provided, then function will prompt user to choose from all available datasets.
#' @param lmk_set for human data, choose dense or sparse landmark sets.
#' @param symmetrized for human data, choose the original (FALSE) or symmetric (TRUE) dataset.
#' @param drop_missing if TRUE, drop individuals/samples without matching landmark and metadata.
#'
#' @return lmks landmarks as p x k x n array.
#' @return meta complete metadata table, as available in BHLAB's servers.
#'
#' @import stringr
#' @import crayon
#'
#' @export

query_data = function(
    species = c("mm", "hs"),
    dataset = NULL,
    lmk_set = "dense",
    symmetrized = TRUE,
    drop_missing = TRUE
    # TODO: add age filtering - important for when embryo landmark data is compatible
    # stage = c("all", "adults", "juvenile"),
    # TODO: Probably dbconnect to SQL when Fernando is done
    # source = NULL,
){

  resolve_source_dir <- function(species) {
    env_key <- if (species == "mm") {
      "EOS_MM_DATA_PATH"
    } else {
      "EOS_HS_DATA_PATH"
    }

    env_path <- trimws(Sys.getenv(env_key, unset = ""))
    if (nzchar(env_path)) {
      return(env_path)
    }

    default_file <- system.file(
      "extdata",
      paste0("pathto", toupper(species), "data.txt"),
      package = "EOS"
    )

    if (!nzchar(default_file)) {
      stop(
        "Could not find default path file for species '", species, "'. ",
        "Please reinstall the package or set the path with set_servers()."
      )
    }

    default_path <- trimws(readLines(default_file, warn = FALSE, n = 1))
    if (!nzchar(default_path)) {
      stop(
        "Default path file is empty for species '", species, "'. ",
        "Please set the path with set_servers()."
      )
    }

    default_path
  }

  # Check if arguments are compatible
  # species
  if(length(species) != 1 | !any(tolower(species) %in% c("mm", "hs"))){
    stop("species must be mm (Mus musculus), or hs (homo sapiens)")
  } else {
    species = tolower(species)
    source_dir = resolve_source_dir(species)
  }

  # set specific objects for hs vs. mm datasets
  # keep it like this and not ifelse to prepare for more differences
  if(species == "hs"){
    if(length(lmk_set) != 1 | !lmk_set %in% c("dense", "sparse"))
      stop("lmk_set must be dense (5k landmarks) or sparse (63 landmarks)")
    dataset_col = "Syndrome"
    ID_col = "Image_Name"
  } else if(species == "mm"){
    dataset_col = "Dataset"
    ID_col = "Biosample"
  }

  # Read meta
  meta = read.csv(paste0(source_dir, ifelse(species == "hs", "//complete_cleaned_sample_cov.csv", "//Adult_Metadata//MusMorph_Adult_Metadata.csv")))
  datasets = unique(meta[,dataset_col])

  # reorder datasets to alphabetical order
  datasets = str_sort(datasets)

  # Confirm if dataset exists in meta
  if(!is.null(dataset)){
    if(length(dataset) != 1)
      stop("Provide only one dataset for ", species)
    if(!any(dataset %in% c(datasets, "all"))){
      cat(red("Dataset (", dataset, ") was not found, please select another or check spelling and case\n")); dataset = NULL}
  }

  # If no dataset is provided, then allow user to select dataset
  while(is.null(dataset)){
    cat(inverse("Please choose from the following datasets to query or type 'all':\n"))
    print(datasets)
    dataset = readline(prompt = "Enter selection here (case sensitive) or \"all\" to get all datasets: \n")
    if(!any(dataset %in% c(datasets, "all"))){cat(red(dataset, "was not found, please select another or check spelling\n")); dataset = NULL}
  }

  # subset meta based on dataset if all was not selected
  if (dataset != "all"){
    meta = meta[meta[,dataset_col] == dataset,]
  }

  # read landmarks
  cat(blue("Reading landmark data...\n"))

  if(species == "mm"){
    # Reads only mice needed
    lmks = import_mm_data(dir = paste0(source_dir, "//Landmarks//Cranium//"), string_del = "_Cranium_Landmarks", propagated = TRUE,
                          names = meta[,ID_col] # specify list to read
    )
    cat("\n")

    # In case .tag were not read in properly
    if(is.character(lmks)){
      if(grepl("[A-Za-z]", lmks[1,1,1]))
        lmks <- lmks[-1,,]
      lmks[,,] <- lmks(as.numeric(lmks), dim = dim(lmks), dimnames = dimnames(lmks))
    }

    # Reduce number of landmarks if using automated 844 set (specific to mm)
    if(dim(lmks)[1] == 844){
      lmks=lmks[-c(337:430,  94, 325),,]
      rownames(lmks) = 1:748
    }
  } else if(species == "hs"){
    # Read entire dataset (since it is a prepared .rds file)
    lmks_path = list.files(source_dir, "*.rds")
    lmks_path = lmks_path[grepl(ifelse(lmk_set == "dense", "5k", "63"), lmks_path)]
    lmks_path = ifelse(symmetrized, lmks_path[grepl("symmetric", lmks_path)], lmks_path[grepl("original", lmks_path)])
    if(is.na(lmks_path)) # stop if no dataset found
      stop("The dataset with the queried parameters does not exist: ", lmk_set, " and ", ifelse(symmetrized, "symmetrized", "nonsymmetrized"))
    if(length(lmks_path) > 1){
      cat("Multiple datasets were found. Please select from the following:\n")
      indexed_path = paste0("[",seq(1:length(lmks_path)),"]: ", lmks_path,"\n")
      cat(indexed_path) # if multiple found, allow user to use
      ans = readline(prompt = "Enter the index specifying the dataset here: ")
      lmks_path[ans]
    }
    # read selected file
    lmks = readRDS(paste0(source_dir, "//", lmks_path))
    # subset if dataset was specified
    if(dataset != "all"){
      lmks = lmks[,,meta[,ID_col]]
    }
  }

  # Drop animals that do not have both lmks and meta
  if(drop_missing){
    cat(blue("Dropping", abs(nrow(meta) - dim(lmks)[[3]]), "samples/participants with missing landmark or metadata.\n"))
    meta = meta[meta[,ID_col] %in% dimnames(lmks)[[3]],]
    lmks = lmks[,,meta[,ID_col]]
  }

  rownames(meta) = meta[,ID_col]

  cat(blue("Dataset query complete! :D\n"))
  return(list("lmks" = lmks, "meta" = meta))
}
