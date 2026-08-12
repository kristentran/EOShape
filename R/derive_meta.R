#' Set up metadata for analysis
#'
#' Automatically detects identifiers, sex, and age. Identifiers are detected as
#' the first column containing only unique values. Sex is detected as the first
#' column containing only F/M (case insensitive) or binary 0/1 (encoded as
#' M = 1, F = 0). Age is detected as the first non-binary numeric-like variable
#' that is not already used as ID or sex. Non-digits is removed from Age.
#'
#' @param borealis borealis object initialised with \code{initialise_borealis()}
#' @param column_definitions optional named list with elements ID, Sex, Age,
#'   Confounders (or Confounder), and Dependants (or Dependant) — column indices
#'   into \code{borealis$meta}. When supplied, these are used as-is with no
#'   detection or automated overwrite.
#' @param automated if TRUE and \code{column_definitions} is NULL, detect ID/Sex/Age,
#'  take the first remaining column as the dependant, and assume no confounders.
#'  If FALSE and \code{column_definitions} is NULL, prompt interactively.
#'  Ignored for column selection when \code{column_definitions} is supplied
#'  (still skips interactive confirmation when TRUE).
#' @param impute if TRUE, impute missing data using \code{mice::mice}. Extra
#'  arguments will be passed to \code{mice::mice()}. if FALSE, missing data will
#'  be removed
#'
#' @return borealis object for subsequent analyses
#'
#' @export
#'

derive_meta = function(borealis, column_definitions = NULL, automated = FALSE, impute = FALSE, ...){
  require(crayon)

  if(!inherits(borealis, "borealis"))
    stop("Object provided is not a borealis object. Please create using initialise_borealis().")
  if(!"meta" %in% names(borealis) || !is.data.frame(borealis$meta))
    stop("Please check meta in borealis object.")
  # check if automated and impute are single value logicals
  if(!is.logical(automated) || !is.logical(impute) || length(automated) > 1 || length(impute) > 1)
    stop("automated and impute parameters must be single value logicals")

  extracted = borealis$meta
  n_col = ncol(extracted)

  is_numeric_like = function(x){
    vals = na.omit(x)
    length(vals) > 0 && all(!is.na(suppressWarnings(as.numeric(as.character(vals)))))
  }

  first_matching_col = function(predicate, exclude = numeric(0)){
    for(i in setdiff(seq_len(n_col), exclude)){
      if(predicate(extracted[[i]]))
        return(i)
    }
    return(numeric(0))
  }

  prompt_cols = function(label, allow_multiple = FALSE){
    if(automated)
      return(numeric(0))

    if(label %in% c("Sex", "Age")){
      cat(bold(red(label, "not detected, please identify the column by number:\n")))
    } else {
      cat(blue("Please identify columns containing", bold(label), "variables by number:\n"))
    }
    cat(paste0("[", seq_len(n_col), "] ", colnames(extracted), "\n"), sep = "")

    if(allow_multiple){
      ans = readline(prompt = "Enter here (if multiple, separate by comma or space. Hit enter if none):")
      ans = as.numeric(unlist(strsplit(ans, "[^[:alnum:]]+")))
    } else {
      ans = as.numeric(readline(prompt = "Enter here (hit enter if none):"))
    }
    if(all(is.na(ans))) numeric(0) else ans
  }

  accepted = FALSE
  while(!accepted){
    # If column_definitions is supplied, use it as-is — no detection or
    # automated overwrite of Dependant / Confounder.
    if(!is.null(column_definitions)){
      ID_col = column_definitions[["ID"]]
      sex_col = column_definitions[["Sex"]]
      age_col = column_definitions[["Age"]]
      confounder_col = column_definitions[["Confounders"]]
      if(is.null(confounder_col))
        confounder_col = column_definitions[["Confounder"]]
      dependant_col = column_definitions[["Dependants"]]
      if(is.null(dependant_col))
        dependant_col = column_definitions[["Dependant"]]

      if(is.null(ID_col)) ID_col = numeric(0)
      if(is.null(sex_col)) sex_col = numeric(0)
      if(is.null(age_col)) age_col = numeric(0)
      if(is.null(confounder_col)) confounder_col = numeric(0)
      if(is.null(dependant_col)) dependant_col = numeric(0)
    } else {
      ID_col = first_matching_col(function(x) !any(duplicated(x)))
      if(!length(ID_col))
        stop("No columns in metadata with unique values. Please include a column containing unique identifiers for each specimen in meta")

      sex_col = first_matching_col(function(x){
        vals = na.omit(x)
        if(!length(vals)) return(FALSE)
        if(all(toupper(as.character(vals)) %in% c("M", "F"))) return(TRUE)
        is_numeric_like(vals) && all(as.numeric(as.character(vals)) %in% c(0, 1))
      }, exclude = ID_col)
      if(!length(sex_col))
        sex_col = prompt_cols("Sex")

      age_col = first_matching_col(function(x){
        vals = na.omit(x)
        is_numeric_like(vals) && !all(as.numeric(as.character(vals)) %in% c(0, 1))
      }, exclude = c(ID_col, sex_col))
      if(!length(age_col))
        age_col = prompt_cols("Age")

      if(automated){
        remaining = setdiff(seq_len(n_col), c(ID_col, sex_col, age_col))
        dependant_col = if(length(remaining)) remaining[[1]] else numeric(0)
        confounder_col = numeric(0)
      } else {
        confounder_col = prompt_cols("confounder", allow_multiple = TRUE)
        dependant_col = prompt_cols("dependant", allow_multiple = TRUE)
      }
    }

    if(!length(ID_col))
      stop("column_definitions must include ID (column index of unique specimen identifiers).")

    derived = data.frame(identifiers = extracted[,ID_col], stringsAsFactors = FALSE)
    for(colu in c("sex_col", "age_col", "confounder_col", "dependant_col")){
      idx = get0(colu)
      if(!length(idx)) next
      df = extracted[, idx, drop = FALSE]
      names(df) = if(any(colu %in% c("confounder_col", "dependant_col"))) paste0(gsub("_col", "", colu), "_", names(df)) else gsub("_col", "", colu)
      for(j in seq_along(df)){
        if(is_numeric_like(df[[j]]))
          df[[j]] = as.numeric(df[[j]])
      }
      derived = cbind(derived, df)
    }

    if(automated){
      accepted = TRUE
    } else {
      cat(bold(blue("Please confirm metadata structure:\n")))
      print(head(derived))
      ans = readline(prompt = "Hit enter to accept or type anything to redo:")
      if(ans == "") accepted = TRUE
    }
  }


  # If metadata is accepted, format:
  if(length(sex_col) && "sex" %in% names(derived)){
    if(is.character(derived$sex) || is.factor(derived$sex))
      derived$sex = as.integer(toupper(as.character(derived$sex)) == "M")
    else
      derived$sex = as.integer(as.numeric(derived$sex) == 1)
  }

  if(length(age_col) && "age" %in% names(derived))
    derived$age = as.numeric(gsub("\\D+", "", derived$age))

  if(is.null(dimnames(borealis$lmks)[[3]]))
    dimnames(borealis$lmks)[[3]] = derived$identifiers

  # if impute is TRUE, then impute missing values. If not, remove NA rows
  if(impute && !automated){
    require(mice)
    cat(crayon::blue("Examine missing data pattern:\n"))
    print(mice::md.pattern(derived))
    ans = readline(prompt = "Hit enter to proceed with data imputation or type anything to remove missing values instead")
    if(ans != "") impute = FALSE
  }

  if(impute){
    derived = mice::mice(derived, ...)
    borealis$analysis_params$imputed_meta = derived$where
    derived = complete(derived, 1)
    borealis$analysis_params$imputed_meta = TRUE
  } else {
    derived = na.omit(derived)
  }

  identifiers = derived$identifiers
  identifiers = identifiers[identifiers %in% dimnames(borealis$lmks)[[3]]]

  derived = derived[match(identifiers, derived$identifiers),]
  borealis$lmks = borealis$lmks[,,identifiers]

  if(any(!is.na(suppressWarnings(as.numeric(identifiers)))))
    identifiers = paste0("ID_", identifiers)

  rownames(derived) = identifiers
  derived$identifiers = NULL
  dimnames(borealis$lmks)[[3]] = identifiers

  borealis$meta = derived
  borealis$analysis_log = c(borealis$analysis_log, "metadata automatically filtered for singular columns")
  borealis$identifiers = identifiers

  return(borealis)
}
