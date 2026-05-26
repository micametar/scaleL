#' Validate and normalize an input data.frame for [scaleL()]
#'
#' Checks required columns, coerces types, fills in optional columns with
#' defaults, and verifies that all rows have `n >= 2` and `sd > 0`.
#'
#' @param d Input data.frame.
#' @param verbose Logical; if `TRUE`, prints a note when `true_L` is detected.
#' @return Validated data.frame, ready for the imputation pipeline.
#' @keywords internal
validate_input <- function(d, verbose = FALSE) {
  required_cols <- c("study_id", "mean", "sd", "n")
  missing_cols <- setdiff(required_cols, names(d))
  if (length(missing_cols) > 0) {
    stop("Input is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!("L" %in% names(d))) d$L <- NA_integer_
  if (!("scale_min" %in% names(d))) d$scale_min <- 1
  if (!("measure_id" %in% names(d))) {
    d$measure_id <- paste0("m", seq_len(nrow(d)))
  }
  if (!("instrument" %in% names(d))) d$instrument <- NA_character_

  if ("true_L" %in% names(d) && isTRUE(verbose)) {
    message("Note: input has a `true_L` column. It is NOT used in ",
            "imputation; it is carried through for validation only.")
  }

  for (col in c("mean", "sd", "n")) {
    if (!is.numeric(d[[col]])) {
      d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
    }
    if (any(is.na(d[[col]]))) {
      stop("Column '", col, "' contains non-numeric or missing values.")
    }
  }
  if (any(d$n < 2)) stop("All rows must have n >= 2.")
  if (any(d$sd <= 0)) stop("All rows must have sd > 0.")

  if (is.numeric(d$L)) {
    non_na <- !is.na(d$L)
    if (any(non_na) && any(d$L[non_na] != as.integer(d$L[non_na]))) {
      stop("Column 'L' contains non-integer values.")
    }
    d$L <- suppressWarnings(as.integer(d$L))
  }
  if ("true_L" %in% names(d) && is.numeric(d$true_L)) {
    non_na <- !is.na(d$true_L)
    if (any(non_na) && any(d$true_L[non_na] != as.integer(d$true_L[non_na]))) {
      stop("Column 'true_L' contains non-integer values.")
    }
    d$true_L <- suppressWarnings(as.integer(d$true_L))
  }
  d
}
