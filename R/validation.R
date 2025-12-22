#' Input Validation Helpers
#'
#' Internal helper functions for validating input data and parameters.
#' These functions emit warnings but allow execution to continue (permissive validation).
#'
#' @name validation
#' @keywords internal
NULL

#' Check that a data.frame has required columns
#'
#' @param df A data.frame to check
#' @param required_cols Character vector of required column names
#' @param df_name Name of the data.frame for error messages
#' @return TRUE if valid, FALSE if columns are missing (with warning)
#' @keywords internal
check_required_columns <- function(df, required_cols, df_name = "data") {
  if (!is.data.frame(df)) {
    warning(sprintf("%s must be a data.frame, got: %s", df_name, class(df)[1]))
    return(FALSE)
  }

  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    warning(sprintf("%s is missing required columns: %s", df_name, paste(missing, collapse = ", ")))
    return(FALSE)
  }
  TRUE
}

#' Check that a value is a probability (between 0 and 1)
#'
#' @param x A numeric value to check
#' @param name Parameter name for error messages
#' @param strict If TRUE, value must be strictly between 0 and 1; if FALSE, endpoints allowed
#' @return TRUE if valid, FALSE otherwise (with warning)
#' @keywords internal
check_probability <- function(x, name, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    warning(sprintf("%s must be a single numeric value, got: %s", name,
                    paste(class(x), collapse = ", ")))
    return(FALSE)
  }

  if (strict) {
    if (x <= 0 || x >= 1) {
      warning(sprintf("%s must be strictly between 0 and 1, got: %s", name, x))
      return(FALSE)
    }
  } else {
    if (x < 0 || x > 1) {
      warning(sprintf("%s must be between 0 and 1, got: %s", name, x))
      return(FALSE)
    }
  }
  TRUE
}

#' Check that a value is a positive integer
#'
#' @param x A numeric value to check
#' @param name Parameter name for error messages
#' @param allow_zero If TRUE, zero is allowed
#' @return TRUE if valid, FALSE otherwise (with warning)
#' @keywords internal
check_positive_integer <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    warning(sprintf("%s must be a single numeric value, got: %s", name,
                    paste(class(x), collapse = ", ")))
    return(FALSE)
  }

  if (x != floor(x)) {
    warning(sprintf("%s must be an integer, got: %s", name, x))
    return(FALSE)
  }

  if (allow_zero) {
    if (x < 0) {
      warning(sprintf("%s must be a non-negative integer, got: %s", name, x))
      return(FALSE)
    }
  } else {
    if (x < 1) {
      warning(sprintf("%s must be a positive integer (>= 1), got: %s", name, x))
      return(FALSE)
    }
  }
  TRUE
}

#' Check that a value is positive (numeric, > 0)
#'
#' @param x A numeric value to check
#' @param name Parameter name for error messages
#' @param allow_zero If TRUE, zero is allowed
#' @return TRUE if valid, FALSE otherwise (with warning)
#' @keywords internal
check_positive_numeric <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
    warning(sprintf("%s must be a single numeric value, got: %s", name,
                    paste(class(x), collapse = ", ")))
    return(FALSE)
  }

  if (allow_zero) {
    if (x < 0) {
      warning(sprintf("%s must be non-negative, got: %s", name, x))
      return(FALSE)
    }
  } else {
    if (x <= 0) {
      warning(sprintf("%s must be positive (> 0), got: %s", name, x))
      return(FALSE)
    }
  }
  TRUE
}

#' Check that a value is one of the allowed options
#'
#' @param x Value to check
#' @param allowed Character vector of allowed values
#' @param name Parameter name for error messages
#' @return TRUE if valid, FALSE otherwise (with warning)
#' @keywords internal
check_one_of <- function(x, allowed, name) {
  if (length(x) != 1) {
    warning(sprintf("%s must be a single value, got length: %d", name, length(x)))
    return(FALSE)
  }

  if (!x %in% allowed) {
    warning(sprintf("%s must be one of: %s. Got: '%s'",
                    name, paste(allowed, collapse = ", "), x))
    return(FALSE)
  }
  TRUE
}

#' Validate CNbins data.frame structure
#'
#' Checks that a CNbins data.frame has the required columns for processing.
#'
#' @param CNbins A data.frame with copy number bin data
#' @return TRUE if valid, FALSE otherwise (with warnings for issues found)
#' @keywords internal
validate_cnbins <- function(CNbins) {
  required_cols <- c("cell_id", "chr", "start", "end", "state", "copy")
  check_required_columns(CNbins, required_cols, "CNbins")
}

#' Validate haplotypes data.frame structure
#'
#' Checks that a haplotypes data.frame has the required columns for processing.
#' This is for formatted haplotypes (after format_haplotypes or format_haplotypes_dlp).
#'
#' @param haplotypes A data.frame with haplotype data
#' @param formatted If TRUE, check for formatted haplotype columns (allele1, allele0);
#'   if FALSE, check for raw haplotype columns (allele_id, readcount)
#' @return TRUE if valid, FALSE otherwise (with warnings for issues found)
#' @keywords internal
validate_haplotypes <- function(haplotypes, formatted = TRUE) {
  base_cols <- c("cell_id", "chr", "start", "end", "hap_label")

  if (formatted) {
    required_cols <- c(base_cols, "allele1", "allele0", "totalcounts")
  } else {
    required_cols <- c(base_cols, "allele_id", "readcount")
  }

  check_required_columns(haplotypes, required_cols, "haplotypes")
}

#' Validate HSCN/ASCN output structure
#'
#' Checks that an HSCN or ASCN result has expected columns.
#'
#' @param hscn A data.frame with haplotype/allele-specific copy number results
#' @return TRUE if valid, FALSE otherwise (with warnings)
#' @keywords internal
validate_hscn <- function(hscn) {
  required_cols <- c("cell_id", "chr", "start", "end", "state", "A", "B",
                     "state_AS_phased", "state_min", "LOH")
  check_required_columns(hscn, required_cols, "hscn")
}
