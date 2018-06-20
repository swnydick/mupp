# ARGUMENT CHECKS #

# a # check to make sure all names is in df and return df in appropriate order
check_names <- function(df, names){

  # pull out argument name AND round argument
  arg_char  <- deparse(substitute(df))
  bad_names <- setdiff(names, names(df))

  if(length(bad_names)){
    stop("column ", bad_names, " is not in ", arg_char, ".",
         call. = FALSE)
  } # END if STATEMENT

  return(df[names])

} # END check_names FUNCTION


# b # ensure argument is within numeric range
check_numeric <- function(arg,
                          min_number = 1){

  # pull out argument name AND round argument
  arg_char <- deparse(substitute(arg))
  arg      <- round(arg)

  # ensure that argument is within bounds
  if(!is.numeric(arg) || length(arg) != 1 || arg < min_number || arg == Inf){
    stop(arg_char, " must be a single number greater than ", min_number,
         call. = FALSE)
  } # END if STATEMENT

  return(arg)
} # END check_numeric FUNCTION


# FIXING DATA FRAMES #

# a # fixing columns to be numeric and in order
sequence_column <- function(df,
                            column,
                            old_values = NULL){

  # all of the values/unique old values/unique new values
  old_values_all <- df[[column]]

  # allow for specifying old values
  if(is.null(old_values)){
    old_values  <- unique(old_values_all)
  } # END if STATEMENT

  new_values     <- seq_along(old_values)

  # determining new names (for all rows)
  new_values_all <- new_values[match(old_values_all, old_values)]

  # adding back to df
  df[[column]]    <- new_values_all

  return(df)

} # END sequence_column FUNCTION
