#####################
# UTILITY FUNCTIONS #
#####################

#####
# A # UTILITY FUNCTIONS FOR USE WITH ESTIMATION ALGORITHMS
#####

# turns hessian from lder2 to matrix
convert_to_hessian <- function(h){

  # total number of diagonal elements
  n <- floor(sqrt(2 * length(h)))

  # replace diags/off-diags
  H               <- diag(h[seq_len(n)])
  H[lower.tri(H)] <- h[-seq_len(n)]
  H[upper.tri(H)] <- H[lower.tri(H)]

  return(H)
} # END convert_to_hessian FUNCTION

# calculate loglik (-1 for minimization) for one person
loglik_mupp_rank_with_prior1 <- function(thetas,
                                         resp,
                                         params,
                                         items,
                                         prior_mean,
                                         prior_sd){
  loglik_mupp_rank_impl(thetas = rbind(thetas),
                        params = params,
                        items  = items,
                        picked_orders = resp) %>%
  sum(na.rm = TRUE) %>%
  multiply_by(-1)
} # END loglik_mupp_rank_with_prior FUNCTION

# calculate lder1 (-1 for minimization) for one person
lder1_mupp_rank_with_prior1 <- function(thetas,
                                        resp,
                                        params,
                                        items,
                                        prior_mean,
                                        prior_sd){
  -1 * c(lder1_mupp_rank_impl(thetas = rbind(thetas),
                              params = params,
                              items  = items,
                              picked_orders = resp))
} # END lder1_mupp_rank_with_prior1 FUNCTION

# calculate lder2 (-1 for minimization) for one person
lder2_mupp_rank_with_prior1 <- function(thetas,
                                        resp,
                                        params,
                                        items,
                                        prior_mean,
                                        prior_sd){
  -1 * convert_to_hessian(lder2_mupp_rank_impl(thetas = rbind(thetas),
                                               params = params,
                                               items  = items,
                                               picked_orders = resp))
} # END lder2_mupp_rank_with_prior1 FUNCTION

# line search to find ideal alpha for estimation algorithm
# (satisfying wolfe conditions, ideally)
wolfe_line_search <- function(fun,
                              dfun,
                              x, p,
                              c1 = .1,
                              c2 = .9,
                              ...){

  is_not_wolfe_f <- function(alpha){
    check <- c(fun(x + alpha * p, ...) > fx + c1 * alpha * fx_dir)
    is.na(check) || check
  } # END is_not_wolfe_f FUNCTION

  is_not_wolfe_d <- function(alpha){
    check <- sum(p * dfun(x + alpha * p, ...)) < c2 * fx_dir
    is.na(check) || check
  } # END is_not_wolfe_d FUNCTION

  # initialize parameter values
  alpha       <- 1
  alpha_left  <- 0
  alpha_right <- Inf

  # calculate function values and direction at starting point
  fx          <- fun(x, ...)
  fx_dir      <- sum(p * dfun(x, ...))
  iter        <- 1

  # ALGORITHM #
  # https://sites.math.washington.edu/~burke/crs/408/notes/nlp/line.pdf (p. 8)
  # https://wiki.math.ntnu.no/_media/tma4180/2015v/bfgs.m (same thing, but less pithy!)
  while(iter < 10){
    if(is_not_wolfe_f(alpha)){
      alpha_right  <- alpha
      alpha        <- (alpha_left + alpha_right) / 2
    } else if(is_not_wolfe_d(alpha)){
      alpha_left   <- alpha
      alpha        <- (alpha_left + c(3 * alpha_left, alpha_right)[is.finite(alpha_right) + 1]) / 2
    } else{
      break
    } # END ifelse STATEMENTS
    iter <- iter + 1
  } # END while LOOP

  return(alpha)

} # END wolfe_line_search FUNCTION

#####
# B # ARGUMENT CHECKS
#####

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

#####
# C # FIXING DATA FRAMES
#####

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

#####
# B # Generic Utility Functions
#####

#' @importFrom roperators
#'             is.irregular_list
#'             any_bad_for_calcs
arrange_by_vars <- function(df,
                            vars       = NULL,
                            decreasing = FALSE){

  # if df is a list, apply arrange_by_vars to the list
  if(is.irregular_list(df)){
    return(lapply(df,
                  FUN         = arrange_by_vars,
                  vars        = vars,
                  decreasing  = decreasing))
  } # END if STATEMENT

  ## checking arguments ##

  # if data frame is NULL, return NULL
  if(!length(df)){
    return(NULL)
  } # END if STATEMENT

  # make sure df is a data.frame #
  stopifnot(is.data.frame(df))

  # if vars is NULL, return df
  if(any_bad_for_calcs(vars)){
    return(df)
  } # END if STATEMENT

  # make sure vars is a numeric or character string
  stopifnot(is.numeric(vars) || is.character(vars))

  # make sure decreasing is a logical vector
  stopifnot(is.logical(decreasing))

  ## determining vars ##
  vars <- unique(vars)

  # if vars is numeric, round it and keep only vars in df
  if(is.numeric(vars)){
    vars <- as.integer(vars)
    vars <- intersect(vars,
                      seq_len(ncol(df)))
  } # END if STATEMENT

  # if vars is character, turn it into indices of df columns
  if(is.character(vars)){
    vars <- intersect(vars,
                      names(df))
    vars <- match(vars, names(df))
  } # END if STATEMENT

  # if we have no variables, return the df
  if(!length(vars)){
    return(df)
  } # END if STATEMENT

  ## determining decreasing ##

  # should be the same length as vars
  decreasing <- rep(decreasing,
                    length.out = length(vars))

  # - pull out columns of df
  # - order all of those columns based on decreasing
  # - arrange df based on order
  lapply(vars,
         function(ind)
           df[[ind]]) %>%
  do.call(what = function(...)
            order(...,
                  decreasing = decreasing,
                  method     = "radix")) %>%
  "["(df, ., , drop = FALSE)
} # END arrange_by_vars FUNCTION
