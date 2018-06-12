#' Simulate MUPP Parameters
#'
#' Generate parameters/thetas that conform to the MUPP model.
#'
#' @param n_persons integer indicating the number of persons
#' @param n_items integer indicating the number of items
#' @param n_dims integer > 2 indicating the total number of dimensions
#'
#' @return a list with items/persons that conform to the MUPP model, as expected
#'         in package functions
#'
#' @details For the purposes of parameter generation, each item will be assumed
#'          to have between 2 and n_dims dimensions, where the number of dimensions
#'          for an item is randomly generated from a discrete uniform distribution.
#'
#' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
#'
#' @importFrom kfhelperfuns arrange_by_vars
#' @importFrom magrittr "%>%"
#' @export

simulate_mupp_params <- function(n_persons = 1,
                                 n_items   = 1,
                                 n_dims    = 2){

  # helper functions #

  # a # function to ensure argument is within numeric range
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

  # argument checks #
  n_persons <- check_numeric(n_persons)
  n_items   <- check_numeric(n_items)
  n_dims    <- check_numeric(n_dims,
                             min_number = 2)

  ## persons ##

  # construct persons df
  persons    <- matrix(r_thetas_prior(n_persons * n_dims),
                       nrow = n_persons,
                       ncol = n_dims) %>%
                as.data.frame() %>%
                setNames(paste0("theta_", seq_len(n_dims)))
  persons    <- cbind(person = seq_len(n_persons), persons)

  # reshaping wide to long
  persons    <- reshape(persons,
                        sep       = "_",
                        direction = "long",
                        varying   = -1,
                        timevar   = "dim",
                        idvar     = names(persons)[1]) %>%
                arrange_by_vars(vars = names(.)[1:2]) %>%
                set_rownames(NULL)

  ## items ##

  # (uses ggum)

  # the dims across all items
  dims         <- seq_len(n_dims)

  # number of dims each item loads on
  items_n_dims <- sample(x       = dims[-1],
                         size    = n_items,
                         replace = TRUE)

  # which dim each item loads on
  items        <- lapply(X   = seq_along(items_n_dims),
                         FUN = function(i){
                           data.frame(item = i,
                                      dim  = sort(sample(x       = dims,
                                                         size    = items_n_dims,
                                                         replace = FALSE)))
                         }) %>%
                  do.call(what = rbind)

  # the total number of parameters
  n_params     <- nrow(items)

  # simulating alpha/delta/eta
  items$alpha  <- r_alpha_prior(n_params)
  items$delta  <- r_delta_prior(n_params)
  items$tau    <- r_tau_prior(n_params)

  ## return ##
  return(list(persons = persons,
              items   = items))
} # END simulate_mupp_params FUNCTION
