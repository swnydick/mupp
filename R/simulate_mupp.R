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
#' @importFrom magrittr "%>%" set_rownames
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
  # (sample.int used in case dims is length 1 to prevent switching methods)
  items_n_dims <- dims[-1] %>%
                  "["(sample.int(length(.),
                                 size    = n_items,
                                 replace = TRUE))

  # which dim each item loads on
  items        <- lapply(X   = seq_along(items_n_dims),
                         FUN = function(i){
                           item_dims  <- items_n_dims[i]
                           data.frame(item      = i,
                                      statement = NA,
                                      dim       = sort(sample(x       = dims,
                                                              size    = item_dims,
                                                              replace = FALSE)))
                         }) %>%
                  do.call(what = rbind)
  items        <- transform(items,
                            statement = seq_along(statement))

  # the total number of parameters
  n_params     <- nrow(items)

  # simulating alpha/delta/eta
  items        <- transform(items,
                            alpha = r_alpha_prior(n_params),
                            delta = r_delta_prior(n_params),
                            tau   = r_tau_prior(n_params))

  ## return ##
  return(list(persons = persons,
              items   = items))
} # END simulate_mupp_params FUNCTION


#' Simulate MUPP Responses
#'
#' Generate responses that can be used for the MUPP model.
#'
#' @param persons persons data.frame with column names
#'        (person, dim, theta)
#' @param items items data.frame with column names
#'        (item, statement, dim, alpha, delta, tau)
#'
#' @return a data.frame of [person x item x response pattern]
#'
#' @details The persons and items df needs to look identical to that coming from
#'          \code{\link{simulate_mupp_params}} or else this function will not work.
#'
#' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
#'
#' @importFrom kfhelperfuns arrange_by_vars "%ni%"
#' @importFrom magrittr "%>%" "%<>%" set_rownames
#' @importFrom data.table dcast as.data.table
#' @export
simulate_mupp_resp <- function(persons,
                               items){

  # helper functions #

  # a # function to ensure argument has appropriate names
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

  # converting to lowercase
  persons %<>% setNames(tolower(names(.)))
  items   %<>% setNames(tolower(names(.)))

  # adding "statement" if it is missing
  if("statement" %ni% names(items)){
    items$statement <- seq_len(nrow(items))
  } # END if STATEMENT

  # pulling out template person/item
  template <- lapply(simulate_mupp_params(),
                     FUN = names)

  # ordering the columns
  persons  <- check_names(persons,
                          template$persons)
  items    <- check_names(items,
                          template$items)

  # fix persons / items #

  # pull out useful names
  item_names     <- names(items)
  item_name      <- item_names[1]
  statement_name <- item_names[2]
  dim_name       <- item_names[3]
  param_names    <- setdiff(item_names, c(item_name, statement_name, dim_name))

  # reshape so that [persons, theta across columns]
  f_vars   <- head(names(persons), 2)
  v_var    <- tail(names(persons), 1)
  persons  <- dcast(data      = as.data.table(persons),
                    formula   = as.formula(paste(f_vars, collapse = "~")),
                    value.var = v_var) %>%
              as.data.frame()

  # split items so that different items are different list elements
  items  %<>% split(.[[item_name]])

  # simulating responses #

  # determining probabilities
  probs      <- lapply(items,
                       FUN            = determine_mupp_probs1,
                       persons        = persons,
                       dimension_name = dim_name,
                       param_names    = param_names)

  # simulating responses
  resp       <- lapply(probs,
                       FUN     = simulate_mupp_resp1)

  # converting to item order
  all_orders <- lapply(X    = seq_len(max(sapply(items, nrow))),
                       FUN  = find_all_permutations,
                       init = 1)
  resp       <- Map(item = items,
                    resp = resp,
                    name = names(items),
                    f = function(item, resp, name){
                      data.frame(persons[1],
                                 item   = type.convert(name),
                                 n_dims = nrow(item),
                                 resp   = resp,
                                 stringsAsFactors = FALSE)
                    }) %>%
               do.call(what = rbind) %>%
               set_rownames(NULL)

  # fixing items
  items      <- do.call(what = rbind,
                        args = items) %>%
                "["(names(.) %ni% param_names) %>%
                set_rownames(NULL)

  return(list(items = items,
              resp  = resp))

} # END simulate_mupp_responses FUNCTION


determine_mupp_probs1 <- function(item,
                                  persons,
                                  dimension_name    = "dim",
                                  param_names       = c("alpha", "delta", "tau"),
                                  picked_order_name = NULL){


  # pull out dimension/params/theta
  dims   <- item[[dimension_name]]
  params <- data.matrix(item[param_names])
  thetas <- data.matrix(persons[-1])

  # pull out picked order
  if(!length(picked_order_name)){
    picked_order <- NA
  } else{
    picked_order <- item[[picked_order_name]]
  } # END if STATEMENT

  # check to make sure picked order name is in data
  if(!length(picked_order)){
    stop("picked_order_name is not in item data.frame",
         call. = FALSE)
  } # END if STATEMENT

  # calculate probability
  probs <- p_mupp_rank_impl(thetas = thetas,
                            params = params,
                            dims   = dims,
                            picked_order_id = picked_order)

  return(probs)

} # END determine_mupp_probs1 FUNCTION

simulate_mupp_resp1 <- function(probs){

  # make sure mat is a data.matrix
  probs <- data.matrix(probs)

  # converting to cumulative sum
  probs <- t(apply(probs, MARGIN = 1, FUN = cumsum))

  # simulating response for everybody
  u     <- runif(n = nrow(probs))

  rowSums(u >= probs) + 1

} # END simulate_mupp_resp FUNCTION
