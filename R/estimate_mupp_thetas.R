#' Estimate MUPP Thetas
#'
#' Estimate MUPP person parameters given item responses, item parameters, and
#' item properties.
#'
#' @param resp a data.frame of (at least) [person, item, resp]
#' @param items a data.frame of (at least) [item, statement, dim, alpha, delta, tau]
#' @param method the estimation method ("bfgs", "MCMC")
#' @param control a list of parameters to control the algorithm. See details.
#' @param ... other parameters to pass to the method
#'
#' @return a list of [thetas, vars, hessians, logliks, iters]
#'
#' @details Method BFGS uses a modified BFGS algorithm, based on Li and Fukushima (2001),
#'          which includes modifications of the BFGS estimated hessian with demonstrated
#'          global convergence properties.
#'
#' \itemize{
#'   \item{For BFGS, control parameters include
#'
#'     \describe{
#'       \item{eps}{convergence criterion, will converge if
#'                  length of gradient is shorter than this value. Defaults
#'                  to 1e-07.}
#'       \item{max_iters}{maximum number of iterations before convergence. Defaults
#'                        to 100.}
#'       \item{starts}{number of random starts at different theta vectors (using
#'                     latin hypercube samples. Defaults to 4.)}
#'     }
#'   }
#' }
#'
#' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
#'
#' @references
#' \itemize{
#'   \item{Li & Dong-Hui (2001). A modified BFGS method and its global
#'         convergence in nonconvex minimization. Journal of Computational and
#'         Applied Mathematics, 129, 15-35.}
#' }
#'
#' @seealso \code{\link{optimumLHS}}
#'
#' @examples
#' set.seed(23523)
#'
#' # simulate parameters and responses to the model
#' # (assumption is that params/resp will follow conventions)
#' params <- simulate_mupp_params(n_persons     = 10,
#'                                n_items       = 300,
#'                                n_dims        = 9,
#'                                max_item_dims = 2,
#'                                unidim_items  = TRUE)
#' resp   <- do.call(simulate_mupp_resp,
#'                   params)
#'
#' # thetas for comparison
#' thetas <- tidyr::spread(params$persons,
#'                         key   = "dim",
#'                         value = "theta")[ , -1]
#'
#' # estimating thetas using algorithm (one start for comparison purposes)
#' est_thetas  <- estimate_mupp_thetas(resp    = resp$resp,
#'                                     items   = params$items,
#'                                     method  = "bfgs",
#'                                     control = list(starts = 1))
#'
#' # correlating (super high correlations!)
#' diag(cor(thetas, est_thetas$thetas))
#'
#' @importFrom kfhelperfuns arrange_by_vars "%ni%"
#' @importFrom magrittr "%>%" "%<>%" set_rownames multiply_by
#' @importFrom data.table dcast as.data.table
#' @importFrom lhs optimumLHS
#' @export
estimate_mupp_thetas <- function(resp,
                                 items,
                                 method  = c("bfgs", "MCMC"),
                                 control = list(),
                                 ...){

  # argument checks #

  # fix method
  method   <- match.arg(method)

  # converting to lowercase
  resp   %<>% setNames(tolower(names(.)))
  items  %<>% setNames(tolower(names(.)))

  # pulling out template resp/item
  template <- lapply(do.call(what = function(...){
                               c(simulate_mupp_resp(...)["resp"], list(...))
                             },
                             args = simulate_mupp_params()),
                     FUN = names)

  # ordering the columns
  resp   %<>% check_names(template$resp)
  items  %<>% check_names(template$items)

  # extracting variable names

  # overall names
  resp_names     <- names(resp)
  item_names     <- names(items)

  # individual names
  person_name    <- head(resp_names, 1)
  item_name      <- head(item_names, 1)
  statement_name <- item_names[2]
  resp_name      <- tail(resp_names, 1)

  # pulling out old item/statement/person and setting new item to index
  items_adj      <- sequence_column(df     = items,
                                    column = item_name) %>%
                    sequence_column(column = statement_name)
  resp_adj       <- sequence_column(df     = resp,
                                    column = person_name) %>%
                    sequence_column(column = item_name,
                                    old_values = unique(items[[item_name]]))

  # transforming response to be of the appropriate form
  cast_resp     <- as.formula(paste0(person_name, "~", item_name))

  resp_adj      <- dcast(data      = as.data.table(resp_adj),
                         formula   = cast_resp,
                         value.var = resp_name)[ , -1] %>%
                   as.matrix()

  # updating to items/params/resp
  params_adj    <- as.matrix(items_adj[4:ncol(items_adj)])
  items_adj     <- as.matrix(items_adj[seq_len(3)])

  # run algorithm
  algorithm_fun <- switch(method,
                          bfgs  = estimate_mupp_thetas_,
                          MCMC  = ,
                          stop(method, " method not implemented at this time."))

  out           <- algorithm_fun(resp    = resp_adj,
                                 params  = params_adj,
                                 items   = items_adj,
                                 method  = method,
                                 control = control,
                                 ...)

  # return
  return(out)

} # END estimate_mupp_thetas FUNCTION


estimate_mupp_thetas_ <- function(resp,
                                  params,
                                  items,
                                  method  = "bfgs",
                                  control = list(),
                                  ...){


  # converting everything to a matrix (to work in C++ algorithm)
  resp   <- as.matrix(resp)
  params <- as.matrix(params)
  items  <- as.matrix(items)

  # indicate basic things
  n_persons <- nrow(resp)
  n_dims    <- max(items[ , 3])

  # create matrix of thetas based on maximum dimension
  out       <- vector(mode   = "list",
                      length = n_persons)

  # update control
  control_default <- list(prior_mean = 0,
                          prior_sd   = 1,
                          eps        = 1e-07,
                          max_iters  = 100,
                          starts     = 4)
  control         <- modifyList(control_default,
                                control)

  # update algorithm
  estimation_fun   <- switch(method,
                             bfgs = estimate_mupp_thetas_bfgs,
                             stop(method, " method not implemented at this time."))

  # progress bar #
  pb <- txtProgressBar(max   = n_persons,
                       char  = "|",
                       style = 3)

  # estimating for each person
  for(person in seq_len(n_persons)){

    # updating arguments
    args          <- c(list(resp   = resp[person, , drop = FALSE],
                            params = params,
                            items  = items,
                            n_dims = n_dims),
                       control)

    # estimating theta
    out[[person]] <- do.call(what = estimation_fun,
                             args = args)

    # progress bar #
    setTxtProgressBar(pb    = pb,
                      value = person)
  } # END for person LOOP

  # progress bar #
  close(pb)

  # update thetas/hessians/iters/logliks to be everything in out bounded together
  thetas   <- lapply(out, "[[", "thetas") %>%
              do.call(what = rbind)
  hessians <- lapply(out, "[[", "hessian") %>%
              lapply(FUN = c) %>%
              do.call(what = rbind)
  iters    <- sapply(out, "[[", "iters")
  logliks  <- sapply(out, "[[", "loglik")

  # determine actual variances
  vars     <- lapply(seq_len(nrow(thetas)),
                     FUN = function(i){
                       H <- lder2_mupp_rank_with_prior1(thetas = thetas[i, ],
                                                        resp   = rbind(resp[i, ]),
                                                        params = params,
                                                        items  = items)
                       c(solve(H))
                     }) %>%
              do.call(what = rbind)

  return(list(thetas   = thetas,
              vars     = vars,
              hessians = hessians,
              logliks  = logliks,
              iters    = iters))

} # END estimate_mupp_thetas_ FUNCTION
