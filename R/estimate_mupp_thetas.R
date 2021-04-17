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
#' @return a list of [estimates, vars, hessian, loglik, iters] for MLE estimation
#'         or [estimates, sds] for MCMC estimation
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
#'       \item{n_starts}{number of random starts at different theta vectors (using
#'                       latin hypercube samples. Defaults to 4.)}
#'     }
#'   }
#'
#'   \item{For MCMC, control parameters include
#'
#'     \describe{
#'       \item{n_iters}{total number of iterations.}
#'       \item{n_burnin}{number of iterations to throw away when calculating
#'                       summary statistics.}
#'       \item{step_size_sd}{the standard deviation of the step size for subsequent
#'                           Metropolis-Hastings draws.}
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
#' \dontrun{
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
#' # estimating thetas using BFGS algorithm (one start for comparison purposes)
#' est_thetas_mle  <- estimate_mupp_thetas(resp    = resp$resp,
#'                                         items   = params$items,
#'                                         method  = "bfgs",
#'                                         control = list(n_starts = 1))
#'
#' # correlating (super high correlations!)
#' diag(cor(thetas, est_thetas_mle$estimates))
#'
#' # estimating thetas using MCMC algorithm
#' est_thetas_mcmc <- estimate_mupp_thetas(resp   = resp$resp,
#'                                         items  = params$items,
#'                                         method  = "mcmc",
#'                                         control = list(n_iters  = 1000,
#'                                                        n_burnin = 500))
#'
#' # correlating with MLE thetas (even higher correlations!)
#' diag(cor(est_thetas_mle$estimates, est_thetas_mcmc$estimates))
#' }
#'
#' @importFrom roperators
#'             "%ni%"
#' @importFrom magrittr
#'             "%>%" "%<>%" set_rownames multiply_by
#' @importFrom data.table
#'             dcast as.data.table
#' @importFrom lhs
#'             optimumLHS
#' @export
estimate_mupp_thetas <- function(resp,
                                 items,
                                 method  = c("bfgs", "MCMC"),
                                 control = list(),
                                 ...){

  # determine resp/items/params for algorithm
  est_args      <- estimate_mupp_header_(resp, items,
                                         type = "person")

  # run algorithm
  algorithm_fun <- switch(method,
                          bfgs = estimate_mupp_thetas_mle,
                          mcmc = , MCMC  = estimate_mupp_thetas_mcmc,
                          stop(method, " method not implemented at this time."))

  out           <- algorithm_fun(resp    = est_args$resp_adj,
                                 params  = est_args$params_adj,
                                 items   = est_args$items_adj,
                                 method  = method,
                                 control = control,
                                 ...)

  # return
  return(out)
} # END estimate_mupp_thetas FUNCTION


# MAXIMUM LIKELIHOOD ESTIMATION #
#' @importFrom utils
#'             modifyList
#'             txtProgressBar setTxtProgressBar
estimate_mupp_thetas_mle <- function(resp,
                                     params,
                                     items,
                                     method  = "bfgs",
                                     control = list(),
                                     ...){


  # converting everything to a matrix (to work in C++ algorithm)
  resp   <- as.matrix(resp)
  params <- as.matrix(params)
  items  <- as.matrix(items)

  # indicate basic things (total number of items/dimensions)
  n_persons     <- nrow(resp)
  n_dims        <- max(items[ , 3])

  # determining maximum number of "dimensions" on a single item
  max_item_dims <- max(table(items[ , 1]))

  # making sure number of dimensions is OK
  if(max_item_dims > 2){
    stop("MLE estimation not yet implemented for MUPP-RANK models. ",
         "Please use MCMC for theta/parameter estimation.",
         call. = FALSE)
  } # END if STATEMENT

  # create matrix of thetas based on maximum dimension
  out       <- vector(mode   = "list",
                      length = n_persons)

  # update control
  control_default <- list(prior_mean = 0,
                          prior_sd   = 1,
                          eps        = 1e-07,
                          max_iters  = 100,
                          n_starts   = 4)
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
  hessian  <- lapply(out, "[[", "hessian") %>%
              lapply(FUN = c) %>%
              do.call(what = rbind)
  iters    <- sapply(out, "[[", "iters")
  loglik   <- sapply(out, "[[", "loglik")

  # determine actual variances
  var     <- lapply(seq_len(nrow(thetas)),
                    FUN = function(i){
                      H <- lder2_mupp_rank_with_prior1(thetas = thetas[i, ],
                                                       resp   = rbind(resp[i, ]),
                                                       params = params,
                                                       items  = items)
                      c(solve(H))
                    }) %>%
             do.call(what = rbind)

  return(list(estimates = thetas,
              vars      = var,
              hessian   = hessian,
              loglik    = loglik,
              iters     = iters))
} # END estimate_mupp_thetas_mle FUNCTION


# MAXIMUM LIKELIHOOD ESTIMATION #
estimate_mupp_thetas_mcmc <- function(resp,
                                      params,
                                      items,
                                      control = list(),
                                      ...){

  # run the mcmc algorithm
  out <- estimate_mupp_params_mcmc(resp    = resp,
                                   items   = items,
                                   control = control,
                                   initial_params = as.list(as.data.frame(params)),
                                   fixed_params   = colnames(params))

  # pull out parameters
  out <- list(estimates = out$means$thetas,
              sds       = out$sds$thetas)


  return(out)
} # END estimate_mupp_thetas_mcmc FUNCTION
