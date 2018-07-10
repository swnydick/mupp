#' Estimate MUPP Parameters
#'
#' Estimate MUPP statement and person parameters given item responses and item
#' properties using MCMC
#'
#' @param resp a data.frame of (at least) [person, item, resp]
#' @param items a data.frame of (at least) [item, statement, dim]
#' @param method the estimation method (MCMC is the only one that works now)
#' @param control a list of parameters to control the algorithm. See details.
#' @param ... other parameters to pass to the estimation algorithm. See details.
#'
#' @return a list of [theta, params] arrays with the third dimension indicating
#'         the iteration as well as means and sds
#'
#' @details
#'
#' \itemize{
#'   \item{For MCMC, additional parameters include
#'     \describe{
#'       \item{initial_params}{a named list with names indicating the parameter
#'                             and values indicating the starting values/initial
#'                             parameter estimates}
#'       \item{fixed_params}{a character vector with elements that are fixed
#'                           to their initial parameters for the entire estimation
#'                           algorithm}
#'     }
#'   }
#' }
#'
#' \itemize{
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
#' @examples
#' \dontrun{
#' set.seed(3452345)
#'
#' # simulate parameters and responses to the model
#' # (assumption is that params/resp will follow conventions)
#' params <- simulate_mupp_params(n_persons     = 100,
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
#' items  <- params$items
#'
#' # estimating thetas using algorithm (one start for comparison purposes)
#' est_params  <- estimate_mupp_params(resp    = resp$resp,
#'                                     items   = resp$items,
#'                                     method  = "MCMC",
#'                                     control = list(n_iters  = 1000,
#'                                                    n_burnin = 500),
#'                                     initial_params = list(delta = sign(items$delta)))
#'
#' # correlating (not great, but small iters and few people)
#' diag(cor(thetas, est_params$means$thetas))
#' cor(items$alpha, est_params$means$alpha)
#' cor(items$delta, est_params$means$delta)
#' cor(items$tau,   est_params$mean$tau)
#' }
#'
#' @importFrom kfhelperfuns arrange_by_vars "%ni%"
#' @importFrom magrittr "%>%" "%<>%" set_rownames
#' @importFrom data.table dcast as.data.table
#' @export
estimate_mupp_params <- function(resp,
                                 items,
                                 method  = "MCMC",
                                 control = list(),
                                 ...){

  # determine resp/items for algorithm
  est_args      <- estimate_mupp_header_(resp, items,
                                         type = "item")

  # run algorithm
  algorithm_fun <- switch(method,
                          mcmc = , MCMC = estimate_mupp_params_mcmc,
                          stop(method, " method not implemented at this time"))
  out           <- algorithm_fun(resp    = est_args$resp_adj,
                                 items   = est_args$items_adj,
                                 control = control,
                                 ...)

  # return
  return(out)

} # END estimate_mupp_params FUNCTION
