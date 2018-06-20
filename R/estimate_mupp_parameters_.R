#' Estimate MUPP Parameters
#'
#' Estimate MUPP statement and person parameters given item responses and item
#' properties using MCMC
#'
#' @param resp a data.frame of (at least) [person, item, resp]
#' @param items a data.frame of (at least) [item, statement, dim]
#' @param method the estimation method (MCMC is the only one that works now)
#' @param ... other parameters to add later
#'
#' @return a list of [theta, params] arrays with the third dimension indicating
#'         the iteration
#'
#' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
#'
#' @importFrom kfhelperfuns arrange_by_vars "%ni%"
#' @importFrom magrittr "%>%" "%<>%" set_rownames
#' @importFrom data.table dcast as.data.table
#' @export
estimate_mupp_params <- function(resp,
                                 items,
                                 method = "MCMC",
                                 ...){

  # argument checks #

  # converting to lowercase
  resp  %<>% setNames(tolower(names(.)))
  items %<>% setNames(tolower(names(.)))

  # pulling out template resp/item
  template <- lapply(do.call(what = simulate_mupp_resp,
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

  # run algorithm
  if(any(tolower(method) == "mcmc")){
    out <- estimate_mupp_params_mcmc(resp  = resp_adj,
                                     items = items_adj,
                                     ...)
  } else{
    stop(method, " method not implemented at this time.")
  } # END ifelse STATEMENT

  return(out)

} # END estimate_mupp_params FUNCTION


# ESTIMATION: MCMC OVERALL #
estimate_mupp_params_mcmc <- function(resp,
                                      items,
                                      n_iters      = 10000,
                                      n_burnin     = 1000,
                                      step_size_sd = .1){

  items          <- as.matrix(items)

  # INITIALIZE THETA/STATEMENT VALUES #
  current_params <- initialize_mupp_params_mcmc(resp  = resp,
                                                items = items)


  # setup output stuff
  thetas_out <- vector(mode   = "list",
                       length = n_iters)
  alphas_out <- deltas_out <- taus_out <- thetas_out

  # iterate ...

  # seems to work ... so cleanup!
  for(iter in seq_len(n_iters)){

    if(iter %% 10 == 0)
      print(iter)

    # generate next parameters
    current_params     <- estimate_mupp_params_mcmc_(resp           = resp,
                                                     items          = items,
                                                     current_params = current_params,
                                                     step_size_sd)

    # assign to output lists (automate?)
    thetas_out[[iter]] <- current_params$thetas
    alphas_out[[iter]] <- current_params$alphas
    deltas_out[[iter]] <- current_params$deltas
    taus_out[[iter]]   <- current_params$taus
  } # END for iter LOOP

  iter_remove <- seq_len(n_burnin)

  all_iters  <- list(thetas = thetas_out,
                     alphas = alphas_out,
                     deltas = deltas_out,
                     taus   = taus_out)
  keep_iters <- lapply(all_iters,
                       "[",
                       -iter_remove)

  mean_sim   <- lapply(keep_iters,
                       function(list){
                          Reduce("+", list) / length(list)
                       })

  list(all  = all_iters,
       mean = mean_sim)

}

# ESTIMATION: MCMC ITER #

estimate_mupp_params_mcmc_ <- function(resp, items, current_params, step_size_sd){

  # extracting parameters
  thetas <- current_params$thetas
  alphas <- current_params$alphas
  deltas <- current_params$deltas
  taus   <- current_params$taus

  # updating parameters
  thetas <- update_mupp_thetas_mcmc(resp, items, thetas, alphas, deltas, taus, step_size_sd)
  alphas <- update_mupp_alphas_mcmc(resp, items, thetas, alphas, deltas, taus, step_size_sd)
  deltas <- update_mupp_deltas_mcmc(resp, items, thetas, alphas, deltas, taus, step_size_sd)
  taus   <- update_mupp_taus_mcmc(  resp, items, thetas, alphas, deltas, taus, step_size_sd)

  # returning
  list(thetas = thetas,
       alphas = alphas,
       deltas = deltas,
       taus   = taus)

} # END estimate_mupp_params_mcmc_

# ESTIMATION: MCMC INDIVIDUAL UPDATE FUNCTIONS #

# Note: FIX
#       pass likelihood so only calculating new likelihood
#       automate better stuffy stuff

# THETAS #
update_mupp_thetas_mcmc <- function(resp,
                                    items,
                                    thetas,
                                    alphas,
                                    deltas,
                                    taus,
                                    step_size_sd){

  # generate new thetas
  thetas   %<>% generate_new_params(step_size_sd)

  # generate likelihoods
  logliks    <- lapply(X      = thetas,
                       FUN    = function(x){
                         loglik_mupp_rank_impl(thetas = x,
                                               params = cbind(alphas,
                                                              deltas,
                                                              taus),
                                               items  = items,
                                               picked_orders = resp)
                       }) %>%
                lapply(FUN    = rowSums)
  priors     <- lapply(X      = thetas,
                       FUN    = function(mat){
                         mat[] <- log(d_thetas_prior(mat))
                         rowSums(mat)
                       })

  # compare
  thetas_new <- update_with_metrop(values_old = thetas[[1]],
                                   values_new = thetas[[2]],
                                   loglik_old = logliks[[1]],
                                   loglik_new = logliks[[2]],
                                   priors_old = priors[[1]],
                                   priors_new = priors[[2]])

  # return
  return(thetas_new)

} # END update_mupp_thetas_mcmc FUNCTION

# ALPHAS #
update_mupp_alphas_mcmc <- function(resp,
                                    items,
                                    thetas,
                                    alphas,
                                    deltas,
                                    taus,
                                    step_size_sd){

  # generate new alphas
  alphas   %<>% generate_new_params(step_size_sd)

  # generate likelihoods
  logliks    <- lapply(X      = alphas,
                       FUN    = function(x){
                         loglik_mupp_rank_impl(thetas = thetas,
                                               params = cbind(x,
                                                              deltas,
                                                              taus),
                                               items  = items,
                                               picked_orders = resp)
                       }) %>%
                lapply(FUN    = colSums)
  priors     <- lapply(X      = alphas,
                       FUN    = function(vec){
                         vec <- log(d_alpha_prior(vec))
                         tapply(vec, items[ , 1], sum)
                       })

  item_list  <- lapply(1:2,
                       FUN = function(i){
                         rep(i, length(priors[[1]]))
                       })

  # compare
  items_new  <- update_with_metrop(values_old = item_list[[1]],
                                   values_new = item_list[[2]],
                                   loglik_old = logliks[[1]],
                                   loglik_new = logliks[[2]],
                                   priors_old = priors[[1]],
                                   priors_new = priors[[2]])
  state_new  <- items_new[items[ , 1]]
  alphas_new <- ifelse(state_new == 1,
                       yes       = alphas[[1]],
                       no        = alphas[[2]])

  # return
  return(alphas_new)

} # END update_mupp_thetas_mcmc FUNCTION


# DELTAS #
update_mupp_deltas_mcmc <- function(resp,
                                    items,
                                    thetas,
                                    alphas,
                                    deltas,
                                    taus,
                                    step_size_sd){

  # generate new deltas
  deltas   %<>% generate_new_params(step_size_sd)

  # generate likelihoods
  logliks    <- lapply(X      = deltas,
                       FUN    = function(x){
                         loglik_mupp_rank_impl(thetas = thetas,
                                               params = cbind(alphas,
                                                              x,
                                                              taus),
                                               items  = items,
                                               picked_orders = resp)
                       }) %>%
                lapply(FUN    = colSums)
  priors     <- lapply(X      = deltas,
                       FUN    = function(vec){
                         vec <- log(d_delta_prior(vec))
                         tapply(vec, items[ , 1], sum)
                       })

  item_list  <- lapply(1:2,
                       FUN = function(i){
                         rep(i, length(priors[[1]]))
                       })

  # compare
  items_new  <- update_with_metrop(values_old = item_list[[1]],
                                   values_new = item_list[[2]],
                                   loglik_old = logliks[[1]],
                                   loglik_new = logliks[[2]],
                                   priors_old = priors[[1]],
                                   priors_new = priors[[2]])
  state_new  <- items_new[items[ , 1]]
  deltas_new <- ifelse(state_new == 1,
                       yes       = deltas[[1]],
                       no        = deltas[[2]])

  # return
  return(deltas_new)

} # END update_mupp_thetas_mcmc FUNCTION


# DELTAS #
update_mupp_taus_mcmc <- function(resp,
                                  items,
                                  thetas,
                                  alphas,
                                  deltas,
                                  taus,
                                  step_size_sd){

  # generate new deltas
  taus     %<>% generate_new_params(step_size_sd)

  # generate likelihoods
  logliks    <- lapply(X      = taus,
                       FUN    = function(x){
                         loglik_mupp_rank_impl(thetas = thetas,
                                               params = cbind(alphas,
                                                              deltas,
                                                              x),
                                               items  = items,
                                               picked_orders = resp)
                       }) %>%
                lapply(FUN    = colSums)
  priors     <- lapply(X      = taus,
                       FUN    = function(vec){
                         vec <- log(d_tau_prior(vec))
                         tapply(vec, items[ , 1], sum)
                       })

  item_list  <- lapply(1:2,
                       FUN = function(i){
                         rep(i, length(priors[[1]]))
                       })

  # compare
  items_new  <- update_with_metrop(values_old = item_list[[1]],
                                   values_new = item_list[[2]],
                                   loglik_old = logliks[[1]],
                                   loglik_new = logliks[[2]],
                                   priors_old = priors[[1]],
                                   priors_new = priors[[2]])
  state_new  <- items_new[items[ , 1]]
  taus_new   <- ifelse(state_new == 1,
                       yes       = taus[[1]],
                       no        = taus[[2]])

  # return
  return(taus_new)

} # END update_mupp_thetas_mcmc FUNCTION

# INITIALIZATION #
initialize_mupp_params_mcmc <- function(resp, items){

  # determine number of persons, dimensions, statements
  n_persons    <- nrow(resp)
  n_dims       <- max(unique(items[ , 3]))
  n_statements <- length(unique(items[ , 2]))

  # generating parameters
  thetas       <- matrix(r_thetas_prior(n = n_persons * n_dims),
                         nrow = n_persons,
                         ncol = n_dims)
  alphas       <- r_alpha_prior(n = n_statements)
  deltas       <- r_delta_prior(n = n_statements)
  taus         <- r_tau_prior(n = n_statements)

  # returnining results
  return(list(thetas = thetas,
              alphas = alphas,
              deltas = deltas,
              taus   = taus))

} # END initialize_mupp_params_mcmc FUNCTION

# HELPER FUNCTIONS #
generate_new_params <- function(x, sd = .2){
  list(old = x,
       new = x + rnorm(length(x), sd = sd))
} # END generate_new_params FUNCTION

update_with_metrop <- function(values_old,
                               values_new,
                               loglik_old,
                               loglik_new,
                               priors_old,
                               priors_new){

  # comparison
  accept_prob <- exp(loglik_new - loglik_old + priors_new - priors_old) %>%
                 "[<-"(is.na(.), value = 0)
  u           <- runif(length(accept_prob))

  accept      <- u < accept_prob

  # replacing
  if(is.matrix(values_old)){
    values_old[accept, ] <- values_new[accept, ]
  } else{
    values_old[accept]   <- values_new[accept]
  } # END ifelse STATEMENT

  return(values_old)
} # END update_with_metrop FUNCTION
