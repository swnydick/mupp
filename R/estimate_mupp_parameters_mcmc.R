# ESTIMATION: MCMC #

# OVERALL #
estimate_mupp_params_mcmc <- function(resp,
                                      items,
                                      n_iters      = 10000,
                                      n_burnin     = 1000,
                                      step_size_sd = .1,
                                      delta_sign   = NULL){

  # converting items to a matrix (so that it works in the C++ algorithm)
  items <- as.matrix(items)

  # INITIALIZE THETA/STATEMENT VALUES #
  current_params  <- initialize_mupp_params_mcmc(resp  = resp,
                                                 items = items,
                                                 delta_sign = delta_sign)

  # environment to store mcmc arguments and current loglikelihood
  mcmc_temp_envir <- new.env()

  # adding argument to the environment
  assign(x     = "arguments",
         value = modifyList(x   = current_params,
                            val = list(items        = items,
                                       resp         = resp,
                                       step_size_sd = step_size_sd,
                                       mcmc_envir   = mcmc_temp_envir)),
         envir = mcmc_temp_envir)

  # adding initial loglikelihood to the environment
  assign(x     = "loglik",
         value = do.call(what = loglik_mupp_rank_mcmc,
                         args = mcmc_temp_envir$arguments),
         envir = mcmc_temp_envir)


  # setup output stuff
  all_params     <- names(current_params)
  all_mcmc_steps <- lapply(all_params %>% setNames(., .),
                           FUN = function(x){
                             vector(mode   = "list",
                                    length = n_iters)
                           })

  # ITERATE #

  # progress bar #
  pb <- txtProgressBar(max   = n_iters,
                       char  = "mcmc",
                       style = 3)

  # iteration: mcmcmcmcmc #
  for(iter in seq_len(n_iters)){

    # generate next parameters
    estimate_mupp_params_mcmc_(mcmc_envir = mcmc_temp_envir)

    # assign to output lists
    for(param in names(current_params)){
      all_mcmc_steps[[param]][[iter]] <- mcmc_temp_envir$arguments[[param]]
    } # END for param LOOP

    # progress bar #
    setTxtProgressBar(pb    = pb,
                      value = iter)
  } # END for iter LOOP

  # progress bar #
  close(pb)

  # removing burnin
  keep_mcmc_steps <- all_mcmc_steps
  if(n_burnin > 0){
    burnin          <- seq_len(n_burnin)
    keep_mcmc_steps <- lapply(keep_mcmc_steps,
                              FUN = "[",
                              -burnin)
  } # END if STATEMENT

  # calcuating mean/sd of iterations (for posterity)
  # ... alternative to Reduce as it's very slow ??
  mean_sims       <- lapply(keep_mcmc_steps,
                            FUN = function(sims){
                              mns <- Reduce("+", sims) / length(sims)
                            })
  sd_sims         <- Map(sims  = keep_mcmc_steps,
                         means = mean_sims,
                         f     = function(sims, means){
                           ss  <- Reduce("+", lapply(sims, "^", 2))
                           len <- length(sims)
                           sds <- sqrt((ss / len - means^2) * (len - 1) / len)
                         })

  list(all  = all_mcmc_steps,
       mean = mean_sims,
       sd   = sd_sims)
} # END estimate_mupp_params_mcmc FUNCTION

# SINGLE ITERATION #

estimate_mupp_params_mcmc_ <- function(mcmc_envir){

  # for each param:
  #  - find the update function
  #  - update the param (using arguments) and assign back to the environment
  for(params in c("thetas", "alphas", "deltas", "taus")){
    mcmc_update_fun                <- paste("update_mupp", params, "mcmc",
                                            sep = "_")
    mcmc_envir$arguments[[params]] <- do.call(what = mcmc_update_fun,
                                              args = list(mcmc_envir = mcmc_envir))
  } # END for params LOOP

} # END estimate_mupp_params_mcmc_

# INDIVIDUAL UPDATE FUNCTIONS #

# THETAS #
update_mupp_thetas_mcmc <- function(mcmc_envir){

  # all arguments
  arguments    <- mcmc_envir$arguments

  # generate new thetas
  thetas_old   <- thetas_new <- arguments$thetas
  thetas_new[] <- generate_new_params(x  = thetas_old,
                                      sd = arguments$step_size_sd)

  # update thetas in argument list
  arguments    <- modifyList(x   = arguments,
                             val = list(thetas = thetas_new))

  # generate likelihoods
  loglik_old   <- mcmc_envir$loglik
  logliks      <- list(old = loglik_old,
                       new = do.call(what = loglik_mupp_rank_mcmc,
                                     args = arguments))

  # calcualate loglik AND logprior for persons
  logliks_pers <- lapply(X   = logliks,
                         FUN = rowSums)
  priors_pers  <- lapply(X   = list(old = thetas_old,
                                    new = thetas_new),
                         FUN = function(mat){
                           mat[] <- log(d_thetas_prior(mat))
                           rowSums(mat)
                         })

  # compare
  persons_new  <- update_with_metrop(loglik_old = logliks_pers$old,
                                     loglik_new = logliks_pers$new,
                                     priors_old = priors_pers$old,
                                     priors_new = priors_pers$new)

  # updating thetas and items
  thetas_old[persons_new, ]        <- thetas_new[persons_new, ]
  mcmc_envir$loglik[persons_new, ] <- logliks$new[persons_new, ]

  # return
  return(thetas_old)

} # END update_mupp_thetas_mcmc FUNCTION

# ALL PARAMS #
update_mupp_params_mcmc <- function(mcmc_envir,
                                    params_name = "alphas"){

  # all arguments and item list
  arguments     <- mcmc_envir$arguments
  items         <- arguments$items

  # singular params name ??
  d_fun         <- gsub(x = params_name,
                        pattern = "s$",
                        replace = "") %>%
                   paste("d", ., "prior",
                         sep = "_")

  # generate new params
  params_old    <- params_new <- arguments[[params_name]]
  params_new[]  <- generate_new_params(x  = params_old,
                                       sd = arguments$step_size_sd)

  # update thetas in argument list
  arguments     <- modifyList(x   = arguments,
                              val = setNames(list(params_new), params_name))

  # generate likelihoods
  loglik_old    <- mcmc_envir$loglik
  logliks       <- list(old = loglik_old,
                        new = do.call(what = loglik_mupp_rank_mcmc,
                                      args = arguments))

  # calcualate loglik AND logprior for persons
  logliks_item <- lapply(X   = logliks,
                         FUN = colSums)
  priors_item  <- lapply(X   = list(old = params_old,
                                    new = params_new),
                         FUN = function(vec){
                           vec   <- log(do.call(what = d_fun, args = list(vec)))
                           tapply(X     = vec,
                                  INDEX = items[ , 1],
                                  FUN   = sum)
                         })

  # compare
  items_new    <- update_with_metrop(loglik_old = logliks_item$old,
                                     loglik_new = logliks_item$new,
                                     priors_old = priors_item$old,
                                     priors_new = priors_item$new)

  # determine the new statements
  # - items     are 1:n_items
  # - items_new is  T/F according to unique items
  # - items     repeat for items in statements
  state_new    <- items[ , 1] %in% seq_along(items_new)[items_new]

  # updating params and items
  params_old[state_new]           <- params_new[state_new]
  mcmc_envir$loglik[ , items_new] <- logliks$new[ , items_new]

  # return
  return(params_old)

} # END update_mupp_params_mcmc FUNCTION

# ALPHA/DELTA/TAU #
update_mupp_alphas_mcmc <- function(mcmc_envir){
  update_mupp_params_mcmc(mcmc_envir  = mcmc_envir,
                          params_name = "alphas")
} # END update_mupp_alphas_mcmc FUNCTION

update_mupp_deltas_mcmc <- function(mcmc_envir){
  update_mupp_params_mcmc(mcmc_envir  = mcmc_envir,
                          params_name = "deltas")
} # END update_mupp_deltas_mcmc FUNCTION

update_mupp_taus_mcmc <- function(mcmc_envir){
  update_mupp_params_mcmc(mcmc_envir  = mcmc_envir,
                          params_name = "taus")
} # END update_mupp_taus_mcmc FUNCTION

# INITIALIZATION #
initialize_mupp_params_mcmc <- function(resp, items,
                                        delta_sign = NULL){

  # determine number of persons, dimensions, statements
  n_persons    <- nrow(resp)
  n_dims       <- max(unique(items[ , 3]))
  n_statements <- length(unique(items[ , 2]))

  if(is.null(delta_sign)){
    delta_sign <- 1
  } # END if STATEMENT

  # generating parameters

  # - thetas all start at 0 (assume nothing)
  thetas <- matrix(0,
                   nrow = n_persons,
                   ncol = n_dims)

  # - alphas all start at 1, taus all start at -1
  alphas <- rep(+1, n_statements)
  taus   <- rep(-1, n_statements)

  # - deltas start at +1 for positive statements and -1 for negative statements
  deltas <- rep(delta_sign,
                length.out = n_statements)

  # returnining results
  return(list(thetas = thetas,
              alphas = alphas,
              deltas = deltas,
              taus   = taus))

} # END initialize_mupp_params_mcmc FUNCTION

# HELPER FUNCTIONS #

# generate new parameters given old parameters
generate_new_params <- function(x,
                                sd = .1){
  rnorm(n = length(x), mean = x, sd = sd)
} # END generate_new_params FUNCTION

# modified loglikelihood function
loglik_mupp_rank_mcmc <- function(thetas,
                                  alphas,
                                  deltas,
                                  taus,
                                  items,
                                  resp,
                                  ...){

  loglik_mupp_rank_impl(thetas = thetas,
                        params = cbind(alphas, deltas, taus),
                        items  = items,
                        picked_orders = resp)

} # END loglik_mupp_rank_mcmc FUNCTION

# determine whether to update the MCMC algorithm
update_with_metrop <- function(values_old,
                               values_new,
                               loglik_old,
                               loglik_new,
                               priors_old,
                               priors_new){

  # comparison

  # - p is the acceptance probability using MCMC
  p           <- exp(loglik_new - loglik_old + priors_new - priors_old)
  p[is.na(p)] <- 0

  # - u is a random deviate indicating whether to accept
  u           <- runif(length(p))

  accept      <- u < p

  return(accept)

} # END update_with_metrop FUNCTION
