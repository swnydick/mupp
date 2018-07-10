# ESTIMATION: MCMC #

# OVERALL #
estimate_mupp_params_mcmc <- function(resp,
                                      items,
                                      control        = list(),
                                      initial_params = list(),
                                      fixed_params   = NULL){

  # converting items to a matrix (so that it works in the C++ algorithm)
  items <- as.matrix(items)

  # update control
  control_default <- list(n_iters      = 10000,
                          n_burnin     = 1000,
                          step_size_sd = .1)
  control         <- modifyList(control_default,
                                control)

  # fix initial_params/fixed_params
  initial_params %<>% setNames(fix_param_names(names(.)))
  fixed_params   %<>% fix_param_names(.)

  # INITIALIZE THETA/STATEMENT VALUES #
  current_params   <- initialize_mupp_params_mcmc(resp  = resp,
                                                  items = items,
                                                  initial_params = initial_params)

  # environment to store mcmc arguments and current loglikelihood
  mcmc_temp_envir  <- new.env()

  # adding argument to the environment
  assign(x     = "arguments",
         value = modifyList(x   = current_params,
                            val = list(items        = items,
                                       resp         = resp,
                                       step_size_sd = control$step_size_sd,
                                       mcmc_envir   = mcmc_temp_envir)),
         envir = mcmc_temp_envir)

  # adding initial loglikelihood to the environment
  assign(x     = "loglik",
         value = do.call(what = loglik_mupp_rank_mcmc,
                         args = mcmc_temp_envir$arguments),
         envir = mcmc_temp_envir)


  # setup output stuff
  n_iters        <- control$n_iters
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
    estimate_mupp_params_mcmc_(mcmc_envir   = mcmc_temp_envir,
                               fixed_params = fixed_params)

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
  n_burnin        <- control$n_burnin
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
                           ss      <- Reduce("+", lapply(sims, "^", 2))
                           len     <- length(sims)
                           means[] <- sqrt(pmax(0, (ss / len - means^2) * (len - 1) / len))
                           return(means)
                         })

  list(all   = all_mcmc_steps,
       means = mean_sims,
       sds   = sd_sims)
} # END estimate_mupp_params_mcmc FUNCTION

# SINGLE ITERATION #

estimate_mupp_params_mcmc_ <- function(mcmc_envir,
                                       fixed_params = NULL){

  # for each param (if NOT FIXED):
  #  - find the update function
  #  - update the param (using arguments) and assign back to the environment
  for(params in c("thetas", "alpha", "delta", "tau")){
    if(params %ni% fixed_params){
      mcmc_update_fun                <- paste("update_mupp", params, "mcmc",
                                              sep = "_")
      mcmc_envir$arguments[[params]] <- do.call(what = mcmc_update_fun,
                                                args = list(mcmc_envir = mcmc_envir))
    } # END if STATEMENT
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
  logliks_pers <- lapply(X     = logliks,
                         FUN   = rowSums,
                         na.rm = TRUE)
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
                                    params_name = "alpha"){

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
  logliks_item <- lapply(X     = logliks,
                         FUN   = colSums,
                         na.rm = TRUE)
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
update_mupp_alpha_mcmc <- function(mcmc_envir){
  update_mupp_params_mcmc(mcmc_envir  = mcmc_envir,
                          params_name = "alpha")
} # END update_mupp_alpha_mcmc FUNCTION

update_mupp_delta_mcmc <- function(mcmc_envir){
  update_mupp_params_mcmc(mcmc_envir  = mcmc_envir,
                          params_name = "delta")
} # END update_mupp_delta_mcmc FUNCTION

update_mupp_tau_mcmc <- function(mcmc_envir){
  update_mupp_params_mcmc(mcmc_envir  = mcmc_envir,
                          params_name = "tau")
} # END update_mupp_tau_mcmc FUNCTION

# INITIALIZATION AND NAMING #

# naming parameters (correctly)
fix_param_names <- function(x){

  # determine old and new names for replacing
  old_names <- c("theta", "alphas", "deltas", "taus")
  new_names <- c("thetas", "alpha", "delta", "tau")

  # fixing names
  if(any(flag <- x %in% old_names)){
    x[flag] %<>% {new_names[match(., old_names)]}
  } # END if STATEMENT

  return(x)

} # END fix_param_names FUNCTION


# initialization
initialize_mupp_params_mcmc <- function(resp, items,
                                        initial_params = NULL){

  # determine number of persons, dimensions, statements
  n_persons    <- nrow(resp)
  n_dims       <- max(unique(items[ , 3]))
  n_statements <- length(unique(items[ , 2]))

  # default parameters #

  # - thetas all start at 0 (assume nothing)
  thetas <- matrix(0,
                   nrow = n_persons,
                   ncol = n_dims)

  # - alphas all start at 1, taus all start at -1
  alpha  <- rep(+1, n_statements)
  delta  <- rep(0,  n_statements)
  tau    <- rep(-1, n_statements)

  # updating parameters (if initial params exists)
  if(!is.null(initial_params$alpha)){
    alpha[] <- initial_params$alpha
  } # END if STATEMENT
  if(!is.null(initial_params$delta)){
    delta[] <- initial_params$delta
  } # END if STATEMENT
  if(!is.null(initial_params$tau)){
    tau[]   <- initial_params$tau
  } # END if STATEMENT
  if(!is.null(initial_params$thetas)){
    thetas[] <- initial_params$thetas
  } # END if STATEMENT

  # returnining results
  return(list(thetas = thetas,
              alpha = alpha,
              delta = delta,
              tau   = tau))

} # END initialize_mupp_params_mcmc FUNCTION

# HELPER FUNCTIONS #

# generate new parameters given old parameters
generate_new_params <- function(x,
                                sd = .1){
  rnorm(n = length(x), mean = x, sd = sd)
} # END generate_new_params FUNCTION

# modified loglikelihood function
loglik_mupp_rank_mcmc <- function(thetas,
                                  alpha,
                                  delta,
                                  tau,
                                  items,
                                  resp,
                                  ...){

  loglik_mupp_rank_impl(thetas = thetas,
                        params = cbind(alpha, delta, tau),
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
