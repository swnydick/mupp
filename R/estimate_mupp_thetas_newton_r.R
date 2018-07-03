# ESTIMATION: Newton R #

# OVERALL #
estimate_mupp_thetas_newton_r <- function(resp,
                                          params,
                                          items,
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
                          max_iters  = 100)
  control         <- modifyList(control_default,
                                control)

  # estimating for each person
  for(person in seq_len(n_persons)){

    # updating arguments
    args          <- c(list(resp   = resp[person, , drop = FALSE],
                            params = params,
                            items  = items,
                            n_dims = n_dims),
                       control)

    # estimating theta
    out[[person]] <- do.call(what = estimate_mupp_thetas_newton_r_,
                             args = args)
  } # END for person LOOP

  # update thetas to be everything in out bounded together
  thetas <- do.call(what = rbind,
                    args = out)

  return(list(thetas = thetas))

} # END estimate_mupp_thetas_newton_r FUNCTION

# ALL ITERATIONS #

estimate_mupp_thetas_newton_r_ <- function(...,
                                           max_iters = 100,
                                           eps       = 1e-07,
                                           n_dims    = 3){

  # initial thetas (assuming 0 for everything)
  thetas       <- rbind(rep(0, n_dims))
  l            <- loglik_mupp_rank_with_prior1(thetas, ...)

  # starting thetas
  start_thetas <- rep(list(c(-1, 0, 1)), times = n_dims) %>%
                  do.call(what = expand.grid) %>%
                  as.matrix()

  # repeating a bunch of times
  for(r in seq_len(nrow(start_thetas))){

    # number of dimensions and temporary storage for thetas
    thetas_new <- start_thetas[r, ]

    # determining theta (running through N/R many times)
    for(i in seq_len(max_iters)){
      thetas_old <- thetas_new
      thetas_new <- estimate_mupp_thetas_newton_r0(thetas_old, ...)

      if(any(is.na(thetas_new))){
        break
      } else if(sum((thetas_new - thetas_old) ^ 2) < eps){
        break
      } # END ifelse STATEMENT
    } # END for i LOOP

    # go to next iteration if we reached bad
    if(any(is.na(thetas_new))){
      next
    } # END if STATEMENT

    # determine whether to update thetas/loglik
    # - might be local max rather than global max
    # - might be local/global min ... ?
    l_new <- loglik_mupp_rank_with_prior1(thetas_new, ...)

    if(l_new > l){
      l      <- l_new
      thetas <- thetas_new
    } # END if STATEMENT
  }

  return(unname(thetas))

} # END estimate_mupp_thetas_newton_r_

# SINGLE ITERATION #

estimate_mupp_thetas_newton_r0 <- function(thetas, ...){

  # bind thetas
  thetas <- rbind(thetas)

  # estimate gradient/hessian
  g      <- lder1_mupp_rank_with_prior1(thetas, ...)
  H      <- lder2_mupp_rank_with_prior1(thetas, ...)
  delta  <- solve(H, -g)
  alpha  <- 1 / (1 + 10 * sqrt(sum(delta ^ 2)))

  # update newton
  thetas - solve(H, g)

} # END estimate_mupp_thetas_newton_r0 FUNCTION


# UTILITY FUNCTIONS (MOVE) #
convert_to_hessian <- function(h){

  # total number of diagonal elements
  n <- floor(sqrt(2 * length(h)))

  # replace diags/off-diags
  H               <- diag(h[seq_len(n)])
  H[lower.tri(H)] <- h[-seq_len(n)]
  H[upper.tri(H)] <- H[lower.tri(H)]

  return(H)

} # END convert_to_hessian FUNCTION

# FUNCTIONS TO CALCULATE LOGLIK FOR ONE PERSON #
loglik_mupp_rank_with_prior1 <- function(thetas,
                                         resp,
                                         params,
                                         items,
                                         prior_mean,
                                         prior_sd){
  sum(loglik_mupp_rank_impl(thetas = thetas,
                            params = params,
                            items  = items,
                            picked_orders = resp)) +
  sum(log(d_thetas_prior(thetas, prior_mean, prior_sd)))
} # END loglik_mupp_rank_with_prior FUNCTION

lder1_mupp_rank_with_prior1 <- function(thetas,
                                        resp,
                                        params,
                                        items,
                                        prior_mean,
                                        prior_sd){
  c(lder1_mupp_rank_impl(thetas = thetas,
                         params = params,
                         items  = items,
                         picked_orders = resp) +
    (thetas - prior_mean) / (2 * prior_sd))
} # END lder1_mupp_rank_with_prior1 FUNCTION

lder2_mupp_rank_with_prior1 <- function(thetas,
                                        resp,
                                        params,
                                        items,
                                        prior_mean,
                                        prior_sd){
  convert_to_hessian(lder2_mupp_rank_impl(thetas = thetas,
                                          params = params,
                                          items  = items,
                                          picked_orders = resp) +
                      1 / (2 * prior_sd))
} # END lder2_mupp_rank_with_prior1 FUNCTION
