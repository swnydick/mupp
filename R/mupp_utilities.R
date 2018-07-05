# UTILITY FUNCTIONS FOR USE WITH ESTIMATION ALGORITHMS #

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
  -1 * sum(loglik_mupp_rank_impl(thetas = rbind(thetas),
                                 params = params,
                                 items  = items,
                                 picked_orders = resp))
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

