# ESTIMATION: MODIFIED BFGS R #
# (Li and Fukushima, 2001, p. 32, Algorithm 2)
estimate_mupp_thetas_bfgs <- function(...,
                                      max_iters = 100,
                                      eps       = 1e-07,
                                      n_dims    = 3,
                                      n_starts  = 3){

  # initial thetas (for multiple iterations)
  random_thetas <- tryCatch(optimumLHS(n = starts - 1,
                                       k = n_dims),
                            error = function(x){
                              matrix(nrow = 0,
                                     ncol = n_dims)
                            })
  thetas      <- rbind(rep(0, n_dims), 2 * (random_thetas  - .5))

  # initial output objects (if we cannot beat our current)
  thetas_out  <- thetas[1, ]
  loglik_out  <- -loglik_mupp_rank_with_prior1(thetas_out, ...)
  hessian_out <- diag(length(thetas_out))
  iters_out   <- 0

  # running through algorithm multiple times ...
  for(init in seq_len(nrow(thetas))){

    # number of dimensions and temporary storage for thetas
    thetas_new <- thetas[init, ]
    lder_new   <- lder1_mupp_rank_with_prior1(thetas_new, ...)
    B          <- diag(length(thetas_new))

    # determining theta (running through M-BFGS many times)
    for(iter in seq_len(max_iters)){

      thetas_old <- thetas_new
      lder_old   <- lder_new

      # obtain p_k by solving B_k * p_k = -delta f(x)
      p          <- solve(B, -lder_old)

      # obtain alpha by doing argmin f(x + alpha * p)
      alpha      <- wolfe_line_search(fun  = loglik_mupp_rank_with_prior1,
                                      dfun = lder1_mupp_rank_with_prior1,
                                      x    = thetas_old,
                                      p    = p,
                                      ...)

      # determine s in algorithm, and update thetas and lder
      s          <- alpha * p
      thetas_new <- thetas_old + s
      lder_new   <- lder1_mupp_rank_with_prior1(thetas_new, ...)

      # update other parameters
      gamma      <- lder_new - lder_old
      eta        <- 1 + max(-sum(gamma * s) / sum(s^2), 0)
      y          <- gamma + eta * sqrt(sum(p^2)) * s

      # update B
      B          <- B - (B %*% s %*% t(s) %*% B) / sum(s * B %*% s) + (y %*% t(y)) / sum(s * y)

      # check convergence and return
      if(sum(lder_new^2) < eps){
        break
      } # END if STATEMENT

    } # END for iter LOOP

    # updating log-likelihood
    loglik_new <- -loglik_mupp_rank_with_prior1(thetas_new, ...)

    if(loglik_new > loglik_out){
      thetas_out  <- thetas_new
      loglik_out  <- loglik_new
      hessian_out <- B
      iters_out   <- iter
    } # END if STATEMENT
  } # END for init LOOP

  return(list(thetas  = thetas_out,
              hessian = hessian_out,
              loglik  = loglik_out,
              iters   = iters_out))
} # END estimate_mupp_thetas_bfgs FUNCTION
