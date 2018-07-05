# ESTIMATION: MODIFIED BFGS R #
# (Li and Fukushima, 2001, p. 32, Algorithm 2)
estimate_mupp_thetas_bfgs <- function(...,
                                      max_iters = 100,
                                      eps       = 1e-07,
                                      n_dims    = 3){

  # initial thetas (assuming 0 for everything)
  thetas       <- rep(0, n_dims)

  # number of dimensions and temporary storage for thetas
  thetas_new   <- thetas
  lder_new     <- lder1_mupp_rank_with_prior1(thetas_new, ...)
  B            <- diag(length(thetas_new))

  # determining theta (running through M-BFGS many times)
  for(i in seq_len(max_iters)){

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
  }

  return(list(thetas  = thetas_new,
              hessian = B,
              iters   = i))
} # END estimate_mupp_thetas_bfgs FUNCTION
