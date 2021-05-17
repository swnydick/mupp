context("test-mupp_probs_3D.R")

test_that("mupp probabilities work: three dimensions", {

  set.seed(9872345)

  # 0 # Misc Objects
  dims      <- 1:3
  n_dims    <- length(dims)

  all_perms <- list(c(1, 2, 3),
                    c(1, 3, 2),
                    c(2, 1, 3),
                    c(2, 3, 1),
                    c(3, 1, 2),
                    c(3, 2, 1))

  # 0 # R functions to calculate probs

  # - GGUM
  p_ggum <- function(theta, alpha, delta, tau){
    e1 <- exp(alpha * ((theta - delta) - tau))
    e2 <- exp(alpha * (2 * (theta - delta) - tau))
    e3 <- exp(alpha * (3 * (theta - delta)))
    (e1 + e2) / (1 + e1 + e2 + e3)
  } # END p_ggum FUNCTION

  # - MUPP (for 3D only) (ONE permutation, intentionally manual and hard coded)
  p_mupp0 <- function(thetas, params, perm){

    # indices
    i <- perm[1]
    j <- perm[2]
    k <- perm[3]

    # all probabilities
    p_0 <- cbind(p_ggum(thetas[ , i], params[i, 1], params[i, 2], params[i, 3]),
                 p_ggum(thetas[ , j], params[j, 1], params[j, 2], params[j, 3]),
                 p_ggum(thetas[ , k], params[k, 1], params[k, 2], params[k, 3]))
    q_0 <- 1 - p_0

    # combined probabilities
    num1 <- p_0[ , 1] * q_0[ , 2] * q_0[ , 3]
    den1 <- p_0[ , 2] * q_0[ , 1] * q_0[ , 3] + p_0[ , 3] * q_0[ , 1] * q_0[ , 2] + num1
    p1   <- num1 / den1

    num2 <- p_0[ , 2] * q_0[ , 3]
    den2 <- p_0[ , 3] * q_0[ , 2] + num2
    p2   <- num2 / den2

    p1 * p2
  } # END p_mupp0 (3D) FUNCTION

  # - MUPP (for 3D only) (ALL permutations)
  p_mupp <- function(thetas, params){
    P <- lapply(all_perms,
                FUN    = p_mupp0,
                thetas = thetas,
                params = params)
    P <- "dimnames<-"(do.call(cbind, P), NULL)

    return(P)
  } # END p_mupp (3D) FUNCTION

  # a # if everything is the same, should be same on all dimensions
  expect_equal(p_mupp_rank_impl(thetas = rbind(c(0, 0, 0)),
                                params = cbind(c(1, 1, 1),
                                               c(0, 0, 0),
                                               c(0, 0, 0)),
                                dims   = dims),
               matrix(1/factorial(n_dims), nrow = 1, ncol = factorial(n_dims)))

  # b # change theta but keep items the same
  thetas <- rbind(c(0, -1, -2))
  alphas <- c(1, 1, 1)
  deltas <- c(0, 0, 0)
  taus   <- c(0, 0, 0)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # c # change alpha but keep items the same
  thetas <- rbind(c(0, 0, 0))
  alphas <- c(1, 2, 2.5)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # d # change delta but keep items the same
  alphas <- c(1, 1, 1)
  deltas <- c(0, -1, 1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # d # change taus but keep items the same
  deltas <- c(0, 0, 0)
  taus   <- c(0, -1, 1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # d # change everything
  thetas <- rbind(c(0, 1, 2))
  alphas <- c(1, 2, 2.5)
  deltas <- c(0, -1, 1)
  taus   <- c(0, -1, 1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # e # multiple people RANDOM, so different each time :)
  n_thetas <- 9
  thetas   <- matrix(r_thetas_prior(n_thetas * n_dims),
                     nrow = n_thetas,
                     ncol = n_dims)
  alphas   <- r_alpha_prior(n_dims)
  deltas   <- r_delta_prior(n_dims)
  taus     <- r_tau_prior(n_dims)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # f # selecting one dimension

  # - first dimension - same for everyone
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = 1),
               p_mupp(thetas, params)[ , 1, drop = FALSE])

  # - second dimension - same for everyone
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = 2),
               p_mupp(thetas, params)[ , 2, drop = FALSE])

  # - third dimension - same for everyone
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = 3),
               p_mupp(thetas, params)[ , 3, drop = FALSE])

  # - first/second dimension altering
  ids    <- rep(dims, length.out = n_thetas)
  id_mat <- cbind(seq_len(n_thetas), ids)
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = ids),
               cbind(p_mupp(thetas, params)[id_mat]))

  # - first/second dimension strange
  ids    <- sample(dims, size = n_thetas, replace = TRUE)
  id_mat <- cbind(seq_len(n_thetas), ids)
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = ids),
               cbind(p_mupp(thetas, params)[id_mat]))
})
