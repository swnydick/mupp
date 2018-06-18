context("test-mupp_probs_2D.R")

test_that("mupp probabilities work: two dimensions", {

  # 0 # Misc Objects
  dims   <- 1:2
  n_dims <- length(dims)

  # 0 # R functions to calculate probs

  # - GGUM
  p_ggum <- function(theta, alpha, delta, tau){
    e1 <- exp(alpha * ((theta - delta) - tau))
    e2 <- exp(alpha * (2 * (theta - delta) - tau))
    e3 <- exp(alpha * (3 * (theta - delta)))
    (e1 + e2) / (1 + e1 + e2 + e3)
  } # END p_ggum FUNCTION

  # - MUPP (for 2D only) (intentionally manual and hard coded)
  p_mupp <- function(thetas, params){
    p_s <- p_ggum(thetas[ , 1], params[1, 1], params[1, 2], params[1, 3])
    p_t <- p_ggum(thetas[ , 2], params[2, 1], params[2, 2], params[2, 3])
    p   <- p_s * (1 - p_t) / (p_s * (1 - p_t) + p_t * (1 - p_s))
    P   <- "dimnames<-"(cbind(p, 1 - p), NULL)
    P
  } # END p_mupp (2D) FUNCTION

  # a # if everything is the same, should be same prob for both orders
  expect_equal(p_mupp_rank_impl(thetas = rbind(c(0, 0)),
                                params = cbind(c(1, 1),
                                               c(0, 0),
                                               c(0, 0)),
                                dims   = dims),
               matrix(1/factorial(n_dims), nrow = 1, ncol = factorial(n_dims)))

  # b # change theta but keep items the same
  thetas <- rbind(c(0, -1))
  alphas <- c(1, 1)
  deltas <- c(0, 0)
  taus   <- c(0, 0)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # c # change alpha but keep items the same
  thetas <- rbind(c(0, 0))
  alphas <- c(1, 2)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # d # change delta but keep items the same
  alphas <- c(1, 1)
  deltas <- c(0, -1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # d # change taus but keep items the same
  deltas <- c(0, 0)
  taus   <- c(0, -1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # d # change everything
  thetas <- rbind(c(0, 1))
  alphas <- c(1, 2)
  deltas <- c(0, -1)
  taus   <- c(0, -1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims),
               p_mupp(thetas, params))

  # e # multiple people RANDOM, so different each time :)
  n_thetas <- 5
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

  # - first/second dimension altering
  ids    <- rep(dims, length.out = n_thetas)
  id_mat <- cbind(seq_len(n_thetas), ids)
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = ids),
               cbind(p_mupp(thetas, params)[id_mat]))

  # - first/second dimension random order
  ids    <- sample(dims, size = n_thetas, replace = TRUE)
  id_mat <- cbind(seq_len(n_thetas), ids)
  expect_equal(p_mupp_rank_impl(thetas = thetas,
                                params = params,
                                dims   = dims,
                                picked_order_id = ids),
               cbind(p_mupp(thetas, params)[id_mat]))
})
