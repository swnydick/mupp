context("test-mupp_d2probs_2D.R")

test_that("mupp probability second derivatives work: two dimensions", {

  set.seed(23453)

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

  # 0 # R functions to calculate derivatives
  pder1_ggum <- function(theta, alpha, delta, tau){
    e1 <- exp(alpha * ((theta - delta) - tau))
    e2 <- exp(alpha * (2 * (theta - delta) - tau))
    e3 <- exp(alpha * (3 * (theta - delta)))
    alpha * (e1 * (1 - 2 * e3) + e2 * (2 - e3)) / (1 + e1 + e2 + e3)^2
  }

  # 0 # R functions to calculate second derivatives
  pder2_ggum <- function(theta, alpha, delta, tau){
    e1 <- exp(alpha * ((theta - delta) - tau))
    e2 <- exp(alpha * (2 * (theta - delta) - tau))
    e3 <- exp(alpha * (3 * (theta - delta)))

    A  <- alpha   * (e1 * (1 - 2 * e3) + e2 * (2 - e3))
    dA <- alpha^2 * (e1 * (1 - 8 * e3) + e2 * (4 - 5 * e3))
    B  <- (1 + e1 + e2 + e3) ^ 2
    dB <- 2 * alpha * (1 + e1 + e2 + e3) * (e1 + 2 * e2 + 3 * e3)
    (B * dA - A * dB) / (B ^ 2)
  } # END pder2_ggum FUNCTION

  # - MUPP (for 2D only) (intentionally manual and hard coded and NOT simplified)
  #                      (same as in paper, even though paper could be simplified more)
  pder2_mupp <- function(thetas, params){
    A <- p_ggum(thetas[ , 1], params[1, 1], params[1, 2], params[1, 3])
    B <- 1 - p_ggum(thetas[ , 2], params[2, 1], params[2, 2], params[2, 3])
    C <- 1 - A
    D <- 1 - B

    Ap <-  pder1_ggum(thetas[ , 1], params[1, 1], params[1, 2], params[1, 3])
    Bp <- -pder1_ggum(thetas[ , 2], params[2, 1], params[2, 2], params[2, 3])

    A2p <-  pder2_ggum(thetas[ , 1], params[1, 1], params[1, 2], params[1, 3])
    B2p <- -pder2_ggum(thetas[ , 2], params[2, 1], params[2, 2], params[2, 3])

    d1  <- B * (1 - B) * (A2p *  (A * B + C * D) - 2 * Ap^2 * (2 * B - 1)) / (A * B + C * D)^3
    d2  <- A * (1 - A) * (B2p *  (A * B + C * D) - 2 * Bp^2 * (2 * A - 1)) / (A * B + C * D)^3
    d12 <- Ap * Bp * (1 - A - B) / (A * B + C * D)^3

    dP <- "dimnames<-"(cbind(d1, d2, d12), NULL)
    dP
  } # END pder2_mupp FUNCTION

  # 0 # R function to put second derivatives in correct order
  pder2_mupp_all <- function(thetas, params){
    d1 <- pder2_mupp(thetas, params)
    d2 <- pder2_mupp(thetas[ , 2:1, drop = FALSE],
                     params[2:1, , drop = FALSE])[ , c(2:1, 3), drop = FALSE]
    dP <- "dimnames<-"(list(d1, d2), NULL)
    dP
  }

  # 0 # Misc Functions
  convert_hessian <- function(mat){
    rbind(c(diag(mat), mat[lower.tri(mat, diag = FALSE)]))
  }

  # b # change theta but keep items the same
  thetas <- rbind(c(0, -1))
  alphas <- c(1, 1)
  deltas <- c(0, 0)
  taus   <- c(0, 0)
  params <- cbind(alphas, deltas, taus)

  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims),
               pder2_mupp_all(thetas, params))

  # c # change alpha but keep items the same
  thetas <- rbind(c(0, 0))
  alphas <- c(1, 2)
  params <- cbind(alphas, deltas, taus)

  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims),
               pder2_mupp_all(thetas, params))

  # d # change delta but keep items the same
  alphas <- c(1, 1)
  deltas <- c(0, -1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims),
               pder2_mupp_all(thetas, params))

  # d # change taus but keep items the same
  deltas <- c(0, 0)
  taus   <- c(0, -1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims),
               pder2_mupp_all(thetas, params))

  # d # change everything
  thetas <- rbind(c(0, 1))
  alphas <- c(1, 2)
  deltas <- c(0, -1)
  taus   <- c(0, -1)
  params <- cbind(alphas, deltas, taus)

  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims),
               pder2_mupp_all(thetas, params))

  # e # multiple people RANDOM, so different each time :)
  n_thetas <- 5
  thetas   <- matrix(r_thetas_prior(n_thetas * n_dims),
                     nrow = n_thetas,
                     ncol = n_dims)
  alphas   <- r_alpha_prior(n_dims)
  deltas   <- r_delta_prior(n_dims)
  taus     <- r_tau_prior(n_dims)

  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims),
               pder2_mupp_all(thetas, params))

  # f # selecting one dimension

  # - first dimension - same for everyone
  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims,
                                    picked_order_id = 1),
               pder2_mupp_all(thetas, params)[1])

  # - second dimension - same for everyone
  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims,
                                    picked_order_id = 2),
               pder2_mupp_all(thetas, params)[2])

  # - first/second dimension altering
  ids    <- rep(dims, length.out = n_thetas)
  pder2  <- pder2_mupp_all(thetas, params)
  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims,
                                    picked_order_id = ids),
               list(do.call(what = rbind,
                            args = lapply(seq_along(ids),
                                          FUN = function(i)
                                            pder2[[ids[i]]][i, , drop = FALSE]))))

  # - first/second dimension random order
  ids    <- sample(dims, size = n_thetas, replace = TRUE)
  pder2  <- pder2_mupp_all(thetas, params)
  expect_equal(pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = dims,
                                    picked_order_id = ids),
               list(do.call(what = rbind,
                            args = lapply(seq_along(ids),
                                          FUN = function(i)
                                            pder2[[ids[i]]][i, , drop = FALSE]))))

  # g # using grad from numDeriv

  # - first dimension, first item picked
  thetas    <- thetas[1, , drop = FALSE]
  pder2     <- numDeriv::hessian(func = function(x){
                                   p_mupp_rank_impl(thetas = x,
                                                    params = params,
                                                    dims   = c(1, 1),
                                                    picked_order_id = 1)
                                 },
                                 x = thetas)
  pder2_a   <- pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = c(1, 1),
                                    picked_order_id = 1)[[1]]
  expect_equivalent(all(abs(convert_hessian(pder2) - pder2_a) < 1e-09),
                    TRUE)

  # - first dimension, second item picked
  pder2     <- numDeriv::hessian(func = function(x){
                                   p_mupp_rank_impl(thetas = x,
                                                    params = params,
                                                    dims   = c(1, 1),
                                                    picked_order_id = 2)
                                 },
                                 x = thetas)
  pder2_a   <- pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = c(1, 1),
                                    picked_order_id = 2)[[1]]
  expect_equivalent(all(abs(convert_hessian(pder2) - pder2_a) < 1e-09),
                    TRUE)

  # - second dimension, first item picked
  pder2     <- numDeriv::hessian(func = function(x){
                                   p_mupp_rank_impl(thetas = x,
                                                    params = params,
                                                    dims   = c(2, 2),
                                                    picked_order_id = 1)
                                 },
                                 x = thetas)
  pder2_a   <- pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = c(2, 2),
                                    picked_order_id = 1)[[1]]
  expect_equivalent(all(abs(convert_hessian(pder2) - pder2_a) < 1e-09),
                    TRUE)

  # - second dimension, second item picked
  pder2     <- numDeriv::hessian(func = function(x){
                                   p_mupp_rank_impl(thetas = x,
                                                    params = params,
                                                    dims   = c(2, 2),
                                                    picked_order_id = 2)
                                 },
                                 x = thetas)
  pder2_a   <- pder2_mupp_rank_impl(thetas = thetas,
                                    params = params,
                                    dims   = c(2, 2),
                                    picked_order_id = 2)[[1]]
  expect_equivalent(all(abs(convert_hessian(pder2) - pder2_a) < 1e-09),
                    TRUE)

})
