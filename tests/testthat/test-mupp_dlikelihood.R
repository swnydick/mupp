context("test-mupp_dlikelihood.R")

test_that("mupp likelihood derivatives work", {

  set.seed(872934)

  # 0 # R functions to calculate likelihoods

  # turn parameters into useful objects for MUPP
  transform_params_to_mupp <- function(params, resp){

    # pulling out objects (thetas, items, params)
    thetas <- as.matrix(tidyr::spread(params$persons,
                                      key   = "dim",
                                      value = "theta")[-1])
    items  <- as.matrix(params$items[1:3])
    params <- as.matrix(params$items[4:6])

    # pulling out objects (resp)
    resp   <- as.matrix(tidyr::spread(resp$resp[c("person", "item", "resp")],
                                      key   = "item",
                                      value = "resp")[-1])

    return(list(thetas = thetas, params = params, items = items, resp = resp))
  } # END transform_params_to_mupp FUNCTION

  # calculating MUPP likelihood (NOT using C++ function)
  estimate_mupp_dlikelihood <- function(thetas, params, items, resp){

    g <- numDeriv::grad(func = function(x){
                          sum(loglik_mupp_rank_impl(rbind(x),
                                                    params,
                                                    items,
                                                    resp))
                        }, x = thetas)
    h <- numDeriv::hessian(func = function(x){
                             sum(loglik_mupp_rank_impl(rbind(x),
                                                       params,
                                                       items,
                                                       resp))
                           }, x = thetas)

    list(gradient = g,
         hessian  = rbind(c(diag(h), h[lower.tri(h, diag = FALSE)])))
  } # END calculate_mupp_dlikelihood FUNCTION

  # creating list of possible responses (nPk permutations)
  create_resp_list <- function(n_dims,
                               n_comb){

    # possible selections for choose matrix (allows permutations!)
    select_opts <- rep(x     = seq_len(n_dims),
                       times = n_comb)

    # determine possible nPk permutations
    mat <- t(combn(x = select_opts,
                   m = n_comb,
                   simplify = TRUE))
    df  <- unique(as.data.frame(mat))

    # arrange in useful order
    df  <- df[do.call(base::order, args = df), ]

    # turn into list and return
    as.list(as.data.frame(t(df)))

  } # END create_resp_list FUNCTION

  # b # one item (one person, two dimensions)

  n_items      <- 1
  n_dims       <- 2

  # calculating dims and stuff
  dims         <- seq_len(n_dims)
  n_statements <- n_items * n_dims
  thetas       <- seq(from = 3, to = -3,
                      length.out = n_dims)
  alphas       <- rep(x = 1:2,
                      length.out = n_statements)
  deltas       <- seq(from = -2, to = 2,
                      length.out = n_statements)
  taus         <- rep(x = c(-1, 0, 1),
                      length.out = n_statements)

  # binding everything together
  thetas <- rbind(thetas)
  params <- cbind(alphas, deltas, taus)
  items  <- cbind(item      = rep(x    = seq_len(n_statements / n_dims),
                                  each = n_dims),
                  statement = seq_len(n_statements),
                  dim       = rep(dims,
                                  times = n_items))
  # creating response options
  resp   <- create_resp_list(n_dims = n_dims,
                             n_comb = n_items)

  for(r in resp){

    # estimating likelihood
    lders   <- estimate_mupp_dlikelihood(thetas = thetas,
                                         params = params,
                                         items  = items,
                                         resp   = rbind(r))
    lder1   <- lders$gradient
    lder2   <- lders$hessian

    # actual likelihood
    lder1_a <- lder1_mupp_rank_impl(thetas, params, items, rbind(r))
    lder2_a <- lder2_mupp_rank_impl(thetas, params, items, rbind(r))

    # comparing
    expect_equivalent(all(abs(lder1 - lder1_a) < 1e-07),
                      TRUE)
    expect_equivalent(all(abs(lder2 - lder2_a) < 1e-07),
                      TRUE)
  } # END for r LOOP

  # c # one item (one person, one dimension)
  for(dim in 1:2){

    n_items      <- 1
    n_dims       <- 2

    # calculating dims and stuff
    dims         <- rep(dim, n_dims)
    n_statements <- n_items * n_dims
    thetas       <- seq(from = 3, to = -3,
                        length.out = n_dims)
    alphas       <- rep(x = 1:2,
                        length.out = n_statements)
    deltas       <- seq(from = -2, to = 2,
                        length.out = n_statements)
    taus         <- rep(x = c(-1, 0, 1),
                        length.out = n_statements)

    # binding everything together
    thetas <- rbind(thetas)
    params <- cbind(alphas, deltas, taus)
    items  <- cbind(item      = rep(x    = seq_len(n_statements / n_dims),
                                    each = n_dims),
                    statement = seq_len(n_statements),
                    dim       = rep(dims,
                                    times = n_items))
    # creating response options
    resp   <- create_resp_list(n_dims = n_dims,
                               n_comb = n_items)

    for(r in resp){

      # estimating likelihood
      lders   <- estimate_mupp_dlikelihood(thetas = thetas,
                                           params = params,
                                           items  = items,
                                           resp   = rbind(r))
      lder1   <- lders$gradient
      lder2   <- lders$hessian

      # actual likelihood
      lder1_a <- lder1_mupp_rank_impl(thetas, params, items, rbind(r))
      lder2_a <- lder2_mupp_rank_impl(thetas, params, items, rbind(r))

      # comparing
      expect_equivalent(all(abs(lder1 - lder1_a) < 1e-07),
                        TRUE)
      expect_equivalent(all(abs(lder2 - lder2_a) < 1e-07),
                        TRUE)
    } # END for r LOOP
  } # END for dim LOOP


  # d # multiple items (one person, two dimensions)
  for(n_items in 2:4){
    n_dims       <- 2

    # calculating dims and stuff
    dims         <- seq_len(n_dims)
    n_statements <- n_items * n_dims
    thetas       <- seq(from = 3, to = -3,
                        length.out = n_dims)
    alphas       <- rep(x = 1:2,
                        length.out = n_statements)
    deltas       <- seq(from = -2, to = 2,
                        length.out = n_statements)
    taus         <- rep(x = c(-1, 0, 1),
                        length.out = n_statements)

    # binding everything together
    thetas <- rbind(thetas)
    params <- cbind(alphas, deltas, taus)
    items  <- cbind(item      = rep(x    = seq_len(n_statements / n_dims),
                                    each = n_dims),
                    statement = seq_len(n_statements),
                    dim       = rep(c(1, 1),
                                    times = n_items))
    # creating response options
    resp   <- create_resp_list(n_dims = n_dims,
                               n_comb = n_items)

    for(r in resp){

      # estimating likelihood
      lders   <- estimate_mupp_dlikelihood(thetas = thetas,
                                           params = params,
                                           items  = items,
                                           resp   = rbind(r))
      lder1   <- lders$gradient
      lder2   <- lders$hessian

      # actual likelihood
      lder1_a <- lder1_mupp_rank_impl(thetas, params, items, rbind(r))
      lder2_a <- lder2_mupp_rank_impl(thetas, params, items, rbind(r))

      # comparing
      expect_equivalent(all(abs(lder1 - lder1_a) < 1e-07),
                        TRUE)
      expect_equivalent(all(abs(lder2 - lder2_a) < 1e-07),
                        TRUE)
    } # END for r LOOP
  } # END for n_items LOOP

  # e # multiple items (one person, one dimension)
  for(n_items in 2:4){

    n_dims       <- 2

    # calculating dims and stuff
    dims         <- list(c(1, 1), c(2, 2))[sample(1:2, size = n_items, replace = TRUE)]
    n_statements <- n_items * n_dims
    thetas       <- seq(from = 3, to = -3,
                        length.out = n_dims)
    alphas       <- rep(x = 1:2,
                        length.out = n_statements)
    deltas       <- seq(from = -2, to = 2,
                        length.out = n_statements)
    taus         <- rep(x = c(-1, 0, 1),
                        length.out = n_statements)

    # binding everything together
    thetas <- rbind(thetas)
    params <- cbind(alphas, deltas, taus)
    items  <- cbind(item      = rep(x    = seq_len(n_statements / n_dims),
                                    each = n_dims),
                    statement = seq_len(n_statements),
                    dim       = unlist(dims))
    # creating response options
    resp   <- create_resp_list(n_dims = n_dims,
                               n_comb = n_items)

    for(r in resp){

      # estimating likelihood
      lders   <- estimate_mupp_dlikelihood(thetas = thetas,
                                           params = params,
                                           items  = items,
                                           resp   = rbind(r))
      lder1   <- lders$gradient
      lder2   <- lders$hessian

      # actual likelihood
      lder1_a <- lder1_mupp_rank_impl(thetas, params, items, rbind(r))
      lder2_a <- lder2_mupp_rank_impl(thetas, params, items, rbind(r))

      # comparing
      expect_equivalent(all(abs(lder1 - lder1_a) < 1e-05),
                        TRUE)
      expect_equivalent(all(abs(lder2 - lder2_a) < 1e-05),
                        TRUE)
    } # END for r LOOP
  } # END for n_items LOOP

  # f # completely random
  n_items       <- 50
  n_dims        <- 5
  max_item_dims <- 2

  # assume these functions work (check elsewhere)
  params_sim   <- simulate_mupp_params(n_persons     = 1,
                                       n_items       = n_items,
                                       n_dims        = n_dims,
                                       max_item_dims = 2,
                                       unidim_items  = TRUE)
  resp_sim     <- do.call(simulate_mupp_resp, params_sim)

  # transform parameters
  params_all   <- transform_params_to_mupp(params_sim, resp_sim)
  thetas       <- params_all$thetas
  params       <- params_all$params
  resp         <- params_all$resp
  items        <- params_all$items

  # estimating likelihood
  lders   <- estimate_mupp_dlikelihood(thetas = thetas,
                                       params = params,
                                       items  = items,
                                       resp   = resp)
  lder1   <- lders$gradient
  lder2   <- lders$hessian

  # actual likelihood
  lder1_a <- lder1_mupp_rank_impl(thetas, params, items, resp)
  lder2_a <- lder2_mupp_rank_impl(thetas, params, items, resp)

  # comparing
  expect_equivalent(all(abs(lder1 - lder1_a) < 1e-05),
                    TRUE)
  expect_equivalent(all(abs(lder2 - lder2_a) < 1e-05),
                    TRUE)

})
