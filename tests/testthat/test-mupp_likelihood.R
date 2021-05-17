context("test-mupp_likelihood.R")

test_that("mupp likelihoods work", {

  set.seed(243534)

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
  calculate_mupp_likelihood <- function(thetas, params, items, resp){

    # pulling out split
    split_vec    <- items[ , 1]

    # splitting params/dimensions/resp
    params_split <- split.data.frame(x = params,
                                     f = split_vec)
    dims_split   <- split(x = items[ , 3],
                          f = split_vec)
    resp_split   <- as.list(as.data.frame(resp))

    # calculating probabilities (assume works, test other places)
    probs <- Map(f      = p_mupp_rank_impl,
                 thetas = list(thetas),
                 params = params_split,
                 dims   = dims_split,
                 picked_order_id = resp_split)

    # calculating likelihood
    loglik <- log(do.call(cbind, probs))

    return(loglik)
  } # END calculate_mupp_likelihood FUNCTION

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

  # a # if everything is the same, should be proportional on all dimensions

  for(n_dims in 2:4){

    # calculating dims and stuff
    dims   <- seq_len(n_dims)
    thetas <- rep(0, n_dims)
    alphas <- rep(1, n_dims)
    deltas <- rep(0, n_dims)
    taus   <- rep(0, n_dims)

    # binding everything together
    thetas <- rbind(thetas)
    params <- cbind(alphas, deltas, taus)
    items  <- cbind(1, dims, dims)

    # running function (in loop)
    for(order in dims){
      expect_equal(loglik_mupp_rank_impl(thetas = rbind(thetas),
                                         params = cbind(alphas,
                                                        deltas,
                                                        taus),
                                         items  = cbind(1, dims, dims),
                                         picked_orders = rbind(order)),
                   matrix(log(1/factorial(n_dims)), nrow = 1, ncol = 1))
    } # END for order LOOP
  } # END for n_dims LOOP

  # b # multiple items (one person)
  for(n_dims in 2:4){
  for(n_items in 2:4){

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
      expect_equal(loglik_mupp_rank_impl(thetas = thetas,
                                         params = params,
                                         items = items,
                                         picked_orders = rbind(r)),
                   calculate_mupp_likelihood(thetas = thetas,
                                             params = params,
                                             items  = items,
                                             resp   = rbind(r)))
    } # END for r LOOP
  }} # END for n_dims/n_items LOOPS

  # c # multiple persons (one item)
  for(n_dims in 2:4){
  for(n_persons in 2:4){

    n_items      <- 1

    # calculating dims and stuff
    dims         <- seq_len(n_dims)
    n_statements <- n_items * n_dims
    thetas       <- matrix(seq(from = 3, to = -3,
                               length.out = n_dims * n_persons),
                           nrow = n_persons,
                           ncol = n_dims)
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
                               n_comb = n_persons)

    for(r in resp){
      expect_equal(loglik_mupp_rank_impl(thetas = thetas,
                                         params = params,
                                         items = items,
                                         picked_orders = cbind(r)),
                   calculate_mupp_likelihood(thetas = thetas,
                                             params = params,
                                             items  = items,
                                             resp   = cbind(r)))
    } # END for r LOOP
  }} # END for n_dims/n_persons LOOPS

  # d # randomize
  n_persons <- 100
  n_items   <- 50

  for(n_dims in 2:4){

    # assumes these functions work (check elsewhere)
    params_sim <- simulate_mupp_params(n_persons = n_persons,
                                       n_items   = n_items,
                                       n_dims    = n_dims)
    resp_sim   <- do.call(simulate_mupp_resp, params_sim)

    # transform parameters
    params     <- transform_params_to_mupp(params_sim, resp_sim)

    # testing
    expect_equal(loglik_mupp_rank_impl(thetas = params$thetas,
                                       params = params$params,
                                       items  = params$items,
                                       picked_orders = params$resp),
                 calculate_mupp_likelihood(thetas = params$thetas,
                                           params = params$params,
                                           items  = params$items,
                                           resp   = params$resp))
  } # END for n_dims LOOP

})
