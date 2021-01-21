context("test-simulate_mupp.R")

test_that("mupp response conversion works", {

  # a # guaranteed responses

  # one person
  probs <- list(c(1, 0, 0, 0),
                c(0, 1, 0, 0),
                c(0, 0, 1, 0),
                c(0, 0, 0, 1))

  for(i in seq_along(probs)){
    prob <- rbind(probs[[i]])
    resp <- apply(prob,
                  MARGIN = 1,
                  FUN    = function(x)
                    which(x == 1))

    expect_equal(simulate_mupp_resp1(prob),
                 resp)
  } # END for i LOOP

  # multiple people
  probs <- list(rbind(c(1, 0, 0, 0),
                      c(0, 1, 0, 0),
                      c(0, 0, 1, 0),
                      c(0, 0, 0, 1)),
                rbind(c(1, 0, 0, 0),
                      c(1, 0, 0, 0),
                      c(1, 0, 0, 0),
                      c(1, 0, 0, 0)),
                rbind(c(0, 1, 0, 0),
                      c(1, 0, 0, 0),
                      c(0, 1, 0, 0),
                      c(1, 0, 0, 0)))

  for(i in seq_along(probs)){
    prob <- rbind(probs[[i]])
    resp <- apply(prob,
                  MARGIN = 1,
                  FUN    = function(x)
                    which(x == 1))

    expect_equal(simulate_mupp_resp1(prob),
                 resp)
  } # END for i LOOP

  # b # probabilistic response (use set.seed to guarantee response)

  # one person (two dimensions)

  # probs
  probs <- list(c(.5, .5),
                c(.3, .7),
                c(.7, .3),
                c(.1, .9),
                c(.9, .1))

  # seed
  seed  <- 768034

  # running loop
  for(i in seq_along(probs)){

    # pull out probabilities
    prob <- rbind(probs[[i]])

    # determine response
    set.seed(seed)
    u    <- runif(nrow(prob))
    resp <- (u >= prob[ , 1]) + 1

    # compare
    set.seed(seed)
    expect_equal(simulate_mupp_resp1(prob),
                 resp)
  } # END for i LOOP

  # one person (four dimensions)

  # probs (resp is after generating u from seed)
  probs <- list(c(.25, .25, .25, .25),
                c(.2, .3, .3, .2),
                c(.1, .9, 0, 0),
                c(0, 0, .1, .9))
  resp  <- c(1, 1, 2, 4)

  # seed
  seed  <- 423422

  # running loop
  for(i in seq_along(probs)){

    # pull out probabilities
    prob <- rbind(probs[[i]])

    # compare
    set.seed(seed)
    expect_equal(simulate_mupp_resp1(prob),
                 resp[i])
  } # END for i LOOP


  # multiple people (two dimensions)

  # probs
  probs <- list(rbind(c(.5, .5),
                      c(.4, .6),
                      c(.2, .8)),
                rbind(c(.3, .7),
                      c(.7, .3)))

  # seed
  seed  <- 123543

  # running loop
  for(i in seq_along(probs)){

    # pull out probabilities
    prob <- rbind(probs[[i]])

    # determine response
    set.seed(seed)
    u    <- runif(nrow(prob))
    resp <- (u >= prob[ , 1]) + 1

    # compare
    set.seed(seed)
    expect_equal(simulate_mupp_resp1(prob),
                 resp)
  } # END for i LOOP

  # multiple people (four dimensions)

  # probs (resp is after generating u from seed)
  probs <- list(rbind(c(.25, .25, .25, .25),
                      c(.2, .3, .3, .2),
                      c(.1, .9, 0, 0)),
                rbind(c(.2, .2, .3, .3),
                      c(0, 0, .1, .9)))
  resp  <- list(c(3, 4, 2),
                c(3, 4))

  # seed
  seed  <- 932434

  # running loop
  for(i in seq_along(probs)){

    # pull out probabilities
    prob <- rbind(probs[[i]])

    # compare
    set.seed(seed)
    expect_equal(simulate_mupp_resp1(prob),
                 resp[[i]])
  } # END for i LOOP
})


test_that("mupp response simulation works", {

  # ONE PERSON #

  # ALL RESPONSES THE SAME
  param_seed <- 73234
  resp_seed  <- 34234

  # - one item/two dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 1,
                                  n_items   = 1,
                                  n_dims    = 2)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = 1,
                        item   = 1,
                        resp   = 2)

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)

  # - one item/three dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 1,
                                  n_items   = 1,
                                  n_dims    = 3)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = 1,
                        item   = 1,
                        resp   = 2)

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)

  # - two items/two dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 1,
                                  n_items   = 2,
                                  n_dims    = 2)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = 1,
                        item   = c(1, 2),
                        resp   = c(2, 1))

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)

  # - two items/three dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 1,
                                  n_items   = 2,
                                  n_dims    = 3)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = 1,
                        item   = c(1, 2),
                        resp   = c(2, 2))

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)


  # TWO PEOPLE #

  # ALL RESPONSES THE SAME
  param_seed  <- 72345
  resp_seed   <- 97324

  # - one item/two dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 2,
                                  n_items   = 1,
                                  n_dims    = 2)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = c(1, 2),
                        item   = 1,
                        resp   = c(1, 1))

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)

  # - one item/three dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 2,
                                  n_items   = 1,
                                  n_dims    = 3)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = c(1, 2),
                        item   = 1,
                        resp   = c(5, 5))

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)

  # - two items/two dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 2,
                                  n_items   = 2,
                                  n_dims    = 2)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = c(1, 1, 2, 2),
                        item   = c(1, 2, 1, 2),
                        resp   = c(2, 2, 2, 1))

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)

  # - two items/three dimensions

  # generate person parameters
  set.seed(param_seed)
  persons <- simulate_mupp_params(n_persons = 2,
                                  n_items   = 2,
                                  n_dims    = 3)

  # specify response (note: set.seed(resp_seed); simulate_mupp_params(); runif(n))
  resp    <- data.frame(person = c(1, 1, 2, 2),
                        item   = c(1, 2, 1, 2),
                        resp   = c(4, 2, 3, 2))

  # set seed (so response matches) and compare
  set.seed(resp_seed)
  expect_equal(do.call(what = simulate_mupp_resp,
                       args = persons)$resp,
               resp)
})
