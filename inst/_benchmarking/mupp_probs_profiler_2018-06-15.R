# PROFILING CODE #
library(mupp)
library(tidyr)
library(magrittr)

# first simulate parameters and responses
params_sim <- simulate_mupp_params(n_persons = 10000,
                                   n_items   = 50,
                                   n_dims    = 2)
resp_sim   <- do.call(simulate_mupp_resp, params_sim)

# need thetas, params, items, picked_orders
thetas <- spread(params_sim$persons,
                 key   = "dim",
                 value = "theta")[ , -1] %>%
          as.matrix()
items  <- as.matrix(params_sim$items[1:3])
params <- as.matrix(params_sim$items[4:6])

resp   <- spread(resp_sim$resp[c("person", "item", "resp")],
                 key   = "item",
                 value = "resp")[ , -1] %>%
          as.matrix()

# now profile code

# a. likelihood
mupp:::start_profiler("/Users/nydicks/Desktop/mupp_profile.out")
mupp:::loglik_mupp_rank_impl(thetas, params, items, resp)
mupp:::stop_profiler()

# b. likelihood vs probability
items_1  <- items[items[ , 1] == 1, ]
params_1 <- params[1:2, ]
resp_1   <- resp[ , 1, drop = FALSE]
microbenchmark::microbenchmark(
  mupp:::loglik_mupp_rank_impl(thetas, params, items, resp),
  mupp:::loglik_mupp_rank_impl(thetas, params_1, items_1, resp_1),
  log(mupp:::p_mupp_rank_impl(thetas, params_1, items_1[ , 3], c(resp_1)))
)

# loglikelihood adding tiny bit of overhead, but not that much ...
