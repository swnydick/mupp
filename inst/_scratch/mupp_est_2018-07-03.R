# SCRATCH CODE #

rm(list = ls())

# specifying parameters
params <- simulate_mupp_params(n_persons = 10,
                               n_items   = 300,
                               n_dims    = 7,
                               max_item_dims = 2,
                               unidim_items  = TRUE)
resp   <- do.call(simulate_mupp_resp,
                  params)
thetas <- tidyr::spread(params$persons,
                        key   = "dim",
                        value = "theta")[ , -1]

# response and items
resp        <- resp$resp
items       <- params$items

est_thetas  <- estimate_mupp_thetas(resp, items, "bfgs")

cor(thetas, est_thetas$ests)
