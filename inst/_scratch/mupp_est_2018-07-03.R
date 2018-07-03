# SCRATCH CODE #

rm(list = ls())

# specifying parameters
params <- simulate_mupp_params(n_persons = 1,
                               n_items   = 100,
                               n_dims    = 5,
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

est_thetas  <- estimate_mupp_thetas(resp, items)

cor(thetas, est_thetas$thetas)
