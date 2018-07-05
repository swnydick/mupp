# SCRATCH CODE #

rm(list = ls())

# specifying parameters
params <- simulate_mupp_params(n_persons = 10,
                               n_items   = 300,
                               n_dims    = 9,
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

diag(cor(thetas, est_thetas$thetas))

# TESTING AGAINST MCMC
act_thetas  <- tidyr::spread(all_params$persons,
                             key   = "dim",
                             value = "theta")[ , -1]
mcmc_thetas <- out$mean$thetas

bfgs_thetas   <- estimate_mupp_thetas(resp   = all_resp$resp,
                                      items  = all_params$items,
                                      method = "bfgs")
bfgs_thetas_2 <- estimate_mupp_thetas(resp   = all_resp$resp,
                                      items  = within(all_params$items, {
                                        alpha = out$mean$alphas
                                        delta = out$mean$deltas
                                        tau   = out$mean$taus
                                      }))
