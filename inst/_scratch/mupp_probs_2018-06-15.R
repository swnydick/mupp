# SCRATCH CODE #

rm(list = ls())

# specifying parameters
params_1 <- c(2, -2, -1)
params_2 <- c(2, 0, 0)
params_3 <- c(2, 1, -1)
params_4 <- c(2, .1, .1)
params_5 <- c(2, -3, 3)
params_6 <- c(1, .8, 0)
params_7 <- c(1, -2, -.2)

# generating parameters/thetas based on specification
params   <- do.call(what = rbind,
                    args = lapply(ls(pattern = "params\\_",
                                     envir   = .GlobalEnv),
                                  FUN = get))
thetas   <- seq(-3, 3, length.out = 1000)
thetas   <- do.call(what = cbind,
                    args = lapply(1:nrow(params),
                                  FUN = function(x) thetas))

# calculating
out      <- mupp:::p_mupp_rank_impl(thetas, params)

# restructuring
comb     <- apply(mupp:::find_all_permutations(nrow(params), 1),
                  MARGIN   = 1,
                  FUN      = paste,
                  collapse = "-")
out      <- setNames(object = as.data.frame(out),
                     nm     = comb)
out      <- cbind(theta = thetas[ , 1], out)
out      <- reshape2::melt(out,
                           id.vars       = "theta",
                           variable.name = "combination",
                           value.name    = "probability")

# plotting
library(ggplot2)
g <- ggplot(out, aes(x        = theta,
                     y        = probability,
                     color    = combination)) +
     geom_line(size = 1) +
     theme_minimal() +
     guides(color = FALSE)
print(g)

# derivatives??
library(numDeriv)

p_mupp_rank <- function(thetas, params, rank_index = 1){

  mupp:::p_mupp_rank_impl(rbind(thetas), rbind(params))[ , rank_index]

}

## TESTING LIKELIHOOD FCT ... ADD TO TESTTHAT?
library(tidyr)
library(magrittr)
params_sim <- mupp::simulate_mupp_params(n_persons = 1000,
                                         n_items   = 50,
                                         n_dims    = 6)
resp_sim   <- do.call(mupp::simulate_mupp_resp, params_sim)

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

microbenchmark::microbenchmark({
lik1   <- mupp:::loglik_mupp_rank_impl(thetas, params, items, resp)
}, {
split_vec    <- items[ , 1]
params_split <- split.data.frame(params, split_vec)
dims_split   <- split(items[ , 3], split_vec)
resp_split   <- as.list(as.data.frame(resp))

probs        <- Map(f      = mupp:::p_mupp_rank_impl,
                    thetas = list(thetas),
                    params = params_split,
                    dims   = dims_split,
                    picked_order_id  = resp_split) %>%
                do.call(what = cbind)
loglik      <- log(probs)
})

all.equal(lik1, loglik)
