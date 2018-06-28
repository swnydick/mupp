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
params_sim <- mupp::simulate_mupp_params(n_persons = 10000,
                                         n_items   = 50,
                                         n_dims    = 2)
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

# LOOKING AT MODEL
p_ggum <- function(theta, alpha, delta, tau){
  e1 <- exp(alpha * ((theta - delta) - tau))
  e2 <- exp(alpha * (2 * (theta - delta) - tau))
  e3 <- exp(alpha * (3 * (theta - delta)))
  (e1 + e2) / (1 + e1 + e2 + e3)
} # END p_ggum FUNCTION

# defaults
theta1 <- 0
alpha1 <- 1
tau1   <- 0
tau2   <- -1
tau3   <- 1
delta1 <- 0
delta2 <- -1
delta3 <- 1

# alpha #
alpha  <- seq(0, 3, by = .01)

prob1  <- p_ggum(theta1, alpha, delta1, tau2)
prob2  <- p_ggum(theta1, alpha, delta2, tau2)

# prob is higher (if theta = delta) for higher alpha
plot(alpha, prob1, type = "l")

# prob peaks around .75 (if theta != delta) and decreases for lower alpha
plot(alpha, prob2, type = "l")


# delta #
delta <- seq(-3, 3, by = .01)

prob1  <- p_ggum(theta1, alpha1, delta, tau1)
prob2  <- p_ggum(theta1, alpha1, delta, tau2)

# prob is higher as delta --> theta (from both sides, all being equal)
plot(delta, prob1, type = "l")

# same as prob1 but prob is higher throughout the range of delta???
plot(delta, prob2, type = "l")


# theta #
theta <- seq(-3, 3, by = .01)

prob1 <- p_ggum(theta, alpha1, delta1, tau1)
prob2 <- p_ggum(theta, alpha1, delta1, tau2)
prob3 <- p_ggum(theta, alpha1, delta1, tau3)

# 1) prob always increases as theta --> delta
plot(theta, prob1, type = "l")

# 2) if tau < delta, then prob is higher throughout the range
plot(theta, prob2, type = "l")

# 3) if tau > delta, then prob is lower throughout the range
plot(theta, prob3, type = "l")


# tau #
theta <- 0
alpha <- 1
tau   <- seq(-3, 3, by = .01)


# TESTING MCMC ALGORITHM
n_persons <- 5000
n_items   <- 60
all_par             <- simulate_mupp_params(n_persons = n_persons, n_items = n_items, n_dims = 2)
all_par$items$delta <- r_delta_prior(n = n_items, min = -2, max = 2)
all_par$items$tau   <- r_tau_prior(n = n_items, min = -2, max = 1)
all_resp            <- do.call(simulate_mupp_resp, all_par)

resp  <- all_resp$resp
items <- all_resp$items

out   <- estimate_mupp_params(resp       = resp,
                              items      = items,
                              n_iters    = 10000,
                              n_burnin   = 2000,
                              delta_sign = sign(all_par$items$delta))


# thetas
thetas_act <- dcast(all_par$persons, person ~ dim, value.var = "theta")
cor(out$mean$thetas, thetas_act[ , -1])
cor(out$mean$alphas, all_par$items$alpha)
cor(out$mean$deltas, all_par$items$delta)
cor(out$mean$taus, all_par$items$tau)

# fixing picked_orders
resp_all <- as.matrix(dcast(data = resp, formula = person ~ item, value.var = "resp")[ , -1])

# comparing bad persons

bad_t2   <- which(out$mean$thetas[ , 2] > 0 & thetas_act[ , 3] <= -1)

for(pers in bad_t2){

  # estimated + true
  thetas1  <- out$mean$thetas[pers, , drop = FALSE]
  thetas2  <- thetas_act[pers, -1, drop = FALSE]

  # always using estimated parameters
  params1  <- do.call(cbind, out$mean[c("alphas", "deltas", "taus")])

  resp1    <- resp_all[pers, , drop = FALSE]

  x <- seq(-4, 4, by = .01)
  y <- sapply(x,
              FUN = function(t2){
                thetas       <- thetas1
                thetas[1, 2] <- t2

                sum(mupp:::loglik_mupp_rank_impl(thetas, params1, items, resp1))
              })

  plot(x, y,
       main = paste0("Person ", pers),
       type = "l")
  abline(v   = c(thetas2[ , 2], thetas1[ , 2]),
         col = c("blue", "red"),
         lwd = 2)
} # END for pers LOOP
