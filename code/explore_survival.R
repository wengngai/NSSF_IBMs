# a null survival model
p     <- 0.5        # per-time survival probability
t_d   <- 1/365.25   # daily time step
t_y   <- 1          # annual time step
t_max <- 22         # maximum number of year to simulate
n0    <- 1e3        # initial population size

# simulation
nsim <- 1000         # number of simulation before taking the mean

# at daily time step
time_d <- seq(0, t_max, t_d)
n_d    <- matrix(NA, nsim, length(time_d) + 1)
n_d[, 1] <- n0
for (i in seq_len(nsim)) {
    for (j in seq_len(length(time_d))) {
        s <- p ^ t_d
        n_d[i, j + 1] <- rbinom(1, n_d[i, j], s)
    }
}
n_d <- colMeans(n_d)

# at annual time step
time_y <- seq(0, t_max, t_y)
n_y    <- matrix(NA, nsim, length(time_y) + 1)
n_y[, 1] <- n0
for (i in seq_len(nsim)) {
    for (j in seq_len(length(time_y))) {
        s <- p ^ t_y
        n_y[i, j + 1] <- rbinom(1, n_y[i, j], s)
    }
}
n_y <- colMeans(n_y)

# plot
plot(c(time_d, time_d[length(time_d)] + t_d), n_d, ylab = "N", xlab = "t")
points(c(time_y, time_y[length(time_y)] + t_y), n_y, col = "red", pch = 19)
