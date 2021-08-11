
set.seed(1)

# Starting conditions -----------------------------------------------------
# initial number of trees
n0 <- 100  
# initial size distribution following a power law
d0 <- poweRlaw::rplcon(n0, 1, 2)



# Define growth function --------------------------------------------------
a_g <- 1
b_g <- 0.5
c_g <- 0.05
e_g <- 0.2  # noise
growth_func <- function(d, t) {
    ( a_g * (d^b_g) * exp(-c_g * d) * exp(rnorm(1, 0, e_g)) ) * t
}
# d_test <- seq(0, 200, length.out = 100)
# plot(d_test, growth_func(d_test, 1))



# Define survival function ------------------------------------------------
a_m <- -2
b_m <- -0.5
c_m <- 0.02
mortality_func <- function(d, t) {
    m <- 1 / ( 1 + exp(-(a_m + b_m * log(d) + c_m * d)) )
    p <- exp(-m * t)
    return(p)
}
# d_test <- seq(0, 200, length.out = 100)
# plot(d_test, mortality_func(d_test, 1))



# Define recruitment function ---------------------------------------------
# constant recruitment for now
r <- 100



# Simulation --------------------------------------------------------------
t_max  <- 200   # years to simulation
t_step <- 1
time   <- seq(0, t_max, t_step)
d_t    <- vector("list", length = length(time))
d_t[[1]] <- d0

# begin simulation
for (t in seq_len(length(time))) {
    # recruitment all at 1cm diameter
    d_r <- rep(1, r)
    # append recruits to the rest of the population
    d <- c(d_r, d_t[[t]])
    # diameter growth
    d_t[[t + 1]] <- d * (1 + growth_func(d, t_step))
    # survival
    s <- mortality_func(d_t[[t]], t_step)
    alive <- sapply(s, function(p) rbinom(1, 1, p)) == 1
    alive[d_t[[t]] <= 0] <- FALSE  # dead if size is non-positive
    d_t[[t + 1]] <- d_t[[t + 1]][sapply(s, function(p) rbinom(1, 1, p)) == 1]
}




# Plots -------------------------------------------------------------------

# population size over time
n_t <- unlist(lapply(d_t, length))
plot(c(time, max(time + t_step)), n_t, type = "l")

# final size distribution
par(mfrow = c(1, 2))
hist(d_t[[length(d_t)]])
hist(d_t[[length(d_t)]], breaks = 100)
