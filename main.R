library(rstan)

# loading data
data = read.csv("results/forecast1992-2014.csv")
cycles = unique(data$cycle)
states = unique(data$state)

# define variables
metadata = list()
mu = list()
sigma = list()
y = list()
nc = c()
counter = 0

# iterate over races
for (cycle in cycles) {
  for (state in states) {
    pmu = data[data$cycle==cycle & data$state==state,c("posteriormean")]
    pstd = data[data$cycle==cycle & data$state==state,c("posteriorstd")]
    vote = data[data$cycle==cycle & data$state==state,c("vote")]
    if(length(pmu)){
      counter = counter + 1
      metadata[[counter]] = c(cycle, state)
      mu[[counter]] = pmu
      sigma[[counter]] = pstd
      y[[counter]] = vote / sum(vote)
      nc = c(nc, length(vote))
    }
  }
}

# build stan data
stan_mu = matrix(0,counter,4)
stan_sigma = matrix(0.0001,counter,4)
stan_y = matrix(0.0001,counter,4)

for (i in 1:counter) {
  stan_mu[i,1:nc[i]] = mu[[i]]
  stan_sigma[i,1:nc[i]] = sigma[[i]]
  stan_y[i,1:nc[i]] = y[[i]]
  stan_y[i,] = stan_y[i,]/sum(stan_y[i,])
}

# define stan data structure
stan_data = list(N = counter, 
                 mu = stan_mu, 
                 sigma = stan_sigma,
                 nc = nc,
                 y = stan_y)

# define stan model
model <- stan_model("model.stan")

# train stan model
fit <- stan(file = "model.stan",
            data = stan_data, 
            warmup = 500, 
            iter = 3000, 
            chains = 3, 
            cores = 3, 
            thin = 4,
            control=list(adapt_delta=.95)
            )

# summary(fit)

# test data
test_mu = matrix(0,1,4)
test_sigma = matrix(0,1,4)
test_mu[1,1:2] = c(0.7,0.3)
test_sigma[1,1:2] = c(0.05,0.05)
test_data = list(N = 1, 
                 mu = test_mu, 
                 sigma = test_sigma,
                 nc = c(2))
pred <- posterior_predict(model, newdata = test_data)
