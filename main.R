library(rstan)

data = read.csv("results/forecast1992-2014.csv")
cycles = unique(data$cycle)
states = unique(data$state)

metadata = list()
mu = list()
sigma = list()
y = list()
nc = c()
counter = 0

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

stan_mu = matrix(0,counter,4)
stan_sigma = matrix(0,counter,4)
stan_y = matrix(0,counter,4)

for (i in 1:counter) {
  stan_mu[i,1:nc[i]] = mu[[i]]
  stan_sigma[i,1:nc[i]] = sigma[[i]]
  stan_y[i,1:nc[i]] = y[[i]]
}

stan_data = list(N = counter, 
                 mu = stan_mu, 
                 sigma = stan_sigma,
                 nc = nc,
                 y = stan_y)

fit <- stan(file = "model.stan",
            data = stan_data, 
            warmup = 500, 
            iter = 10000, 
            chains = 4, 
            cores = 4, 
            thin = 4)
