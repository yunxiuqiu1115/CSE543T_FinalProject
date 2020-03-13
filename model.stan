data {
  int<lower=0> N; // number of races
  matrix[N, 4] mu; // vector of prior Gaussian means
  matrix[N, 4] sigma; // vector of prior Gaussian stds
  int<lower=0> nc[N]; // number of candidates each race
  vector<lower=0, upper=1>[4] y[N]; // vector of true votings
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  vector<lower=0, upper=1>[4] gamma[N]; // vector of underlying true support rates
}

transformed parameters{
  vector<lower=0>[4] p[N]; // vector of underlying dirichlet parameters
  for(i in 1:N){
    for(j in 1:4){
      if(j<=nc[i]){
        // assume a simple linear relation
        p[i,j] = alpha + beta*gamma[i,j];
      }
      else{
        // hardcode unocurring candidates
        p[i,j]=0.0001;
      }
    }
  }
}

model {
  // very flat priors on linear model parameters
  alpha ~ normal(0, 1);
  beta ~ normal(1, 1);
  // generate underlying true support rates
  for(i in 1:N){
    for(j in 1:4){
      gamma[i,j] ~ normal(mu[i,j], sigma[i,j]);
    }
  }

  // generate voting rates from support rates
  for(i in 1:N){
    y[i] ~ dirichlet(p[i]);
  }
}

