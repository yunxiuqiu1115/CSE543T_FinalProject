data {
  int<lower=0> N; // number of races
  matrix[N, 4] mu; // vector of prior Gaussian means
  matrix[N, 4] sigma; // vector of prior Gaussian stds
  int<lower=0> nc[N]; // number of candidates each race
  vector<lower=0, upper=1>[4] y[N]; // vector of true votings
}

parameters {
  real alpha;
  real beta;
  real<lower=0> ols_sigma;
  matrix[N, 4] gamma; // vector of underlying true support rate on election days
  vector<lower=0>[4] p[N];
}

model {
  alpha ~ normal(0, 0.1);
  beta ~ normal(1, 0.1);
  for(i in 1:N){
    for(j in 1:nc[i]){
      gamma[i,j] ~ normal(mu[i,j], sigma[i,j]);
    }
  }
  
  for(i in 1:N){
    for(j in 1:nc[i]){
      p[i,j] ~ normal(alpha + beta*gamma[i,j], ols_sigma); 
    }
  }
  
  for(i in 1:N){
    y[i,1:nc[i]] ~ dirichlet(p[i,1:nc[i]]);
  }
}

