data {
  int<lower=0> N; // number of races
  matrix[N, 4] mu; // vector of prior Gaussian means
  matrix[N, 4] sigma; // vector of prior Gaussian stds
  int<lower=0> nc[N]; // number of candidates each race
  vector<lower=0, upper=1>[4] y[N]; // vector of true votings
  int<lower=0> test_N; // number of races
  matrix[test_N, 4] test_mu; // vector of prior Gaussian means
  matrix[test_N, 4] test_sigma; // vector of prior Gaussian stds
  int<lower=0> test_nc[test_N]; // number of candidates each race
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

generated quantities {
  vector[4] test_y[test_N]; // vector of true votings
  vector[4] test_p[test_N]; // vector of underlying dirichlet parameters
  vector[4] test_gamma[test_N]; // vector of underlying true support rates
  for(i in 1:test_N){
    for(j in 1:4){
      test_gamma[i,j] = normal_rng(test_mu[i,j], test_sigma[i,j]);
      if(test_gamma[i,j]<0){
        test_gamma[i,j] = 0;
      }
      else if(test_gamma[i,j]>1){
        test_gamma[i,j]=1;
      }
    }
  }
  for(i in 1:test_N){
    for(j in 1:4){
      if(j<=test_nc[i]){
        // assume a simple linear relation
        test_p[i,j] = alpha + beta*test_gamma[i,j];
        if(test_p[i,j]<0){
          test_p[i,j]=0;
        }
      }
      else{
        // hardcode unocurring candidates
        test_p[i,j]=0.0001;
      }
    }
  }
  for (i in 1:test_N)
    test_y[i] = dirichlet_rng(test_p[i]);
}

