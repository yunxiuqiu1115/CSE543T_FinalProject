data {
  int n_candidate; //  number of candidates
  int N_max; // max number of polls
  int N[n_candidate]; // number of polls for each race;
  int days[n_candidate,N_max]; // days until election
  vector[N_max] ns[n_candidate]; // poll samplesizes
  vector[N_max] ys[n_candidate]; // polling proportions
  vector<lower=0,upper=1>[n_candidate] h; // prior estimates
  int J[n_candidate]; // max value of days until election
  int J_max;
  
}

parameters {
  vector[J_max] betas[n_candidate];
  real<lower=0> sigma_J;
  real<lower=0> sigma_beta;
  
}

transformed parameters{
  vector<lower=0,upper=1>[N_max] ps[n_candidate];
  vector<lower=0>[N_max] vs[n_candidate];
  for(i in 1:n_candidate){
    for(j in 1:N[i]){
      ps[i,j] = inv_logit(betas[i, days[i,j]]);
      vs[i,j] = sqrt(ps[i,j]*(1-ps[i,j])/ns[i,j]);
    }
  }
  
}

model {
  sigma_beta ~ normal(0,10);
  sigma_J ~ normal(0,10);
  for(i in 1:n_candidate){
     betas[i,1] ~ normal(logit(h[i]), sigma_J);
     for(j in 2:J[i]){
       betas[i,j] ~ normal(betas[i,j-1], sigma_beta);
      }
      for(j in 1:J[i]){
        ys[i,j] ~ normal(ps[i,j], vs[i,j]);
      }
      
  }
  
}

generated quantities{
  vector<lower=0,upper=1>[n_candidate] y;
  for(i in 1:n_candidate){
      y[i] = inv_logit(betas[i,1]);
  }
  
}

