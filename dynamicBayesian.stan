data {
  int N;
  int days[N];
  vector[N] ns;
  vector[N] ys;
  real h;
  real v;
  int J;
  
}

parameters {
  real betas[J];
  real<lower=0> sigma_J;
  real<lower=0> sigma_beta;
  
}

transformed parameters{
  real ps[N];
  real vs[N];
  for(i in 1:N){
    ps[i] = inv_logit(betas[days[i]]);
    vs[i] = sqrt(ps[i]*(1-ps[i])/ns[i]);
  }
  
}

model {
  sigma_beta ~ normal(0,10);
  sigma_J ~ normal(0,10);
  betas[1] ~ normal(logit(h), sigma_J);
  for(i in 2:J){
    betas[i] ~ normal(betas[i-1], sigma_beta);
  }
  ys ~ normal(ps, vs);
  
}

generated quantities{
  real y;
  real ll;
  y = inv_logit(betas[1]);
  ll = lognormal_lpdf(v|betas[1], sigma_beta);
  
}

