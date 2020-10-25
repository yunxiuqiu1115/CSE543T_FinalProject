data{
	int<lower=0> I;  // number of states
	int<lower=0> J;  // number of days
	int<lower=0> K;  // number of polls
	int<lower=1> n[K]; // # of two-party voters in poll k
	int<lower=0> y[K]; // Democratic voters in poll k
	int<lower=1,upper=I> state[K]; // state of poll k
	int<lower=1,upper=J> day[K]; // day of poll k
	vector[I] h; // prior for each state
	vector[I] tau; // prior precision for each state
	vector[I] v; // vote share for each state
}

transformed data{
	vector[I] s; // prior standard deviation for each state
	for (i in 1:I){
	   s[i] = 1/sqrt(tau[i]);
	}
}

parameters{
	matrix[I,J] beta;
	real<lower=0> sigma_beta;
}

model{	
	sigma_beta ~ normal(0, 1);
	for (k in 1:K)
		y[k] ~ binomial_logit(n[k], beta[state[k], day[k]]);
		
	for (i in 1:I){
		beta[i,J] ~ normal(logit(h[i]), s[i]);
		for (j in 1:(J-1)){
			beta[i,j] ~ normal(beta[i,j+1], sigma_beta);
		}
	}
	
}

generated quantities{
  vector[I] ll;
  for (i in 1:I){
    ll[i] = lognormal_lpdf(v[i]| beta[i,J], sigma_beta);
	}

}


