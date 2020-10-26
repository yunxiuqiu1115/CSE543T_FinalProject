data{
	int<lower=0> I;  // number of candidate
	int<lower=0> J;  // number of days
	int<lower=0> K;  // number of polls
	int<lower=1> n[K]; // # of two-party voters in poll k
	int<lower=0> y[K]; // Democratic voters in poll k
	int<lower=1,upper=I> candidate[K]; // candidate of poll k
	int<lower=1,upper=J> day[K]; // day of poll k
	vector[I] h; // prior for each candidate
	vector[I] tau; // prior precision for each candidate
	vector[I] v; // vote share for each candidate
}

transformed data{
	vector[I] s; // prior standard deviation for each candidate
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
		y[k] ~ binomial_logit(n[k], beta[candidate[k], day[k]]);
		
	for (i in 1:I){
		beta[i,J] ~ normal(logit(h[i]), s[i]);
		for (j in 1:(J-1)){
			beta[i,j] ~ normal(beta[i,j+1], sigma_beta);
		}
	}
	
}

// generated quantities{
//   // vector[I] f;
//   vector[I] ll;
//   for (i in 1:I){
//     // f[i] = inv_logit(beta[i,J]);
//     ll[i] = lognormal_lpdf(v[i]| beta[i,J], sigma_beta);
// 	}
// 
// }
// 

